"""
Fit eta restringido al regimen Rydberg X<100, reutilizando exactamente
el modelo de calibracion del repo:  err_lin = A0 + Gamma*X^2 + eta*(d*X)

Se replica run_calibration.calibrate_and_run pero (a) sin winsorizar,
(b) recortando el train set a X < X_CUT en vez de cap. Devuelve eta_hat
y banda bootstrap del 95%.
"""
import sys, os
import numpy as np
import pandas as pd

SRC = os.path.join(os.path.dirname(__file__), "..", "..", "..", "src")
sys.path.insert(0, os.path.abspath(SRC))

from doftstudy.run_calibration import (
    f_from_K, f_from_eV, find_best_lock, generate_prime_products,
    read_input_csv, N_BOOTSTRAPS,
)
from scipy.optimize import lsq_linear
from sklearn.utils import resample

INPUT = os.path.join(os.path.dirname(__file__), "..", "..", "raw",
                     "materials_clusters_real_v6.csv")
FIT_CATEGORIES = ["SC_TypeI", "SC_TypeII"]   # default del pipeline
PRIME_MAX = 10000
X_CUT = 100.0
kB = 1.380649e-23
h = 6.62607015e-34


def build_records(df):
    generate_prime_products(max_val=PRIME_MAX)
    df_single = df[df['sub_network'].str.lower() == 'single'].copy()
    recs = []
    for _, row in df_single.iterrows():
        name = str(row['name']); cat = str(row['category'])
        lock_family = str(row['lock_family']).lower()
        f_th = f_from_K(row['Tc_K']); f_debye = f_from_K(row['ThetaD_K'])
        f_elec = f_from_eV(row['EF_eV'])
        f_gap = np.nan
        if pd.notnull(row['Gap_meV']) and float(row['Gap_meV']) > 0:
            f_gap = f_from_eV(row['Gap_meV'] * 1e-3)
        elif pd.notnull(row['Tc_K']) and float(row['Tc_K']) > 0:
            f_gap = (1.76 * kB * float(row['Tc_K'])) / h
        anchors = [("thermal@Tc", f_th), ("gap", f_gap), ("Debye", f_debye), ("EF", f_elec)]
        jump_defs = [(0, 1, "t->gap", 3), (1, 2, "gap->Debye", 2), (2, 3, "Debye->EF", 1)]
        X = (float(row['ThetaD_K']) / float(row['Tc_K'])
             if (pd.notnull(row['Tc_K']) and float(row['Tc_K']) > 0
                 and pd.notnull(row['ThetaD_K']) and float(row['ThetaD_K']) > 0) else np.nan)
        for i, j, label, d in jump_defs:
            fa = anchors[i][1]; fb = anchors[j][1]
            if not (np.isfinite(fa) and np.isfinite(fb) and fa > 0):
                continue
            pe, pv, chosen = find_best_lock(fb / fa, lock_family, cat)
            recs.append({"name": name, "category": cat, "d": d, "X": X, "err_before": pe})
    return pd.DataFrame(recs)


def r2(Xv, dv, yv, coef):
    A0, G, E = coef
    yhat = A0 + G * Xv**2 + E * (dv * Xv)
    ss_res = np.sum((yv - yhat) ** 2)
    ss_tot = np.sum((yv - yv.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def fit(Xv, dv, yv):
    A = np.vstack([np.ones_like(Xv), Xv**2, dv * Xv]).T
    bounds = ([-np.inf, 0.0, 0.0], [np.inf, np.inf, np.inf])
    etas, gammas = [], []
    idx = np.arange(len(Xv))
    rng = np.random.RandomState(0)
    for _ in range(N_BOOTSTRAPS):
        b = resample(idx, random_state=rng.randint(1 << 30))
        r = lsq_linear(np.vstack([np.ones_like(Xv[b]), Xv[b]**2, dv[b] * Xv[b]]).T,
                       yv[b], bounds=bounds)
        if r.success:
            gammas.append(r.x[1]); etas.append(r.x[2])
    point = lsq_linear(A, yv, bounds=bounds).x
    return point, np.array(etas), np.array(gammas)


def main():
    df = read_input_csv(os.path.abspath(INPUT))
    out = build_records(df)
    mask = (out['category'].isin(FIT_CATEGORIES)
            & np.isfinite(out['X']) & np.isfinite(out['err_before']) & np.isfinite(out['d']))
    train = out[mask]
    nmat = train['name'].nunique()
    print(f"Train jumps totales (SC_TypeI/II, single): {len(train)}  | materiales: {nmat}")
    HERE = os.path.dirname(__file__)
    summary_rows = []  # se vuelca a CSV al final

    # --- Exportar el per-salto canonico ---
    cols = ['name', 'category', 'd', 'X', 'err_before']
    expo = train[cols].copy()
    expo['dX'] = expo['d'] * expo['X']
    expo['X2'] = expo['X'] ** 2
    out_csv = os.path.join(os.path.dirname(__file__), "calib_perjump_SCTypeI_II.csv")
    expo.sort_values(['name', 'd']).to_csv(out_csv, index=False)
    print(f"Per-salto canonico exportado -> {out_csv}")

    # --- Fit canonico winsorizado (cap) con R^2 en su propio set ---
    for cap in [400.0, 600.0, 800.0]:
        Xv = np.minimum(train['X'].values, cap)
        dv = train['d'].values.astype(float); yv = train['err_before'].values
        point, etas, gammas = fit(Xv, dv, yv)
        R2_point = r2(Xv, dv, yv, point)
        R2_boot = r2(Xv, dv, yv, [point[0], gammas.mean(), etas.mean()])
        lo, hi = np.percentile(etas, [2.5, 97.5])
        print(f"[canonico cap={cap:.0f}] eta_boot={etas.mean():.4e}  "
              f"R2_point={R2_point:.4f}  R2_boot={R2_boot:.4f}  N={len(Xv)}")
        summary_rows.append({
            "config": f"canonico_winsor_cap{cap:.0f}", "winsor_cap": cap, "X_cut": "",
            "N_saltos": len(Xv), "N_materiales": train['name'].nunique(),
            "A0": point[0], "gamma_point": point[1], "eta_point": point[2],
            "eta_boot_mean": etas.mean(), "eta_ci_low": lo, "eta_ci_high": hi,
            "R2_point": R2_point})

    for cut, lbl in [(np.inf, "TODO (sin corte)"), (X_CUT, "X<100 (Rydberg)")]:
        sub = train[train['X'] < cut]
        Xv = sub['X'].values; dv = sub['d'].values.astype(float); yv = sub['err_before'].values
        point, etas, gammas = fit(Xv, dv, yv)
        lo, hi = np.percentile(etas, [2.5, 97.5])
        R2_point = r2(Xv, dv, yv, point)
        print(f"\n=== {lbl} ===")
        print(f"  N saltos={len(sub)}  materiales={sub['name'].nunique()}  Xmax={Xv.max():.1f}")
        print(f"  eta_point = {point[2]:.4e}   gamma_point={point[1]:.4e}  A0={point[0]:.4e}")
        print(f"  eta_boot  = {etas.mean():.4e}  (95% CI [{lo:.4e}, {hi:.4e}])")
        print(f"  R2_point  = {R2_point:.4f}")
        summary_rows.append({
            "config": f"X_lt_{cut}" if np.isfinite(cut) else "todo_sin_corte",
            "winsor_cap": "", "X_cut": cut if np.isfinite(cut) else "",
            "N_saltos": len(sub), "N_materiales": sub['name'].nunique(),
            "A0": point[0], "gamma_point": point[1], "eta_point": point[2],
            "eta_boot_mean": etas.mean(), "eta_ci_low": lo, "eta_ci_high": hi,
            "R2_point": R2_point})

    # --- LOMO: leave-one-material-out sobre el fit X<100 ---
    sub = train[train['X'] < X_CUT]
    mats = sorted(sub['name'].unique())
    base = lsq_linear(np.vstack([np.ones(len(sub)), sub['X'].values**2,
                                 sub['d'].values.astype(float) * sub['X'].values]).T,
                      sub['err_before'].values, bounds=([-np.inf, 0, 0], [np.inf, np.inf, np.inf])).x
    R2_full = r2(sub['X'].values, sub['d'].values.astype(float), sub['err_before'].values, base)
    print(f"\n=== LOMO (X<100)  R2_full={R2_full:.4f}  materiales={len(mats)} ===")
    rows = []
    for m in mats:
        s = sub[sub['name'] != m]
        Xv = s['X'].values; dv = s['d'].values.astype(float); yv = s['err_before'].values
        coef = lsq_linear(np.vstack([np.ones_like(Xv), Xv**2, dv * Xv]).T, yv,
                          bounds=([-np.inf, 0, 0], [np.inf, np.inf, np.inf])).x
        R2_m = r2(Xv, dv, yv, coef)
        n_dropped = (sub['name'] == m).sum()
        rows.append((m, n_dropped, R2_m, coef[2]))
    rows.sort(key=lambda r: r[2])  # peor R2 primero
    print(f"  {'material_removido':<22} {'n_saltos':>8} {'R2_LOMO':>9} {'eta_LOMO':>11}  dR2")
    for m, n, R2_m, eta_m in rows:
        print(f"  {m:<22} {n:>8} {R2_m:>9.4f} {eta_m:>11.3e}  {R2_m - R2_full:+.4f}")
    r2s = np.array([r[2] for r in rows])
    print(f"  --- R2_LOMO: min={r2s.min():.4f}  max={r2s.max():.4f}  "
          f"mediana={np.median(r2s):.4f}  (full={R2_full:.4f})")

    # --- Volcado a CSV ---
    lomo_df = pd.DataFrame(rows, columns=["material_removido", "n_saltos_removidos",
                                          "R2_LOMO", "eta_LOMO"])
    lomo_df["R2_full"] = R2_full
    lomo_df["dR2"] = lomo_df["R2_LOMO"] - R2_full
    lomo_df["eta_full"] = base[2]
    lomo_path = os.path.join(HERE, "lomo_rydberg_X100.csv")
    lomo_df.to_csv(lomo_path, index=False)
    print(f"\nLOMO -> {lomo_path}")

    summary_path = os.path.join(HERE, "fit_summary_eta.csv")
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"Resumen fits -> {summary_path}")


if __name__ == "__main__":
    main()
