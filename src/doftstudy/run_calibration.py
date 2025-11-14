import os, math, argparse, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Importar dependencias de scikit-learn y scipy
try:
    from scipy.optimize import lsq_linear
    from sklearn.preprocessing import StandardScaler
    from sklearn.utils import resample
    from numpy.linalg import cond
except ImportError:
    print("ERROR: Faltan dependencias. Ejecuta: pip install scipy scikit-learn")
    exit()

# --- Constantes físicas ---
kB = 1.380649e-23; h = 6.62607015e-34; eV = 1.602176634e-19
N_BOOTSTRAPS = 500
CONFIG_FILE = "doft_config.json"
DEFAULT_X_CAP = 600.0 # Límite para Winsorizing (tu regla de He3)

# --- Funciones de conversión ---
def f_from_K(T_K): 
    try: T = float(T_K); return (kB*T)/h if T>0 else np.nan
    except Exception: return np.nan

def f_from_eV(E_eV): 
    try: E = float(E_eV); return (E*eV)/h if E>0 else np.nan
    except Exception: return np.nan

# --- Funciones de Locking ---
def best_integer_lock(R, max_exp=9):
    if not np.isfinite(R) or R<=0: return (np.nan, np.nan, "integer")
    best_err, best_val = None, None
    for a in range(0, max_exp+1):
        val2 = (2**a);
        if val2>R*10: break
        for b in range(0, max_exp+1):
            val23 = val2 * (3**b);
            if val23>R*10: break
            for c in range(0, max_exp+1):
                val235 = val23 * (5**c);
                if val235>R*10: break
                for d in range(0, max_exp+1):
                    val = val235 * (7**d);
                    if val>R*10: break
                    err = abs(val-R)/R
                    if (best_err is None) or (err<best_err): best_err, best_val = err, val
    return (best_err, best_val, "integer")

SMALL_DENOMINATORS = np.array([1, 2, 3, 4, 5, 6, 7, 8]) # Tu regla Q<=8
PRIME_VALS = None # Se generará una sola vez

def generate_prime_products(max_val=10000):
    global PRIME_VALS
    # --- ¡CAMBIO! No regenerar si ya existe con la misma valencia ---
    # Esto es útil si el wrapper lo llama múltiples veces.
    if PRIME_VALS is not None and len(PRIME_VALS) > 0 and PRIME_VALS.max() >= max_val:
        # print("--- Usando lista 'prime values' existente. ---")
        return
        
    print(f"--- Generando lista de 'prime values' (P, Q) hasta {max_val} ---")
    vals = set()
    for a in range(0, int(math.log2(max_val)) + 2):
        val2 = 2**a;
        if val2 > max_val: break
        for b in range(0, int(math.log(max_val, 3)) + 2):
            val23 = val2 * (3**b);
            if val23 > max_val: break
            for c in range(0, int(math.log(max_val, 5)) + 2):
                val235 = val23 * (5**c);
                if val235 > max_val: break
                for d in range(0, int(math.log(max_val, 7)) + 2):
                    val = val235 * (7**d);
                    if val > max_val: break
                    if val > 0: vals.add(val)
    PRIME_VALS = np.array(sorted(list(vals)))
    print(f"--- Lista generada con {len(PRIME_VALS)} valores únicos. ---")

def best_rational_lock(R):
    if not np.isfinite(R) or R <= 0: return (np.nan, np.nan, "rational")
    best_err, best_val = None, None
    for Q in SMALL_DENOMINATORS:
        target_P = R * Q; idx = np.searchsorted(PRIME_VALS, target_P)
        indices_to_check = [idx, idx - 1]
        for i in indices_to_check:
            if 0 <= i < len(PRIME_VALS):
                P = PRIME_VALS[i]
                if Q == 0 or P == 0: continue
                val = P / Q; err = abs(val - R) / R
                if (best_err is None) or (err < best_err): best_err, best_val = err, val
    return (best_err, best_val, "rational")

def find_best_lock(R, lock_family, category):
    """
    Implementa tu regla 'lock_family' del CSV Y los 'Asserts de coherencia'
    """
    if not np.isfinite(R):
        return (np.nan, np.nan, lock_family)
    
    # --- Assert de coherencia DOFT (automático en el pipeline) ---
    if category in ['SC_TypeI', 'SC_TypeII'] and lock_family == 'integer':
        # Type I/II: prohibir p/q en la elección final.
        return best_integer_lock(R)
        
    if category in ['Superfluid', 'SC_Molecular'] and lock_family == 'rational':
        # Superfluid/Molecular: prohibir entero en la elección final.
        return best_rational_lock(R)
    
    # --- Lógica estándar para 'mixed' o casos no forzados ---
    if lock_family == 'integer':
        return best_integer_lock(R)
        
    if lock_family == 'rational':
        return best_rational_lock(R)
        
    if lock_family == 'mixed':
        # Probar p/q primero si R < 1, como pediste
        if R < 1.0:
            err_r, val_r, _ = best_rational_lock(R)
            err_i, val_i, _ = best_integer_lock(R)
        else:
            err_i, val_i, _ = best_integer_lock(R)
            err_r, val_r, _ = best_rational_lock(R)
            
        if np.isnan(err_i) and np.isnan(err_r): return (np.nan, np.nan, "mixed")
        
        # Manejar el caso donde uno es nan
        if np.isnan(err_r): return (err_i, val_i, "integer")
        if np.isnan(err_i): return (err_r, val_r, "rational")
        
        # Comparar si ambos son válidos
        if (err_i <= err_r):
            return (err_i, val_i, "integer")
        else:
            return (err_r, val_r, "rational")
            
    return (np.nan, np.nan, lock_family) # Default


# --- ¡CAMBIO! Se añadió 'prime_max_val' como argumento ---
def calibrate_and_run(df, outdir, run_label, fit_categories, winsor_X_cap, prime_max_val=10000):
    records = []
    
    # --- ¡CAMBIO! Se pasa el argumento a la función ---
    generate_prime_products(max_val=prime_max_val) # Genera la lista global PRIME_VALS
    
    print("\n--- 1. Pre-procesando todos los materiales 'single' ---")
    df_single = df[df['sub_network'].str.lower() == 'single'].copy()
    
    for _, row in df_single.iterrows():
        try:
            name = str(row['name'])
            cat = str(row['category'])
            lock_family = str(row['lock_family']).lower()
            
            f_th = f_from_K(row['Tc_K'])
            f_debye = f_from_K(row['ThetaD_K'])
            f_elec = f_from_eV(row['EF_eV']) # Null-safe
            
            f_gap = np.nan
            if pd.notnull(row['Gap_meV']) and float(row['Gap_meV']) > 0:
                f_gap = f_from_eV(row['Gap_meV'] * 1e-3)
            elif pd.notnull(row['Tc_K']) and float(row['Tc_K']) > 0:
                f_gap = (1.76 * kB * float(row['Tc_K'])) / h
            
            anchors = [
                ("thermal@Tc", f_th), ("gap(Δ)", f_gap),
                ("Debye", f_debye), ("electronic(E_F)", f_elec),
            ]
            jump_defs = [
                (0,1,"thermal@Tc→gap(Δ)", 3), (1,2,"gap(Δ)→Debye", 2),
                (2,3,"Debye→electronic(E_F)", 1),
            ]
            
            X = float(row['ThetaD_K']) / float(row['Tc_K']) if (pd.notnull(row['Tc_K']) and float(row['Tc_K']) > 0 and pd.notnull(row['ThetaD_K']) and float(row['ThetaD_K']) > 0) else np.nan

            for i,j,label,d in jump_defs:
                fa = anchors[i][1]; fb = anchors[j][1]
                
                if not (np.isfinite(fa) and np.isfinite(fb) and fa > 0):
                    records.append({
                        "name": name, "category": cat, "sub_network": 'single',
                        "jump_desc": label, "d": d, "X": X, "lock_family": lock_family,
                        "R_obs": np.nan, "chosen_lock": "skip (nan)", "prime_value": np.nan, 
                        "err_before": np.nan, "R_corr_eta": np.nan, "err_after_eta": np.nan,
                        "R_corr_kappa": np.nan, "err_after_kappa": np.nan, "C_AB": 0.0
                    })
                    continue

                R_obs = fb/fa
                pe, pv, chosen_lock = find_best_lock(R_obs, lock_family, cat)
                
                records.append({
                    "name": name, "category": cat, "sub_network": 'single',
                    "jump_desc": label, "d": d, "X": X, "lock_family": lock_family,
                    "R_obs": R_obs, "chosen_lock": chosen_lock, "prime_value": pv, 
                    "err_before": pe, "R_corr_eta": np.nan, "err_after_eta": np.nan,
                    "R_corr_kappa": np.nan, "err_after_kappa": np.nan, "C_AB": 0.0
                })
        
        except Exception as e:
            print(f"Error procesando fila: {row['name']}. Error: {e}")

    out = pd.DataFrame(records)
    out['err_lin'] = out['err_before'] 

    print(f"\n--- 2. Calibrando (Gamma, Eta) usando solo: {fit_categories} ---")
    
    fit_mask = (out['category'].isin(fit_categories)) & (out['sub_network'] == 'single')
    train_mask = fit_mask & np.isfinite(out["X"]) & np.isfinite(out["err_lin"]) & np.isfinite(out["d"])
    
    Xv_fit = out.loc[train_mask,"X"].values
    dv_fit = out.loc[train_mask,"d"].values.astype(float)
    yv_fit = out.loc[train_mask,"err_lin"].values
    fit_names = out.loc[train_mask, "name"].unique()
    
    CALIBRATED_GAMMA = 0.0; CALIBRATED_ETA = 0.0; eta_ci_low = 0.0; eta_ci_high = 0.0
    
    if len(yv_fit) < 3:
        print(f"ADVERTENCIA: No hay suficientes datos de calibración (N={len(yv_fit)}). Se usarán Gamma=0, Eta=0.")
    else:
        print(f"--- 2a. Winsorizando X (límite superior = {winsor_X_cap}) ---")
        Xv_fit_winsor = np.minimum(Xv_fit, winsor_X_cap)
        
        print("--- 2b. Chequeo de Estabilidad del Ajuste ---")
        A_fit = np.vstack([np.ones_like(Xv_fit_winsor), Xv_fit_winsor**2, dv_fit*Xv_fit_winsor]).T
        scaler = None
        try:
            cond_num = cond(A_fit)
            print(f"Número de Condición de la Matriz: {cond_num:.2e}")
            if cond_num > 1e6:
                print("¡ADVERTENCIA! Número de Condición > 1e6. Activando auto-re-escalado (z-score).")
                scaler = StandardScaler()
                A_features_scaled = scaler.fit_transform(A_fit[:, 1:])
                A_fit = np.hstack([A_fit[:, 0:1], A_features_scaled])
        except Exception as ex:
            print(f"Error al calcular el Número de Condición: {ex}. Omitiendo escalado.")
            
        print(f"--- 2c. Calculando CIs (Bootstrap N={N_BOOTSTRAPS}) ---")
        print("Aplicando ajuste con restricciones (Gamma >= 0, Eta >= 0)...")
        
        bootstrap_gammas = []; bootstrap_etas = []
        bounds = ([-np.inf, 0.0, 0.0], [np.inf, np.inf, np.inf]) # A0 (libre), Gf (>=0), Ef (>=0)
        indices = np.arange(len(Xv_fit_winsor))

        for i in range(N_BOOTSTRAPS):
            boot_indices = resample(indices)
            Xv_boot, dv_boot, yv_boot = Xv_fit_winsor[boot_indices], dv_fit[boot_indices], yv_fit[boot_indices]
            A_boot = np.vstack([np.ones_like(Xv_boot), Xv_boot**2, dv_boot*Xv_boot]).T
            if scaler:
                A_boot_features_scaled = scaler.transform(A_boot[:, 1:])
                A_boot = np.hstack([A_boot[:, 0:1], A_boot_features_scaled])
            result = lsq_linear(A_boot, yv_boot, bounds=bounds)
            if result.success:
                A0_s, Gf_s, Ef_s = result.x
                if scaler:
                    Gf_orig = Gf_s / scaler.scale_[0] if scaler.scale_[0] != 0 else 0
                    Ef_orig = Ef_s / scaler.scale_[1] if scaler.scale_[1] != 0 else 0
                else:
                    Gf_orig, Ef_orig = Gf_s, Ef_s
                bootstrap_gammas.append(Gf_orig); bootstrap_etas.append(Ef_orig)
        
        CALIBRATED_GAMMA = np.mean(bootstrap_gammas)
        CALIBRATED_ETA = np.mean(bootstrap_etas)
        eta_ci_low = np.percentile(bootstrap_etas, 2.5)
        eta_ci_high = np.percentile(bootstrap_etas, 97.5)
        
        print("--- Resultados del Bootstrap ---")
        print(f"Gamma (g): Media={CALIBRATED_GAMMA:.2e}, StdDev={np.std(bootstrap_gammas):.2e}")
        print(f"           95% CI=[{np.percentile(bootstrap_gammas, 2.5):.2e}, {np.percentile(bootstrap_gammas, 97.5):.2e}]")
        print(f"Eta (e):   Media={CALIBRATED_ETA:.2e}, StdDev={np.std(bootstrap_etas):.2e}")
        print(f"           95% CI=[{eta_ci_low:.2e}, {eta_ci_high:.2e}]")

        print("\n--- 2d. Reporte de Influencia (Leave-One-Out) ---")
        print("Calculando cambio en Eta (e) al remover cada material de calibración:")
        influence_report = {}; eta_baseline = CALIBRATED_ETA
        for name_to_remove in fit_names:
            loo_mask = train_mask & (out["name"] != name_to_remove)
            Xv_loo = out.loc[loo_mask,"X"].values; dv_loo = out.loc[loo_mask,"d"].values.astype(float); yv_loo = out.loc[loo_mask,"err_lin"].values
            Xv_loo = np.minimum(Xv_loo, winsor_X_cap) # Aplicar Winsorizing también aquí
            if len(yv_loo) < 3:
                influence_report[name_to_remove] = float("nan"); continue
            A_loo = np.vstack([np.ones_like(Xv_loo), Xv_loo**2, dv_loo*Xv_loo]).T
            if scaler:
                A_loo_features_scaled = scaler.transform(A_loo[:, 1:])
                A_loo = np.hstack([A_loo[:, 0:1], A_loo_features_scaled])
            result_loo = lsq_linear(A_loo, yv_loo, bounds=bounds)
            if result_loo.success:
                A0_loo_s, Gf_loo_s, Ef_loo_s = result_loo.x
                if scaler:
                    Ef_loo_orig = Ef_loo_s / scaler.scale_[1] if scaler.scale_[1] != 0 else 0
                else:
                    Ef_loo_orig = Ef_loo_s
                rel_change = (Ef_loo_orig - eta_baseline) / eta_baseline if eta_baseline != 0 else 0.0
                influence_report[name_to_remove] = rel_change
            else:
                influence_report[name_to_remove] = float("nan")
        print("Material (Removido) | Cambio Relativo en Eta (e)")
        print("----------------------|---------------------------")
        sorted_influence = sorted(influence_report.items(), key=lambda item: item[1], reverse=True)
        for name, change in sorted_influence:
            print(f"{name:<21} | {change:+.2%}")

    print(f"\n--- 3. Guardando parámetros en {CONFIG_FILE} ---")
    config_data = {
        'CALIBRATED_GAMMA': CALIBRATED_GAMMA,
        'CALIBRATED_ETA': CALIBRATED_ETA,
        'METADATA': {
            'winsor_X_cap': winsor_X_cap,
            'eta_95_CI': [eta_ci_low, eta_ci_high],
            'fit_families': fit_categories,
            'prime_max_val': prime_max_val # ¡CAMBIO! Guardar el metadato
        }
    }
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config_data, f, indent=4)

    print("\n--- 4. Aplicando corrección final a todos los materiales 'single'... ---")
    
    g = CALIBRATED_GAMMA; e = CALIBRATED_ETA
    out_single = out[out['sub_network'] == 'single'].copy()
    
    def apply_correction(row, g, e):
        if not (np.isfinite(row['R_obs']) and np.isfinite(row['X']) and np.isfinite(row['d'])):
            return np.nan
        X = row['X']; d = row['d']
        # Tu regla: R_corr = R_obs * max(0.2, 1 - η*d*X)
        correction_factor = 1 - g*(X**2) - e*(d*X)
        clamped_factor = max(0.2, correction_factor)
        return row['R_obs'] * clamped_factor

    out_single['R_corr_eta'] = out_single.apply(apply_correction, args=(g, e), axis=1)
    out_single['err_after_eta'] = np.abs(out_single['R_corr_eta'] - out_single['prime_value']) / np.abs(out_single['R_corr_eta'])
    out_single['R_corr_kappa'] = np.nan
    out_single['err_after_kappa'] = np.nan
    
    print(f"--- 5. Guardando archivos de resultados en: {outdir} ---")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    cols_order = [
        "name", "sub_network", "category", "lock_family", "jump_desc", "d", "X", "C_AB",
        "R_obs", "chosen_lock", "prime_value", "err_before", 
        "R_corr_eta", "err_after_eta", 
        "R_corr_kappa", "err_after_kappa"
    ]
    final_cols = [c for c in cols_order if c in out_single.columns]
    out_single = out_single[final_cols]
    
    out_single.to_csv(os.path.join(outdir, f"results_{run_label}.csv"), index=False, float_format='%.6g')
    
    by_sys = out_single.groupby("name")[["err_before","err_after_eta"]].sum(min_count=1).reset_index()
    by_sys.rename(columns={"err_before": "sum_err_before", "err_after_eta": "sum_err_after_full"}, inplace=True)
    by_sys.to_csv(os.path.join(outdir, f"summary_{run_label}.csv"), index=False, float_format='%.6g')
    
    print("--- ¡Calibración y análisis 'single' completos! ---")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DOFT Study - Calibrator (v-Master 1.2)")
    parser.add_argument("-i", "--input", dest="input_csv", type=str, required=True, help="Input CSV file (e.g., materials_clusters_real_v2.csv)")
    parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True, help="Output directory (e.g., 'results_calibration')")
    parser.add_argument("-l", "--label", dest="run_label", type=str, default="calibration", help="Run label for output files")
    parser.add_argument("-f", "--fit_families", dest="fit_families", type=str, default="SC_TypeI,SC_TypeII", help="Comma-separated list of categories to use for fitting")
    parser.add_argument("--winsor_X", dest="winsor_X_cap", type=float, default=DEFAULT_X_CAP, help=f"Límite superior para X en el ajuste de Eta (default: {DEFAULT_X_CAP})")
    
    # --- ¡NUEVO ARGUMENTO! ---
    parser.add_argument("--prime_max_val", dest="prime_max_val", type=int, default=10000, help="Valor máximo para la generación de 'prime values' (default: 10000)")

    args = parser.parse_args()
    
    try:
        df = pd.read_csv(args.input_csv)
        fit_categories = args.fit_families.split(',')
        
        # --- ¡CAMBIO! Se pasa el nuevo argumento ---
        calibrate_and_run(df, args.outdir, args.run_label, fit_categories, args.winsor_X_cap, args.prime_max_val)
            
    except FileNotFoundError:
        print(f"ERROR: No se pudo encontrar el archivo de entrada: {args.input_csv}")
    except Exception as e:
        print(f"Ha ocurrido un error inesperado: {e}")