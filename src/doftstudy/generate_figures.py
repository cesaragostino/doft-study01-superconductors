"""Script para generar las figuras solicitadas del estudio DOFT."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")  # Generación headless

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.stats import gaussian_kde
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DEFAULT_FIGURES_DIR: Optional[Path] = None
DEFAULT_BASE_TAG = "w600_p10000"
DEFAULT_FIT_CATEGORIES = ("SC_TypeI", "SC_TypeII")

CATEGORY_LABELS = {
    "SC_TypeI": "Type-I",
    "SC_TypeII": "Type-II",
    "SC_HighPressure": "High-pressure",
    "SC_Binary": "Binary",
    "SC_Molecular": "Molecular",
    "Superfluid": "Superfluid",
    "SC_Oxide": "Oxide",
    "SC_IronBased": "Fe-based",
}

plt.style.use("seaborn-v0_8-whitegrid")


def _scale_axis(data: np.ndarray) -> Tuple[np.ndarray, str]:
    finite = data[np.isfinite(data) & (data != 0)]
    if finite.size == 0:
        return data, " (dimensionless)"
    median_abs = np.median(np.abs(finite))
    if median_abs == 0:
        return data, " (dimensionless)"
    exponent = int(np.floor(np.log10(median_abs)))
    if -2 <= exponent <= 2:
        return data, " (dimensionless)"
    scale = 10 ** exponent
    scaled = data / scale
    return scaled, f" (dimensionless, ×10^{exponent})"


def _build_design_matrix(
    x_vals: np.ndarray,
    d_vals: np.ndarray,
    scaler: Optional[StandardScaler],
) -> np.ndarray:
    features = np.column_stack((x_vals ** 2, d_vals * x_vals))
    if scaler is not None:
        features = scaler.transform(features)
    return np.column_stack((np.ones_like(x_vals), features))


def _maybe_scale_features(x_vals: np.ndarray, d_vals: np.ndarray) -> Tuple[np.ndarray, Optional[StandardScaler], float]:
    features = np.column_stack((x_vals ** 2, d_vals * x_vals))
    design = np.column_stack((np.ones_like(x_vals), features))
    cond_number = np.linalg.cond(design)
    scaler: Optional[StandardScaler] = None
    if cond_number > 1e6:
        scaler = StandardScaler()
        scaled_features = scaler.fit_transform(features)
        design = np.column_stack((np.ones_like(x_vals), scaled_features))
    return design, scaler, cond_number


def _solve_ls(A: np.ndarray, y: np.ndarray, scaler: Optional[StandardScaler]) -> Tuple[float, float]:
    bounds = ([-np.inf, 0.0, 0.0], [np.inf, np.inf, np.inf])
    result = lsq_linear(A, y, bounds=bounds)
    if not result.success:
        raise RuntimeError("No se pudo resolver el sistema lineal con restricciones.")
    intercept, gamma_coef, eta_coef = result.x
    if scaler is not None:
        if scaler.scale_[0] != 0:
            gamma_coef = gamma_coef / scaler.scale_[0]
        else:
            gamma_coef = 0.0
        if scaler.scale_[1] != 0:
            eta_coef = eta_coef / scaler.scale_[1]
        else:
            eta_coef = 0.0
    return gamma_coef, eta_coef


def compute_bootstrap_calibration(
    calib_df: pd.DataFrame,
    winsor_cap: float,
    fit_categories: Sequence[str],
    n_bootstrap: int,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    mask = (
        calib_df["sub_network"].str.lower() == "single"
    ) & calib_df["category"].isin(fit_categories)
    work_df = calib_df.loc[mask].dropna(subset=["X", "d", "err_before"]).copy()
    if len(work_df) < 3:
        raise ValueError("No hay suficientes datos para recalcular el bootstrap de calibración.")

    x_vals = np.minimum(work_df["X"].astype(float).to_numpy(), winsor_cap)
    d_vals = work_df["d"].astype(float).to_numpy()
    y_vals = work_df["err_before"].astype(float).to_numpy()

    _, scaler, cond_number = _maybe_scale_features(x_vals, d_vals)

    bounds = ([-np.inf, 0.0, 0.0], [np.inf, np.inf, np.inf])
    idx = np.arange(len(work_df))
    bootstrap_gammas: List[float] = []
    bootstrap_etas: List[float] = []

    for _ in range(n_bootstrap):
        boot_idx = resample(idx)
        xb = np.minimum(work_df.iloc[boot_idx]["X"].astype(float).to_numpy(), winsor_cap)
        db = work_df.iloc[boot_idx]["d"].astype(float).to_numpy()
        yb = work_df.iloc[boot_idx]["err_before"].astype(float).to_numpy()
        A_boot = _build_design_matrix(xb, db, scaler) if scaler else _build_design_matrix(xb, db, None)
        result = lsq_linear(A_boot, yb, bounds=bounds)
        if result.success:
            _, gamma_coef, eta_coef = result.x
            if scaler is not None:
                if scaler.scale_[0] != 0:
                    gamma_coef = gamma_coef / scaler.scale_[0]
                else:
                    gamma_coef = 0.0
                if scaler.scale_[1] != 0:
                    eta_coef = eta_coef / scaler.scale_[1]
                else:
                    eta_coef = 0.0
            bootstrap_gammas.append(gamma_coef)
            bootstrap_etas.append(eta_coef)

    if not bootstrap_gammas or not bootstrap_etas:
        raise RuntimeError("El bootstrap no generó resultados válidos.")

    eta_baseline = float(np.mean(bootstrap_etas))

    loo_records: List[dict] = []
    unique_names = work_df["name"].unique()
    for name in unique_names:
        df_loo = work_df[work_df["name"] != name]
        if len(df_loo) < 3:
            continue
        xl = np.minimum(df_loo["X"].astype(float).to_numpy(), winsor_cap)
        dl = df_loo["d"].astype(float).to_numpy()
        yl = df_loo["err_before"].astype(float).to_numpy()
        A_loo = _build_design_matrix(xl, dl, scaler)
        try:
            gamma_loo, eta_loo = _solve_ls(A_loo, yl, scaler)
        except RuntimeError:
            continue
        delta = eta_loo - eta_baseline
        rel = delta / eta_baseline if eta_baseline != 0 else np.nan
        loo_records.append(
            {
                "name": name,
                "delta_eta": delta,
                "relative_delta": rel,
            }
        )

    loo_df = pd.DataFrame(loo_records).sort_values("relative_delta", ascending=False)
    return np.array(bootstrap_gammas), np.array(bootstrap_etas), loo_df


def _hist_with_kde(ax: plt.Axes, data: np.ndarray, color: str, title: str, xlabel: str) -> None:
    scaled_data, units_label = _scale_axis(data)
    ax.hist(scaled_data, bins=30, color=color, alpha=0.6, edgecolor="black", density=True)
    if np.ptp(scaled_data) > 0:
        xs = np.linspace(scaled_data.min(), scaled_data.max(), 200)
        kde = gaussian_kde(scaled_data)
        ax.plot(xs, kde(xs), color="black", linewidth=1.5)
    ax.set_title(title)
    ax.set_xlabel(f"{xlabel}{units_label}")
    ax.axvline(np.mean(scaled_data), color="black", linestyle="--", linewidth=1)


def fig1_calibration(
    calib_df: pd.DataFrame,
    out_path: Path,
    winsor_cap: float,
    fit_categories: Sequence[str],
    n_bootstrap: int,
):
    gammas, etas, loo_df = compute_bootstrap_calibration(
        calib_df, winsor_cap, fit_categories, n_bootstrap
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    _hist_with_kde(axes[0], etas, "#4C6EF5", "Bootstrap de η", "η")
    _hist_with_kde(axes[1], gammas, "#20C997", "Bootstrap de Γ", "Γ")

    if not loo_df.empty:
        loo_top = loo_df.copy()
        axes[2].barh(
            loo_top["name"],
            loo_top["relative_delta"] * 100,
            color="#F59F00",
        )
        axes[2].invert_yaxis()
        axes[2].set_xlabel("Δη (%) relative to full fit")
        axes[2].set_title("Influence LOO en η")
    else:
        axes[2].text(0.5, 0.5, "Sin datos LOO", ha="center", va="center")
        axes[2].set_axis_off()

    fig.suptitle("Fig. 1 – Calibración de Γ y η", fontsize=18)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def _violin_by_category(
    ax: plt.Axes,
    df: pd.DataFrame,
    categories: Sequence[Tuple[str, str]],
    column: str,
    palette: Sequence[str],
) -> None:
    prepared = [
        (label, df[df["category"] == raw][column].dropna().to_numpy())
        for raw, label in categories
    ]
    prepared = [(cat, values) for cat, values in prepared if len(values) > 0]
    if not prepared:
        ax.text(0.5, 0.5, "Sin datos", ha="center", va="center")
        ax.set_xticks([])
        return
    categories = [item[0] for item in prepared]
    data = [item[1] for item in prepared]
    positions = np.arange(1, len(categories) + 1)
    parts = ax.violinplot(data, positions=positions, showmeans=True, showextrema=False)
    colors = (
        palette
        if len(palette) >= len(categories)
        else plt.cm.viridis(np.linspace(0.2, 0.9, len(categories)))
    )
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(0.7)
    if "cmeans" in parts:
        parts["cmeans"].set_color("black")
    ax.set_xticks(positions)
    ax.set_xticklabels(categories, rotation=30, ha="right")


def fig2_integer_fingerprint(fingerprint_df: pd.DataFrame, out_path: Path):
    families = ["SC_TypeI", "SC_TypeII", "SC_HighPressure", "SC_Binary", "SC_Molecular"]
    cat_pairs = [(cat, CATEGORY_LABELS.get(cat, cat)) for cat in families]
    subset = fingerprint_df[
        (fingerprint_df["lock_family"] == "integer")
        & (fingerprint_df["category"].isin(families))
    ].copy()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    palette1 = plt.cm.viridis(np.linspace(0.2, 0.9, len(families)))
    palette2 = plt.cm.magma(np.linspace(0.2, 0.9, len(families)))

    _violin_by_category(axes[0], subset, cat_pairs, "exp_d_7", palette1)
    axes[0].set_ylabel("exp_d_7")
    axes[0].set_title("Exponentes de 7 (exp_d_7)")

    _violin_by_category(axes[1], subset, cat_pairs, "exp_a_2", palette2)
    axes[1].set_ylabel("exp_a_2")
    axes[1].set_title("Exponentes de 2 (exp_a_2)")

    fig.suptitle("Fig. 2 – Integer fingerprint por familia", fontsize=18)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def fig3_rational_q(fingerprint_df: pd.DataFrame, out_path: Path):
    families = ["Superfluid", "SC_HighPressure", "SC_Oxide", "SC_IronBased"]
    subset = fingerprint_df[
        (fingerprint_df["lock_family"] == "rational")
        & (fingerprint_df["category"].isin(families))
        & fingerprint_df["q"].notna()
    ].copy()
    cat_pairs = [(cat, CATEGORY_LABELS.get(cat, cat)) for cat in families]

    fig, ax = plt.subplots(figsize=(10, 5))
    palette = plt.cm.viridis(np.linspace(0.2, 0.9, len(families)))
    _violin_by_category(ax, subset, cat_pairs, "q", palette)
    ax.set_ylabel("q")
    ax.set_title("Fig. 3 – Distribución de q por familia")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def fig4_residuals(fingerprint_df: pd.DataFrame, out_path: Path):
    df = fingerprint_df.dropna(subset=["log_residual_eta", "category", "sub_network", "d"])

    grouped = (
        df.groupby(["category", "sub_network"], as_index=False)["log_residual_eta"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    grouped["label"] = grouped.apply(
        lambda row: f"{CATEGORY_LABELS.get(row['category'], row['category'])}\n({row['sub_network']})",
        axis=1,
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    colors = ["#74C0FC" if cnt >= 5 else "#FFA94D" for cnt in grouped["count"]]
    axes[0].bar(
        grouped["label"],
        grouped["mean"],
        yerr=grouped["std"],
        color=colors,
        capsize=4,
    )
    axes[0].set_title("Fig. 4a – Residuales medios (log)")
    axes[0].set_ylabel("mean ± std log_residual_eta")
    axes[0].tick_params(axis="x", rotation=60)
    axes[0].legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color="#74C0FC"),
            plt.Rectangle((0, 0), 1, 1, color="#FFA94D"),
        ],
        labels=["N ≥ 5", "N < 5"],
        title="Sample size",
    )

    families = ["SC_HighPressure", "SC_IronBased", "Superfluid"]
    scatter_df = df[df["category"].isin(families)].copy()
    colors = {
        "SC_HighPressure": "#FF6B6B",
        "SC_IronBased": "#5C7CFA",
        "Superfluid": "#51CF66",
    }
    for cat, group in scatter_df.groupby("category"):
        axes[1].scatter(
            group["d"],
            group["log_residual_eta"],
            label=f"{CATEGORY_LABELS.get(cat, cat)} (N={len(group)})",
            alpha=0.8,
            s=50,
            color=colors.get(cat, None),
            edgecolor="black",
        )
    axes[1].set_title("Fig. 4b – Residuo vs d")
    axes[1].set_xlabel("d")
    axes[1].set_ylabel("log_residual_eta")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def fig5_kappa_delta_hist(fingerprint_cluster_df: pd.DataFrame, out_path: Path):
    df = fingerprint_cluster_df.copy()
    if {"err_after_kappa", "err_after_eta"}.issubset(df.columns):
        df["delta_error"] = df["err_after_kappa"].fillna(0) - df["err_after_eta"].fillna(0)
    else:
        df["delta_error"] = np.nan
    delta = df["delta_error"].dropna()
    if delta.empty:
        print("ADVERTENCIA: No se pudieron calcular delta_error para la figura 5.")
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(delta, bins=40, color="#4C6EF5", alpha=0.75, edgecolor="white")
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title("Fig. 5 – Distribución de Δerror (κ − no κ)")
    ax.set_xlabel("Δerror = err_after_kappa − err_after_eta")
    ax.set_ylabel("Número de jumps")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def fig6_kappa_topdelta(fingerprint_cluster_df: pd.DataFrame, out_path: Path):
    df = fingerprint_cluster_df.copy()
    if {"err_after_kappa", "err_after_eta"}.issubset(df.columns):
        df["delta_error"] = df["err_after_kappa"].fillna(0) - df["err_after_eta"].fillna(0)
    else:
        df["delta_error"] = np.nan
    df = df.dropna(subset=["delta_error"])
    if df.empty:
        print("ADVERTENCIA: No se pudieron calcular delta_error para la figura 6.")
        return
    df["abs_delta"] = df["delta_error"].abs()
    top = df.sort_values("abs_delta", ascending=False).head(10).copy()
    def short_label(row):
        base = f"{row['name']} ({row['sub_network']})"
        if isinstance(row["jump_desc"], str) and len(row["jump_desc"]) > 0:
            return f"{base}: {row['jump_desc'].split('→')[-1]}"
        return base
    top["is_mgb2"] = top["name"].str.contains("MgB2", case=False, na=False)
    top["label"] = top.apply(
        lambda r: f"{r['name']} | {r['sub_network']} | {r['jump_desc']}" if r["is_mgb2"] else short_label(r),
        axis=1,
    )
    colors = ["#FF6B6B" if flag else "#8D99AE" for flag in top["is_mgb2"]]
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.barh(top["label"], top["delta_error"], color=colors)
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Δerror = err_after_kappa − err_after_eta")
    ax.set_title("Fig. 6 – Top 10 |Δerror| (κ vs. no κ)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo requerido: {path}")
    return pd.read_csv(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Genera las figuras principales del estudio DOFT")
    parser.add_argument("--base_tag", default=DEFAULT_BASE_TAG, help="Etiqueta base (ej. w600_p10000)")
    parser.add_argument(
        "--results_dir",
        type=str,
        default=str(DATA_PROCESSED),
        help="Directorio con los resultados procesados",
    )
    parser.add_argument(
        "--figures_dir",
        type=str,
        default=None,
        help="Directorio donde guardar las figuras (default: results_<tag>/figures)",
    )
    parser.add_argument(
        "--winsor_cap",
        type=float,
        default=600.0,
        help="Límite superior de Winsor para la calibración",
    )
    parser.add_argument(
        "--fit_categories",
        type=str,
        default=",".join(DEFAULT_FIT_CATEGORIES),
        help="Categorías usadas para calibración (separadas por coma)",
    )
    parser.add_argument(
        "--n_bootstrap",
        type=int,
        default=500,
        help="Número de iteraciones para el bootstrap",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_results_dir = Path(args.results_dir) / f"results_{args.base_tag}"
    if args.figures_dir:
        figures_dir = Path(args.figures_dir)
    else:
        figures_dir = base_results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    calib_path = base_results_dir / "calib" / f"results_calib_{args.base_tag}.csv"
    fingerprint_path = (
        base_results_dir
        / "fingerprint"
        / f"results_fp_kappa_{args.base_tag}_full_factorized.csv"
    )
    cluster_path = base_results_dir / "cluster" / f"results_cluster_kappa_{args.base_tag}.csv"

    calib_df = load_csv(calib_path)
    fingerprint_df = load_csv(fingerprint_path)
    cluster_df = load_csv(cluster_path)

    fit_categories = [item.strip() for item in args.fit_categories.split(",") if item.strip()]

    fig1_calibration(
        calib_df,
        figures_dir / "fig01_calibration.png",
        args.winsor_cap,
        fit_categories,
        args.n_bootstrap,
    )
    fig2_integer_fingerprint(
        fingerprint_df,
        figures_dir / "fig02_integer_fingerprint.png",
    )
    fig3_rational_q(
        fingerprint_df,
        figures_dir / "fig03_rational_q.png",
    )
    fig4_residuals(
        fingerprint_df,
        figures_dir / "fig04_residuals.png",
    )
    fig5_kappa_delta_hist(
        cluster_df,
        figures_dir / "fig05_kappa_delta_hist.png",
    )
    fig6_kappa_topdelta(
        cluster_df,
        figures_dir / "fig06_kappa_topdelta.png",
    )

    print(f"Figuras guardadas en: {figures_dir}")


if __name__ == "__main__":
    main()
