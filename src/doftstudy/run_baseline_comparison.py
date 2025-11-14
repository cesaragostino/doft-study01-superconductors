"""Comparación simple contra un modelo baseline basado en power-laws."""

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


DEFAULT_GROUPINGS = ["lock_family", "jump_type", "category", "family_type"]
EPS = 1e-12


def detect_error_column(df: pd.DataFrame) -> str:
    """Selecciona la columna de error DOFT más informativa disponible."""

    for candidate in ("err_after_kappa", "err_after_eta", "err_before"):
        if candidate in df.columns:
            return candidate
    raise ValueError("No se encontró ninguna columna de error DOFT en el CSV de cluster.")


def fit_power_law(d_vals: np.ndarray, r_vals: np.ndarray) -> Tuple[float, float, np.ndarray]:
    """Ajusta R ≈ A * d^alpha en log-log y retorna (alpha, A, predicciones)."""

    log_d = np.log(d_vals)
    log_r = np.log(r_vals)

    if np.allclose(log_d, log_d[0]):
        # Si todos los 'd' son iguales el mejor ajuste es constante (alpha=0).
        alpha = 0.0
        intercept = log_r.mean()
    else:
        alpha, intercept = np.polyfit(log_d, log_r, 1)

    scale = np.exp(intercept)
    preds = scale * np.power(d_vals, alpha)
    return alpha, scale, preds


def summarize_group(group_name: str, df_group: pd.DataFrame, error_col: str, min_samples: int) -> Optional[dict]:
    """Computa métricas de baseline vs DOFT para un subconjunto."""

    valid_mask = (df_group['R_obs'] > 0) & (df_group['d'] > 0)
    df_valid = df_group.loc[valid_mask].copy()

    if len(df_valid) < max(min_samples, 2):
        return None

    d_vals = df_valid['d'].astype(float).to_numpy()
    r_vals = df_valid['R_obs'].astype(float).to_numpy()

    try:
        alpha, scale, preds = fit_power_law(d_vals, r_vals)
    except np.linalg.LinAlgError:
        return None

    baseline_mare = np.mean(np.abs(r_vals - preds) / np.maximum(np.abs(r_vals), EPS))
    baseline_log_mae = np.mean(np.abs(np.log(r_vals) - np.log(preds)))

    doft_mask = df_valid['prime_value'].notna() if 'prime_value' in df_valid else pd.Series(False, index=df_valid.index)
    if doft_mask.any():
        doft_vals = df_valid.loc[doft_mask]
        doft_mare = np.mean(
            np.abs(doft_vals['prime_value'] - doft_vals['R_obs']) /
            np.maximum(np.abs(doft_vals['R_obs']), EPS)
        )
    else:
        doft_mare = np.nan

    return {
        "group_value": group_name,
        "n_total": int(len(df_group)),
        "n_valid": int(len(df_valid)),
        "baseline_alpha": float(alpha),
        "baseline_scale": float(scale),
        "baseline_mare": float(baseline_mare),
        "baseline_log_mae": float(baseline_log_mae),
        "doft_mare": float(doft_mare) if pd.notnull(doft_mare) else np.nan,
        "baseline_params": 2,
        "doft_params": 0,
        "err_delta_baseline_minus_doft": float(baseline_mare - doft_mare) if pd.notnull(doft_mare) else np.nan,
        "error_metric": "mean_absolute_relative_error",
        "doft_error_col": error_col,
    }


def build_summary(df: pd.DataFrame, groupings: List[str], min_samples: int) -> pd.DataFrame:
    """Genera un DataFrame con métricas por agrupación."""

    error_col = detect_error_column(df)
    records: List[dict] = []

    for grouping in groupings:
        if grouping not in df.columns:
            print(f"ADVERTENCIA: La columna '{grouping}' no existe en el CSV. Se omite esa agrupación.")
            continue

        for group_value, df_group in df.groupby(grouping):
            if pd.isna(group_value):
                continue
            metrics = summarize_group(str(group_value), df_group, error_col, min_samples)
            if not metrics:
                continue
            metrics.update({
                "grouping": grouping,
            })
            records.append(metrics)

    if not records:
        return pd.DataFrame()

    summary_df = pd.DataFrame.from_records(records)
    cols_order = [
        "grouping", "group_value", "error_metric",
        "n_total", "n_valid",
        "baseline_scale", "baseline_alpha",
        "baseline_mare", "baseline_log_mae",
        "doft_mare", "doft_error_col",
        "baseline_params", "doft_params", "err_delta_baseline_minus_doft",
    ]
    summary_df = summary_df[cols_order]
    summary_df.sort_values(["grouping", "baseline_mare"], inplace=True)
    return summary_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Comparación DOFT vs baseline power-law")
    parser.add_argument("--cluster_csv", "-i", required=True, help="CSV de resultados de cluster (usar el archivo *_kappa)")
    parser.add_argument("--output", "-o", required=True, help="Ruta del CSV resumen a generar")
    parser.add_argument(
        "--group_by", nargs='+', default=DEFAULT_GROUPINGS,
        help=f"Columnas para agrupar (default: {' '.join(DEFAULT_GROUPINGS)})"
    )
    parser.add_argument("--min_samples", type=int, default=3, help="Mínimo de registros válidos por grupo para ajustar el baseline")

    args = parser.parse_args()

    cluster_path = Path(args.cluster_csv)
    if not cluster_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de cluster: {cluster_path}")

    df_cluster = pd.read_csv(cluster_path)
    if 'R_obs' not in df_cluster.columns:
        raise ValueError("El CSV de cluster no contiene la columna 'R_obs'.")
    if 'd' not in df_cluster.columns:
        raise ValueError("El CSV de cluster no contiene la columna 'd'.")

    summary = build_summary(df_cluster, args.group_by, args.min_samples)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)

    print(f"--- Comparación baseline guardada en: {output_path} ---")


if __name__ == "__main__":
    main()
