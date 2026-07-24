"""Genera archivos CSV condensados (digest) a partir de los resultados procesados."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def parse_base_tag(base_dir: Path) -> str:
    name = base_dir.name
    if name.startswith("results_"):
        name = name[len("results_") :]
    return name


def collect_calibration_metadata(base_dirs: List[Path]) -> Tuple[List[dict], Dict[Tuple[float, int], dict]]:
    records: List[dict] = []
    lookup: Dict[Tuple[float, int], dict] = {}

    for base_dir in base_dirs:
        base_tag = parse_base_tag(base_dir)
        meta_path = base_dir / "calib" / f"calibration_metadata_{base_tag}.json"
        if not meta_path.exists():
            alt_path = base_dir / "calib" / f"calibration_metadata_calib_{base_tag}.json"
            if alt_path.exists():
                meta_path = alt_path
        if not meta_path.exists():
            print(f"ADVERTENCIA: No se encontró metadata en {base_dir / 'calib'}.")
            continue
        with open(meta_path, "r") as f:
            meta = json.load(f)
        winsor = float(meta.get("winsor_X_cap"))
        prime = int(meta.get("prime_max_val"))
        eta = meta.get("eta", {})
        gamma = meta.get("gamma", {})
        key = (round(winsor, 6), prime)
        lookup[key] = meta
        records.append(
            {
                "scenario": f"baseline_{base_tag}",
                "winsor": winsor,
                "prime": prime,
                "eta_hat": eta.get("mean"),
                "eta_ci_low": eta.get("ci_low"),
                "eta_ci_high": eta.get("ci_high"),
                "gamma_hat": gamma.get("mean"),
                "gamma_ci_low": gamma.get("ci_low"),
                "gamma_ci_high": gamma.get("ci_high"),
            }
        )
    return records, lookup


def collect_sensitivity_entries(sensitivity_dir: Path, lookup: Dict[Tuple[float, int], dict]) -> List[dict]:
    records: List[dict] = []
    if not sensitivity_dir.exists():
        return records
    for csv_path in sorted(sensitivity_dir.glob("sensitivity_*.csv")):
        scenario = csv_path.stem.replace("sensitivity_", "")
        df = pd.read_csv(csv_path)
        if "prime_max_val" not in df.columns:
            print(f"ADVERTENCIA: El archivo {csv_path} no contiene 'prime_max_val'. Omitiendo.")
            continue
        prime_vals = df["prime_max_val"].dropna().unique()
        winsor_vals = df.get("winsor")
        winsor_val = None
        if winsor_vals is not None:
            winsor_vals = winsor_vals.dropna().unique()
            if len(winsor_vals) > 0:
                winsor_val = float(winsor_vals[0])
        if len(prime_vals) == 0 or winsor_val is None:
            print(f"ADVERTENCIA: No se pudo inferir prime/winsor para {csv_path}.")
            continue
        prime_val = int(prime_vals[0])
        key = (round(winsor_val, 6), prime_val)
        meta = lookup.get(key)
        if not meta:
            print(f"ADVERTENCIA: No se encontró metadata de calibración para winsor={winsor_val}, prime={prime_val}.")
            continue
        eta = meta.get("eta", {})
        gamma = meta.get("gamma", {})
        records.append(
            {
                "scenario": scenario,
                "winsor": winsor_val,
                "prime": prime_val,
                "eta_hat": eta.get("mean"),
                "eta_ci_low": eta.get("ci_low"),
                "eta_ci_high": eta.get("ci_high"),
                "gamma_hat": gamma.get("mean"),
                "gamma_ci_low": gamma.get("ci_low"),
                "gamma_ci_high": gamma.get("ci_high"),
            }
        )
    return records


def summarize_integer_fingerprint(fp_file: Path) -> pd.DataFrame:
    if not fp_file.exists():
        print(f"ADVERTENCIA: No se encontró {fp_file}.")
        return pd.DataFrame()
    df = pd.read_csv(fp_file)
    df_int = df[df["chosen_lock"] == "integer"].copy()
    if df_int.empty:
        return pd.DataFrame()
    records = []
    groups = df_int.groupby(["category", "sub_network"], dropna=False)
    for (category, sub_network), grp in groups:
        record = {
            "category": category,
            "sub_network": sub_network,
            "N": len(grp),
        }
        for col, label in [
            ("exp_a_2", "exp2"),
            ("exp_b_3", "exp3"),
            ("exp_c_5", "exp5"),
            ("exp_d_7", "exp7"),
        ]:
            values = pd.to_numeric(grp[col], errors="coerce")
            record[f"{label}_mean"] = values.mean()
            record[f"{label}_std"] = values.std(ddof=0)
        records.append(record)
    return pd.DataFrame(records)


def summarize_rational_q(fp_file: Path) -> pd.DataFrame:
    if not fp_file.exists():
        return pd.DataFrame()
    df = pd.read_csv(fp_file)
    df_rat = df[df["chosen_lock"] == "rational"].copy()
    if df_rat.empty:
        return pd.DataFrame()
    records = []
    for category, grp in df_rat.groupby("category", dropna=False):
        q_vals = pd.to_numeric(grp["q"], errors="coerce").dropna()
        if q_vals.empty:
            continue
        records.append(
            {
                "category": category,
                "N": len(q_vals),
                "q_mean": q_vals.mean(),
                "q_std": q_vals.std(ddof=0),
                "q_median": q_vals.median(),
                "q_min": q_vals.min(),
                "q_max": q_vals.max(),
            }
        )
    return pd.DataFrame(records)


def summarize_kappa_impact(cluster_file: Path) -> pd.DataFrame:
    if not cluster_file.exists():
        print(f"ADVERTENCIA: No se encontró {cluster_file}.")
        return pd.DataFrame()
    df = pd.read_csv(cluster_file)
    mask = df["err_after_eta"].notna() & df["err_after_kappa"].notna()
    df = df[mask].copy()
    if df.empty:
        return pd.DataFrame()
    df["Delta_error"] = df["err_after_kappa"] - df["err_after_eta"]
    return df[["name", "sub_network", "jump_desc", "category", "Delta_error"]]


def main():
    parser = argparse.ArgumentParser(description="Genera CSVs de digest para los resultados DOFT.")
    parser.add_argument("--output_root", required=True, help="Directorio raíz donde viven los resultados (results_w*_p*).")
    parser.add_argument("--fingerprint_tag", required=True, help="Etiqueta base (wXXX_pYYY) para extraer los fingerprints/kappa.")
    parser.add_argument("--label", default="digest", help="Etiqueta a usar en los nombres de los archivos (default: digest).")
    args = parser.parse_args()

    root = Path(args.output_root).expanduser().resolve()
    digest_dir = root / "digest"
    digest_dir.mkdir(parents=True, exist_ok=True)

    base_dirs = sorted([d for d in root.iterdir() if d.is_dir() and d.name.startswith("results_w")])
    if not base_dirs:
        print(f"ADVERTENCIA: No se encontraron directorios results_w*_p* en {root}.")

    eta_gamma_records, meta_lookup = collect_calibration_metadata(base_dirs)
    eta_gamma_records.extend(collect_sensitivity_entries(root / "sensitivity", meta_lookup))

    label = args.label

    eta_columns = ["scenario", "winsor", "prime", "eta_hat", "eta_ci_low", "eta_ci_high", "gamma_hat", "gamma_ci_low", "gamma_ci_high"]
    eta_gamma_df = pd.DataFrame(eta_gamma_records, columns=eta_columns)
    eta_gamma_df.to_csv(digest_dir / f"eta_gamma_summary_{label}.csv", index=False)
    print(f"--- Archivo 'eta_gamma_summary_{label}.csv' generado con {len(eta_gamma_df)} filas. ---")

    tag = args.fingerprint_tag
    base_dir = root / f"results_{tag}"
    if not base_dir.exists():
        print(f"ADVERTENCIA: No existe el directorio {base_dir}. Se intentará usar el primero disponible.")
        if base_dirs:
            base_dir = base_dirs[0]
            tag = parse_base_tag(base_dir)
        else:
            base_dir = None

    if base_dir is not None and base_dir.exists():
        fp_file = base_dir / "fingerprint" / f"results_fp_kappa_{tag}_full_factorized.csv"
        integer_df = summarize_integer_fingerprint(fp_file)
        integer_df.to_csv(digest_dir / f"integer_fingerprint_summary_{label}.csv", index=False)
        print(f"--- Archivo 'integer_fingerprint_summary_{label}.csv' generado ({len(integer_df)} filas). ---")

        rational_df = summarize_rational_q(fp_file)
        rational_df.to_csv(digest_dir / f"rational_q_summary_{label}.csv", index=False)
        print(f"--- Archivo 'rational_q_summary_{label}.csv' generado ({len(rational_df)} filas). ---")

        cluster_file = base_dir / "cluster" / f"results_cluster_kappa_{tag}.csv"
        kappa_df = summarize_kappa_impact(cluster_file)
        kappa_df.to_csv(digest_dir / f"kappa_impact_{label}.csv", index=False)
        print(f"--- Archivo 'kappa_impact_{label}.csv' generado ({len(kappa_df)} filas). ---")
    else:
        print("ADVERTENCIA: No se pudo determinar un directorio base para fingerprint/kappa.")


if __name__ == "__main__":
    main()
