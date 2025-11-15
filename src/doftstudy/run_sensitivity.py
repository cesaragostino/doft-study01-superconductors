import os
import subprocess
import re
import numpy as np
import pandas as pd
import argparse
import sys  # <-- ¡AÑADIDO! Necesario para sys.executable
from pathlib import Path
from tqdm import tqdm

# Definir la constante que faltaba
CONFIG_FILE = "doft_config.json"

def run_sensitivity_analysis(input_csv, n_runs, jitter_percent, prime_max_val=10000, output_dir=None, winsor_value=None):
    """
    Ejecuta el script de análisis de clústeres N veces con jitter
    y reporta las estadísticas de C_AB y kappa.
    """
    
    print(f"--- Iniciando Análisis de Sensibilidad ---")
    print(f"Archivo de entrada: {input_csv}")
    print(f"Iteraciones:        {n_runs}")
    print(f"Jitter:             +/- {jitter_percent:.1f}%")
    print(f"Prime Max Val:      {prime_max_val}")
    print("------------------------------------------")

    TARGET_MATERIALS = [
        'LaH10', 
        'MgB2', 
        '2H-NbSe2', 
        'FeSe', 
        'Helium-3 (He-3) B-phase', 
        'Helium-4 (He-4)'
    ]
    
    KAPPA_REGEX = re.compile(r"--- Coeficiente de Cluster \(kappa\) para (.+?) = ([\d\.e\+-nan]+)")
    
    TEMP_DIR = "temp_sensitivity_run"
    TEMP_CSV = os.path.join(TEMP_DIR, "results_temp.csv")
    os.makedirs(TEMP_DIR, exist_ok=True)
    
    kappa_results = {name: [] for name in TARGET_MATERIALS}
    c_ab_results = {name: [] for name in TARGET_MATERIALS}

    # --- ¡BLOQUE DEFECTUOSO ELIMINADO! ---
    # El bloque 'base_cmd = [...]' que contenía el 'pass'
    # ha sido completamente eliminado.

    # --- Loop de Monte Carlo ---
    print(f"\nEjecutando {n_runs} simulaciones de Monte Carlo (Jittering)...")
    
    # Obtener el comando de python que se está usando para ejecutar ESTE script
    # (ej. 'python3') para llamar al sub-script.
    python_cmd = sys.executable or "python3"
    
    base_cmd = [
        python_cmd, "run_cluster_analysis.py",
        "-i", input_csv,
        "-o", TEMP_DIR,
        "-l", "temp",
        "--estimate_kappa",
        "--jitter_percent", str(jitter_percent),
        "--prime_max_val", str(prime_max_val)
    ]
    
    for _ in tqdm(range(n_runs), desc="Simulaciones"):
        try:
            result = subprocess.run(base_cmd, capture_output=True, text=True, encoding='utf-8', check=False)
            
            if result.returncode != 0:
                print(f"ADVERTENCIA: La simulación falló.")
                print(result.stderr)
                continue

            # 1. Parsear Kappa (κ) desde stdout
            matches = KAPPA_REGEX.findall(result.stdout)
            found_kappas = {}
            for name, kappa_val_str in matches:
                try:
                    kappa_val = float(kappa_val_str)
                    if np.isfinite(kappa_val):
                        found_kappas[name] = kappa_val
                except ValueError:
                    pass # Ignorar 'nan'
            
            for name in kappa_results.keys():
                if name in found_kappas:
                    kappa_results[name].append(found_kappas[name])

            # 2. Parsear C_AB desde el CSV de salida
            if os.path.exists(TEMP_CSV):
                df_temp = pd.read_csv(TEMP_CSV)
                for name in c_ab_results.keys():
                    c_ab_val = df_temp[
                        (df_temp['name'] == name) & (df_temp['jump_type'] == 'intra')
                    ]['C_AB'].max()
                    
                    if pd.notnull(c_ab_val) and np.isfinite(c_ab_val):
                        c_ab_results[name].append(c_ab_val)
            
        except Exception as e:
            print(f"Error durante la simulación: {e}")
            if 'result' in locals():
                print(result.stderr)

    # --- Reporte de Resultados ---
    print("\n--- Resultados del Análisis de Sensibilidad (Jitter) ---")
    
    summary_rows = []
    for name in TARGET_MATERIALS:
        kappas = np.array(kappa_results[name])
        c_abs = np.array(c_ab_results[name])
        
        print(f"\nMaterial: {name}")
        
        if len(kappas) > 0:
            k_mean, k_std = np.mean(kappas), np.std(kappas)
            k_low, k_high = np.percentile(kappas, 2.5), np.percentile(kappas, 97.5)
            print(f"  Kappa (κ) [N={len(kappas)}]:")
            print(f"    Media:   {k_mean:.4g}")
            print(f"    StdDev:  {k_std:.4g}")
            print(f"    95% CI:  [{k_low:.4g}, {k_high:.4g}]")
            summary_rows.append({
                "material": name,
                "metric": "kappa",
                "n": len(kappas),
                "mean": k_mean,
                "std": k_std,
                "ci_low": k_low,
                "ci_high": k_high,
                "prime_max_val": prime_max_val,
                "jitter_percent": jitter_percent,
                "winsor": winsor_value,
            })
        else:
            print("  Kappa (κ): No se pudieron calcular estadísticas (N=0).")
            
        if len(c_abs) > 0:
            c_mean, c_std = np.mean(c_abs), np.std(c_abs)
            c_low, c_high = np.percentile(c_abs, 2.5), np.percentile(c_abs, 97.5)
            print(f"  C_AB [N={len(c_abs)}]:")
            print(f"    Media:   {c_mean:.4f}")
            print(f"    StdDev:  {c_std:.4f}")
            print(f"    95% CI:  [{c_low:.4f}, {c_high:.4f}]")
            summary_rows.append({
                "material": name,
                "metric": "C_AB",
                "n": len(c_abs),
                "mean": c_mean,
                "std": c_std,
                "ci_low": c_low,
                "ci_high": c_high,
                "prime_max_val": prime_max_val,
                "jitter_percent": jitter_percent,
                "winsor": winsor_value,
            })
        else:
            print("  C_AB: No se pudieron calcular estadísticas (N=0).")

    if output_dir and summary_rows:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        tag_j = str(jitter_percent).replace(".", "p").replace("-", "neg")
        p_tag = str(int(prime_max_val))
        summary_file = output_path / f"sensitivity_j{tag_j}_p{p_tag}.csv"
        pd.DataFrame(summary_rows).to_csv(summary_file, index=False)
        print(f"\n--- Resumen de sensibilidad guardado en: {summary_file} ---")

    print("\n--- Limpiando archivos temporales ---")
    try:
        os.remove(TEMP_CSV)
        os.rmdir(TEMP_DIR)
    except Exception as e:
        print(f"No se pudieron limpiar los archivos temporales: {e}")

    print("\n--- ¡Análisis de Sensibilidad Completo! ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DOFT Study - Sensitivity/Jitter Wrapper (v1.2)")
    parser.add_argument("-i", "--input", dest="input_csv", type=str, required=True, help="Input CSV file (e.g., materials_clusters_real_v5.csv)")
    parser.add_argument("-n", "--n_runs", dest="n_runs", type=int, default=500, help="Número de simulaciones de Monte Carlo (default: 500)")
    parser.add_argument("-j", "--jitter", dest="jitter_percent", type=float, default=5.0, help="Porcentaje de jitter +/- aleatorio (e.g., 5.0 para +/- 5%)")
    parser.add_argument("--prime_max_val", dest="prime_max_val", type=int, default=10000, help="Valor máximo para la generación de 'prime values' (default: 10000)")
    parser.add_argument("--output_dir", dest="output_dir", type=str, default=None, help="Directorio donde guardar el resumen CSV (opcional)")
    parser.add_argument("--winsor", dest="winsor_value", type=float, default=None, help="Valor de Winsor usado como baseline (opcional)")
    
    args = parser.parse_args()
    
    if not os.path.exists(CONFIG_FILE):
        print(f"ERROR: No se encontró {CONFIG_FILE}.")
        print("Por favor, ejecuta 'run_calibration.py' (sin jitter) primero para generar este archivo.")
        exit()

    run_sensitivity_analysis(args.input_csv, args.n_runs, args.jitter_percent, args.prime_max_val, args.output_dir, args.winsor_value)
