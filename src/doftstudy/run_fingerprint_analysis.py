import pandas as pd
import numpy as np
import argparse
from fractions import Fraction
import os 
import matplotlib.pyplot as plt

# Importar dependencias
try:
    from sklearn.utils import resample
except ImportError:
    print("ERROR: Faltan dependencias. Ejecuta: pip install scikit-learn")
    exit()

# --- Constantes ---
N_BOOTSTRAPS = 500

# --- Funciones de Factorización (Tu nueva lógica) ---

def get_prime_exponents(n):
    """
    Descompone un entero N en sus exponentes de primos {2, 3, 5, 7}.
    Ej: N=24 -> {a:3, b:1, c:0, d:0} (porque 2^3 * 3^1)
    """
    if not np.isfinite(n) or n <= 0:
        return {'exp_a_2': np.nan, 'exp_b_3': np.nan, 'exp_c_5': np.nan, 'exp_d_7': np.nan, 'exp_other': np.nan}
    n = int(round(n))
    exponents = {'exp_a_2': 0, 'exp_b_3': 0, 'exp_c_5': 0, 'exp_d_7': 0}
    for prime, key in [(2, 'exp_a_2'), (3, 'exp_b_3'), (5, 'exp_c_5'), (7, 'exp_d_7')]:
        while n % prime == 0:
            exponents[key] += 1
            n = n // prime
    exponents['exp_other'] = n if n > 1 else np.nan
    return exponents

def get_rational_pq(val):
    """
    Descompone un float V en su fracción P/Q más simple, con Q <= 8.
    Ej: V=2.625 -> {p: 21, q: 8}
    """
    if not np.isfinite(val): return {'p': np.nan, 'q': np.nan}
    try:
        frac = Fraction(val).limit_denominator(8)
        return {'p': frac.numerator, 'q': frac.denominator}
    except Exception:
        return {'p': np.nan, 'q': np.nan}

def analyze_fingerprints(calib_results_file, cluster_results_file, outdir, run_label):
    """
    Script Maestro 3:
    1. Carga los resultados de los Scripts 1 y 2.
    2. Los combina.
    3. Factoriza los 'prime_value' según 'chosen_lock'.
    4. Genera reportes de 'fingerprint' promediando exponentes/denominadores.
    """
    print(f"--- 1. Cargando archivos de resultados ---")
    try:
        df_calib = pd.read_csv(calib_results_file)
    except FileNotFoundError:
        print(f"ERROR: No se encontró {calib_results_file}. Ejecuta 'run_calibration.py' primero.")
        return
        
    try:
        df_cluster = pd.read_csv(cluster_results_file)
    except FileNotFoundError:
        print(f"ERROR: No se encontró {cluster_results_file}. Ejecuta 'run_cluster_analysis.py' primero.")
        return
        
    df = pd.concat([df_calib, df_cluster], ignore_index=True)
    
    print("\n--- 2. Factorizando 'prime_values' y calculando residuos ---")
    
    factor_records = []
    log_residuals = []
    
    for _, row in df.iterrows():
        pv = row['prime_value']
        lock = row['chosen_lock']
        
        # Calcular Fingerprint de Factores
        if lock == 'integer':
            factor_records.append(get_prime_exponents(pv))
        elif lock == 'rational':
            factor_records.append(get_rational_pq(pv))
        else:
            factor_records.append({'p': np.nan, 'q': np.nan, 'exp_a_2': np.nan, 'exp_b_3': np.nan, 'exp_c_5': np.nan, 'exp_d_7': np.nan, 'exp_other': np.nan})
        
        # Fingerprint Residual (log(R_corr_eta) - log(L*))
        if np.isfinite(row['R_corr_eta']) and row['R_corr_eta'] > 0 and np.isfinite(pv) and pv > 0:
            log_residuals.append(np.log(row['R_corr_eta']) - np.log(pv))
        else:
            log_residuals.append(np.nan)
    
    df_factors = pd.DataFrame(factor_records, index=df.index)
    df['log_residual_eta'] = log_residuals
    df_full = pd.concat([df, df_factors], axis=1)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    results_path = os.path.join(outdir, f"results_{run_label}_full_factorized.csv")
    df_full.to_csv(results_path, index=False, float_format='%.6g')
    print(f"Reporte factorizado completo guardado en: {results_path}")

    # --- 3. Fingerprint Residual (Estacionariedad) ---
    print("\n--- 3. Reporte de Fingerprint Residual (log(R_corr_eta) - log(prime_value)) ---")
    df_residual = df_full.groupby(['category', 'sub_network'])['log_residual_eta'].agg(['mean', 'std', 'count']).reset_index()
    residual_path = os.path.join(outdir, f"fingerprint_{run_label}_log_residual.csv")
    df_residual.to_csv(residual_path, index=False, float_format='%.4g')
    
    # Arreglado: Quitado el floatfmt de to_markdown
    print(df_residual.to_markdown(index=False))

    # --- 4. Bootstrap de Promedios de Familia ---
    print("\n--- 4. Reporte de Fingerprint (Bootstrap CIs por Familia) ---")
    
    groups = df_full.groupby(['category', 'sub_network'])
    bootstrap_results = []
    
    int_cols = ['exp_a_2', 'exp_b_3', 'exp_c_5', 'exp_d_7']
    rat_cols = ['q']
    
    for (category, sub_network), group in groups:
        # A. Bootstrap para 'integer'
        int_data = group[group['chosen_lock'] == 'integer'][int_cols].dropna()
        if not int_data.empty:
            boot_means = {col: [] for col in int_cols}
            for _ in range(N_BOOTSTRAPS):
                sample = resample(int_data)
                for col in int_cols:
                    boot_means[col].append(sample[col].mean())
            
            for col in int_cols:
                ci = np.percentile(boot_means[col], [2.5, 97.5])
                bootstrap_results.append({
                    "family_type": "category+sub_network",
                    "group": f"{category}_{sub_network}",
                    "lock_type": "integer",
                    "metric": col,
                    "n": len(int_data),
                    "mean": np.mean(boot_means[col]),
                    "std_dev": np.std(boot_means[col]),
                    "ci_low": ci[0],
                    "ci_high": ci[1]
                })

        # B. Bootstrap para 'rational'
        rat_data = group[group['chosen_lock'] == 'rational'][rat_cols].dropna()
        if not rat_data.empty:
            boot_means_q = []
            for _ in range(N_BOOTSTRAPS):
                sample = resample(rat_data)
                boot_means_q.append(sample['q'].mean())
            
            ci = np.percentile(boot_means_q, [2.5, 97.5])
            bootstrap_results.append({
                "family_type": "category+sub_network",
                "group": f"{category}_{sub_network}",
                "lock_type": "rational",
                "metric": "q_avg",
                "n": len(rat_data),
                "mean": np.mean(boot_means_q),
                "std_dev": np.std(boot_means_q),
                "ci_low": ci[0],
                "ci_high": ci[1]
            })

    df_boot = pd.DataFrame(bootstrap_results)
    boot_path = os.path.join(outdir, f"fingerprint_{run_label}_bootstrap_CIs.csv")
    df_boot.to_csv(boot_path, index=False, float_format='%.4g')
    print(f"Reporte de Bootstrap CIs guardado en: {boot_path}")
    
    # Arreglado: Quitado el floatfmt de to_markdown
    print(df_boot.to_markdown(index=False))

    # --- 5. Gráficos de Barras con Error ---
    print("\n--- 5. Generando gráficos de barras con error ---")
    
    # Gráfico 1: Exponentes Enteros
    int_plot_data = df_boot[df_boot['lock_type'] == 'integer']
    if not int_plot_data.empty:
        fig, ax = plt.subplots(figsize=(max(10, len(int_plot_data) * 0.5), 6))
        
        x_labels = int_plot_data['group'] + " (" + int_plot_data['metric'] + ")"
        y_vals = int_plot_data['mean']
        y_err = (int_plot_data['ci_high'] - int_plot_data['ci_low']) / 2.0
        
        ax.bar(x_labels, y_vals, yerr=y_err, capsize=5, color='royalblue', alpha=0.7)
        
        ax.set_ylabel("Promedio de Exponente (Bootstrap)")
        ax.set_title("Fingerprint de Locking de Enteros (Media y 95% CI)")
        plt.xticks(rotation=90)
        
        plot_path_int = os.path.join(outdir, f"plot_{run_label}_fingerprint_INTEGER.png")
        fig.savefig(plot_path_int, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Gráfico de 'integer' guardado en: {plot_path_int}")

    # Gráfico 2: Denominadores Racionales
    rat_plot_data = df_boot[df_boot['lock_type'] == 'rational']
    if not rat_plot_data.empty:
        fig, ax = plt.subplots(figsize=(max(8, len(rat_plot_data) * 0.8), 5))
        
        x_labels = rat_plot_data['group']
        y_vals = rat_plot_data['mean']
        y_err = (rat_plot_data['ci_high'] - rat_plot_data['ci_low']) / 2.0
        
        ax.bar(x_labels, y_vals, yerr=y_err, capsize=5, color='firebrick', alpha=0.7)
        
        ax.set_ylabel("Promedio de Denominador 'q' (Bootstrap)")
        ax.set_title("Fingerprint de Locking Racional (Media y 95% CI)")
        plt.xticks(rotation=45, ha='right')
        
        plot_path_rat = os.path.join(outdir, f"plot_{run_label}_fingerprint_RATIONAL.png")
        fig.savefig(plot_path_rat, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Gráfico de 'rational' guardado en: {plot_path_rat}")

    print("\n--- ¡Análisis de Fingerprint completo! ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DOFT Study - Fingerprint Analyzer (v-Master 3.2)")
    
    # --- ¡ARREGLO! Se añadieron los flags largos '--' ---
    parser.add_argument("-i_calib", "--input_calib", dest="input_calib", type=str, required=True, help="Input CSV from run_calibration.py (e.g., results_final_calib.csv)")
    parser.add_argument("-i_cluster", "--input_cluster", dest="input_cluster", type=str, required=True, help="Input CSV from run_cluster_analysis.py (e.g., results_cluster_fingerprints.csv)")
    
    parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True, help="Output directory (e.g., 'results_fingerprints')")
    parser.add_argument("-l", "--label", dest="run_label", type=str, default="fingerprint", help="Run label for output files")
    
    args = parser.parse_args()
    
    try:
        analyze_fingerprints(args.input_calib, args.input_cluster, args.outdir, args.run_label)
            
    except Exception as e:
        print(f"Ha ocurrido un error inesperado: {e}")