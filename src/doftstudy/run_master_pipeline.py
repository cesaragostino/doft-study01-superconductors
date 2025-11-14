import subprocess
import os
import sys
import argparse
from itertools import product

# -----------------------------------------------------------------
# --- CONFIGURACIÓN "RUN ALL" ---
# Define aquí las listas completas para el análisis
# -----------------------------------------------------------------
ALL_WINSOR_VALUES = [400.0, 600.0, 800.0]
ALL_PRIME_VALUES = [10000, 7919]
ALL_JITTER_VALUES = [5.0, 10.0]
DEFAULT_JITTER_N_RUNS = 500
# -----------------------------------------------------------------

def execute_command(cmd, step_name):
    """
    Ejecuta un comando en el shell, comprueba si hay errores
    e imprime información útil.
    """
    print(f"\n[... {step_name} ...]")
    print(" ".join(cmd))
    
    # Ejecuta el comando y captura la salida
    result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
    
    # --- ¡CAMBIO v3.1! ---
    # Imprimir stdout y stderr SIEMPRE, para que el log los capture.
    
    # Imprimir stdout si no está vacío
    if result.stdout:
        print("--- Salida Estándar (stdout): ---")
        print(result.stdout)
    # Imprimir stderr si no está vacío
    if result.stderr:
        print("--- Salida de Error (stderr): ---")
        print(result.stderr)
    
    # Ahora, comprobar el código de retorno
    if result.returncode != 0:
        print(f"¡ERROR! El pipeline falló en: [{step_name}]")
        # El error ya fue impreso arriba
        sys.exit(f"Deteniendo el script debido a un error en {step_name}.")
    
    print(f"[¡OK! {step_name}]")


# -----------------------------------------------------------------
# FASE 1: PIPELINE DE ROBUSTEZ (Workflow Bifurcado)
# -----------------------------------------------------------------
def run_robustness_pipeline(w_val, p_val, input_file, python_cmd):
    """
    Ejecuta el pipeline COMPLETO de 5 pasos (bifurcado con/sin kappa)
    para una combinación de Winsor (w_val) y Prime Max (p_val).
    """
    print(f"=====================================================")
    print(f"   INICIANDO PIPELINE ROBUSTEZ (W={w_val}, P={p_val})   ")
    print(f"=====================================================")
    
    w_tag = f"w{int(w_val)}"
    p_tag = f"p{int(p_val)}"
    base_tag = f"{w_tag}_{p_tag}"
    out_dir_base = f"results_{base_tag}"
    
    out_dir_calib = os.path.join(out_dir_base, "calib")
    out_dir_cluster = os.path.join(out_dir_base, "cluster")
    out_dir_fp = os.path.join(out_dir_base, "fingerprint")
    
    os.makedirs(out_dir_calib, exist_ok=True)
    os.makedirs(out_dir_cluster, exist_ok=True)
    os.makedirs(out_dir_fp, exist_ok=True)

    # Paso 1: Calibración
    label_calib = f"calib_{base_tag}"
    csv_calib = os.path.join(out_dir_calib, f"results_{label_calib}.csv")
    cmd1 = [
        python_cmd, "run_calibration.py",
        "-i", input_file, "-o", out_dir_calib, "-l", label_calib,
        "--winsor_X", str(w_val), "--prime_max_val", str(p_val)
    ]
    execute_command(cmd1, f"Paso 1: Calibración (W={w_val}, P={p_val})")
    
    # Paso 2a: Cluster (SIN Kappa)
    label_cluster_nokappa = f"cluster_nokappa_{base_tag}"
    csv_cluster_nokappa = os.path.join(out_dir_cluster, f"results_{label_cluster_nokappa}.csv")
    cmd2a = [
        python_cmd, "run_cluster_analysis.py",
        "-i", input_file, "-o", out_dir_cluster, "-l", label_cluster_nokappa,
        "--prime_max_val", str(p_val)
    ]
    execute_command(cmd2a, f"Paso 2a: Cluster SIN Kappa (W={w_val}, P={p_val})")

    # Paso 2b: Cluster (CON Kappa)
    label_cluster_kappa = f"cluster_kappa_{base_tag}"
    csv_cluster_kappa = os.path.join(out_dir_cluster, f"results_{label_cluster_kappa}.csv")
    cmd2b = [
        python_cmd, "run_cluster_analysis.py",
        "-i", input_file, "-o", out_dir_cluster, "-l", label_cluster_kappa,
        "--estimate_kappa", "--prime_max_val", str(p_val)
    ]
    execute_command(cmd2b, f"Paso 2b: Cluster CON Kappa (W={w_val}, P={p_val})")

    # Paso 3a: Fingerprint (SIN Kappa)
    label_fp_nokappa = f"fp_nokappa_{base_tag}"
    cmd3a = [
        python_cmd, "run_fingerprint_analysis.py",
        "-i_calib", csv_calib, "-i_cluster", csv_cluster_nokappa,
        "-o", out_dir_fp, "-l", label_fp_nokappa
    ]
    execute_command(cmd3a, f"Paso 3a: Fingerprint SIN Kappa (W={w_val}, P={p_val})")

    # Paso 3b: Fingerprint (CON Kappa)
    label_fp_kappa = f"fp_kappa_{base_tag}"
    cmd3b = [
        python_cmd, "run_fingerprint_analysis.py",
        "-i_calib", csv_calib, "-i_cluster", csv_cluster_kappa,
        "-o", out_dir_fp, "-l", label_fp_kappa
    ]
    execute_command(cmd3b, f"Paso 3b: Fingerprint CON Kappa (W={w_val}, P={p_val})")
    
    print(f"--- Pipeline Robustez (W={w_val}, P={p_val}) ¡COMPLETADO! ---")

# -----------------------------------------------------------------
# FUNCIÓN AUXILIAR PARA FASE 2
# -----------------------------------------------------------------
def run_calibration_only(w_val, p_val, input_file, python_cmd):
    """
    Ejecuta SOLO el Paso 1 (Calibración) para generar el 
    archivo 'doft_config.json' necesario para el jittering.
    Lo ejecuta en el directorio base para esa config.
    """
    print(f"--- Asegurando 'doft_config.json' base (W={w_val}, P={p_val}) ---")
    
    w_tag = f"w{int(w_val)}"
    p_tag = f"p{int(p_val)}"
    
    out_dir_calib = os.path.join(f"results_{w_tag}_{p_tag}", "calib")
    label_calib = f"calib_{w_tag}_{p_tag}"
    os.makedirs(out_dir_calib, exist_ok=True)

    cmd_calib = [
        python_cmd, "run_calibration.py",
        "-i", input_file, "-o", out_dir_calib, "-l", label_calib,
        "--winsor_X", str(w_val), "--prime_max_val", str(p_val)
    ]
    execute_command(cmd_calib, f"Calibración Base (W={w_val}, P={p_val})")
    print("--- 'doft_config.json' generado/actualizado. ---")


# -----------------------------------------------------------------
# FASE 2: PIPELINE DE SENSIBILIDAD (JITTER)
# -----------------------------------------------------------------
def run_sensitivity_analysis(j_val, p_val, input_file, n_runs, python_cmd):
    """
    Ejecuta el script de sensibilidad 'run_sensitivity.py'
    para una combinación de Jitter (j_val) y Prime Max (p_val).
    """
    print(f"=====================================================")
    print(f"   INICIANDO PIPELINE SENSIBILIDAD (J={j_val}%, P={p_val})   ")
    print(f"=====================================================")
    
    cmd_jitter = [
        python_cmd, "run_sensitivity.py",
        "-i", input_file,
        "-n", str(n_runs),
        "-j", str(j_val),
        "--prime_max_val", str(p_val)
    ]
    execute_command(cmd_jitter, f"Paso S: Jittering (J={j_val}, P={p_val})")
    print(f"--- Pipeline Sensibilidad (J={j_val}, P={p_val}) ¡COMPLETADO! ---")


# -----------------------------------------------------------------
# SCRIPT PRINCIPAL
# -----------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Script 'Maestro' (v3.1) para ejecutar todos los análisis (Robustez y Sensibilidad).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # --- Flag "RUN ALL" ---
    parser.add_argument(
        '--run_all',
        action='store_true',
        help="¡ACTIVA EL MODO 'RUN ALL'!\nEjecuta TODAS las combinaciones predefinidas (ignora -w, -p, -j)."
    )
    
    # --- Argumentos de Configuración General ---
    parser.add_argument(
        '-i', '--input_file',
        type=str,
        default="materials_clusters_real_v5.csv",
        help="Archivo CSV de entrada (default: %(default)s)"
    )
    parser.add_argument(
        '--python_cmd',
        type=str,
        default="python3",
        help="Comando para ejecutar Python (ej. 'python3' en Mac, 'python' en Win)\n(default: %(default)s)"
    )
    
    # --- Argumentos para corridas CUSTOM (se ignoran con --run_all) ---
    parser.add_argument(
        '-w', '--winsor_values',
        nargs='+', type=float, default=[600.0],
        help="Valores de Winsor (default: [600.0])"
    )
    parser.add_argument(
        '-p', '--prime_values',
        nargs='+', type=int, default=[10000],
        help="Valores 'prime set' (default: [10000])"
    )
    parser.add_argument(
        '-j', '--jitter_values',
        nargs='+', type=float, default=[],
        help="Valores de Jitter (default: [])"
    )
    parser.add_argument(
        '-n', '--n_runs',
        type=int, default=DEFAULT_JITTER_N_RUNS,
        help="Número de iteraciones para el jitter (default: %(default)s)"
    )
    
    # --- Flags para saltar fases ---
    parser.add_argument(
        '--skip_robustness',
        action='store_true',
        help="Saltar la Fase 1 (Robustez: Calib+Cluster+FP)"
    )
    parser.add_argument(
        '--skip_sensitivity',
        action='store_true',
        help="Saltar la Fase 2 (Sensibilidad: Jittering)"
    )

    args = parser.parse_args()
    
    # --- Lógica principal: Decidir qué listas usar ---
    
    if args.run_all:
        print("============================================")
        print("     ¡MODO 'RUN ALL' ACTIVADO!     ")
        print("Se ejecutarán todas las combinaciones predefinidas.")
        print("============================================")
        winsor_list   = ALL_WINSOR_VALUES
        prime_list    = ALL_PRIME_VALUES
        jitter_list   = ALL_JITTER_VALUES
        n_runs        = DEFAULT_JITTER_N_RUNS
    else:
        print("============================================")
        print("     Modo de ejecución 'Custom'     ")
        print("Ejecutando con los parámetros provistos (o defaults).")
        print("============================================")
        winsor_list   = args.winsor_values
        prime_list    = args.prime_values
        jitter_list   = args.jitter_values
        n_runs        = args.n_runs

    print(f"Comando Python:   {args.python_cmd}")
    print(f"Archivo de entrada: {args.input_file}")
    print(f"Valores Winsor:   {winsor_list}")
    print(f"Valores Prime:    {prime_list}")
    print(f"Valores Jitter:   {jitter_list}")
    print(f"Jitter N-Runs:    {n_runs}")
    print("--------------------------------------------")

    # --- FASE 1: ROBUSTEZ (Workflow Bifurcado) ---
    if not args.skip_robustness:
        print("\n\n--- FASE 1: INICIANDO ANÁLISIS DE ROBUSTEZ (Workflow Completo) ---")
        robustness_jobs = list(product(winsor_list, prime_list))
        print(f"Se ejecutarán {len(robustness_jobs)} pipelines de robustez (de 5 pasos c/u)...")
        
        for w_val, p_val in robustness_jobs:
            try:
                run_robustness_pipeline(w_val, p_val, args.input_file, args.python_cmd)
            except Exception as e:
                print(f"¡¡ERROR INESPERADO!! El pipeline (W={w_val}, P={p_val}) falló.")
                print(f"Detalle del error: {e}")
                
    else:
        print("\n--- FASE 1: Saltando Análisis de Robustez (por flag). ---")

    # --- FASE 2: SENSIBILIDAD (JITTER) ---
    if not args.skip_sensitivity and jitter_list:
        print("\n\n--- FASE 2: INICIANDO ANÁLISIS DE SENSIBILIDAD ---")
        
        baseline_w = winsor_list[0]
        baseline_p = prime_list[0]
        config_path = "doft_config.json"
        
        if not os.path.exists(config_path) or args.skip_robustness:
             try:
                print("El 'doft_config.json' no existe o se saltó la Fase 1. Generando config base...")
                run_calibration_only(baseline_w, baseline_p, args.input_file, args.python_cmd)
             except Exception as e:
                print("¡¡ERROR FATAL!! No se pudo generar el 'doft_config.json' base.")
                print(f"Detalle del error: {e}")
                sys.exit(1)
        else:
            print(f"--- 'doft_config.json' ya existe. Usando como base para Jitter. ---")
            
        sensitivity_jobs = list(product(jitter_list, prime_list))
        print(f"Se ejecutarán {len(sensitivity_jobs)} pipelines de sensibilidad...")

        for j_val, p_val in sensitivity_jobs:
            try:
                run_sensitivity_analysis(j_val, p_val, args.input_file, n_runs, args.python_cmd)
            except Exception as e:
                print(f"¡¡ERROR INESPERADO!! El pipeline (J={j_val}, P={p_val}) falló.")
                print(f"Detalle del error: {e}")
                
    elif args.skip_sensitivity:
        print("\n--- FASE 2: Saltando Análisis de Sensibilidad (por flag). ---")
    else:
        print("\n--- FASE 2: No se especificaron --jitter_values, saltando análisis de sensibilidad. ---")
        
    print("\n\n============================================")
    print("   ¡Master Pipeline finalizado!   ")
    print("============================================")


if __name__ == "__main__":
    main()