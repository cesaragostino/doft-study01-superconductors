import subprocess
import os
import sys
import argparse
from itertools import product
from pathlib import Path
from datetime import datetime
from contextlib import contextmanager
import pandas as pd

# -----------------------------------------------------------------
# --- CONFIGURACIÓN "RUN ALL" ---
# Define aquí las listas completas para el análisis
# -----------------------------------------------------------------
ALL_WINSOR_VALUES = [400.0, 600.0, 800.0]
ALL_PRIME_VALUES = [10000, 7919]
ALL_JITTER_VALUES = [5.0, 10.0]
DEFAULT_JITTER_N_RUNS = 500
# -----------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_INPUT_FILE = PROJECT_ROOT / "data" / "raw" / "materials_clusters_real_v6.csv"
RUN_CALIBRATION_SCRIPT = SCRIPT_DIR / "run_calibration.py"
RUN_CLUSTER_SCRIPT = SCRIPT_DIR / "run_cluster_analysis.py"
RUN_FP_SCRIPT = SCRIPT_DIR / "run_fingerprint_analysis.py"
RUN_SENSITIVITY_SCRIPT = SCRIPT_DIR / "run_sensitivity.py"
RUN_BASELINE_SCRIPT = SCRIPT_DIR / "run_baseline_comparison.py"
RUN_FIGURES_SCRIPT = SCRIPT_DIR / "generate_figures.py"
RUN_DIGEST_SCRIPT = SCRIPT_DIR / "generate_digest.py"

def sanitize_tag(value):
    value_str = str(value)
    if value_str.endswith(".0"):
        value_str = value_str[:-2]
    return value_str.replace(".", "p").replace("-", "neg")

def infer_digest_label(input_path):
    stem = Path(input_path).stem.lower()
    for token in ["v6", "v5", "v4"]:
        if token in stem:
            return token
    return stem


def sensitivity_summary_complete(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8") as f:
            header = f.readline().strip().split(",")
        return "prime_max_val" in header and "winsor" in header
    except Exception:
        return False


class TeeStream:
    """Duplica stdout/stderr hacia un archivo de log."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for stream in self._streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self._streams:
            stream.flush()


@contextmanager
def tee_output(log_path):
    """Envía stdout/stderr tanto a consola como a un archivo de log."""

    original_stdout, original_stderr = sys.stdout, sys.stderr
    log_file = open(log_path, "a", encoding="utf-8")
    sys.stdout = TeeStream(original_stdout, log_file)
    sys.stderr = TeeStream(original_stderr, log_file)
    try:
        yield
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()

def execute_command(cmd, step_name, cwd=SCRIPT_DIR):
    """
    Ejecuta un comando en el shell, comprueba si hay errores
    e imprime información útil.
    """
    print(f"\n[... {step_name} ...]")
    print(" ".join(cmd))

    # Ejecuta el comando y captura la salida
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding='utf-8',
        cwd=str(cwd) if cwd else None,
    )
    
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
def run_robustness_pipeline(w_val, p_val, input_file, python_cmd, output_root):
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
    out_dir_base = output_root / f"results_{base_tag}"
    
    out_dir_calib = out_dir_base / "calib"
    out_dir_cluster = out_dir_base / "cluster"
    out_dir_fp = out_dir_base / "fingerprint"
    out_dir_baseline = out_dir_base / "baseline"
    out_dir_figures = out_dir_base / "figures"
    
    out_dir_calib.mkdir(parents=True, exist_ok=True)
    out_dir_cluster.mkdir(parents=True, exist_ok=True)
    out_dir_fp.mkdir(parents=True, exist_ok=True)
    out_dir_baseline.mkdir(parents=True, exist_ok=True)
    out_dir_figures.mkdir(parents=True, exist_ok=True)

    # Rutas de salida esperadas
    label_calib = f"calib_{base_tag}"
    csv_calib = out_dir_calib / f"results_{label_calib}.csv"
    summary_calib = out_dir_calib / f"summary_{label_calib}.csv"
    meta_calib = out_dir_calib / f"calibration_metadata_{label_calib}.json"
    label_cluster_nokappa = f"cluster_nokappa_{base_tag}"
    csv_cluster_nokappa = out_dir_cluster / f"results_{label_cluster_nokappa}.csv"
    label_cluster_kappa = f"cluster_kappa_{base_tag}"
    csv_cluster_kappa = out_dir_cluster / f"results_{label_cluster_kappa}.csv"
    label_fp_nokappa = f"fp_nokappa_{base_tag}"
    fp_nokappa_csv = out_dir_fp / f"results_{label_fp_nokappa}_full_factorized.csv"
    label_fp_kappa = f"fp_kappa_{base_tag}"
    fp_kappa_csv = out_dir_fp / f"results_{label_fp_kappa}_full_factorized.csv"
    baseline_summary = out_dir_baseline / f"baseline_summary_{base_tag}.csv"
    figure_targets = [
        out_dir_figures / "fig01_calibration.png",
        out_dir_figures / "fig02_integer_fingerprint.png",
        out_dir_figures / "fig03_rational_q.png",
        out_dir_figures / "fig04_residuals.png",
        out_dir_figures / "fig05_kappa_delta_hist.png",
        out_dir_figures / "fig06_kappa_topdelta.png",
    ]

    existing_outputs = [
        csv_calib,
        summary_calib,
        meta_calib,
        csv_cluster_nokappa,
        csv_cluster_kappa,
        fp_nokappa_csv,
        fp_kappa_csv,
        baseline_summary,
        *figure_targets,
    ]
    if all(path.exists() for path in existing_outputs):
        print(f"--- Resultados detectados para {base_tag}. Saltando pipeline completo. ---")
        return

    # Paso 1: Calibración
    if csv_calib.exists() and summary_calib.exists() and meta_calib.exists():
        print(f"--- Calibración detectada para {base_tag}. Saltando Paso 1. ---")
    else:
        cmd1 = [
            python_cmd, str(RUN_CALIBRATION_SCRIPT),
            "-i", input_file, "-o", str(out_dir_calib), "-l", label_calib,
            "--winsor_X", str(w_val), "--prime_max_val", str(p_val)
        ]
        execute_command(cmd1, f"Paso 1: Calibración (W={w_val}, P={p_val})")
    
    # Paso 2a: Cluster (SIN Kappa)
    if csv_cluster_nokappa.exists():
        print(f"--- Cluster sin κ detectado para {base_tag}. Saltando Paso 2a. ---")
    else:
        cmd2a = [
            python_cmd, str(RUN_CLUSTER_SCRIPT),
            "-i", input_file, "-o", str(out_dir_cluster), "-l", label_cluster_nokappa,
            "--prime_max_val", str(p_val)
        ]
        execute_command(cmd2a, f"Paso 2a: Cluster SIN Kappa (W={w_val}, P={p_val})")

    # Paso 2b: Cluster (CON Kappa)
    if csv_cluster_kappa.exists():
        print(f"--- Cluster con κ detectado para {base_tag}. Saltando Paso 2b. ---")
    else:
        cmd2b = [
            python_cmd, str(RUN_CLUSTER_SCRIPT),
            "-i", input_file, "-o", str(out_dir_cluster), "-l", label_cluster_kappa,
            "--estimate_kappa", "--prime_max_val", str(p_val)
        ]
        execute_command(cmd2b, f"Paso 2b: Cluster CON Kappa (W={w_val}, P={p_val})")

    # Paso 3a: Fingerprint (SIN Kappa)
    if fp_nokappa_csv.exists():
        print(f"--- Fingerprint sin κ detectado para {base_tag}. Saltando Paso 3a. ---")
    else:
        cmd3a = [
            python_cmd, str(RUN_FP_SCRIPT),
            "-i_calib", str(csv_calib), "-i_cluster", str(csv_cluster_nokappa),
            "-o", str(out_dir_fp), "-l", label_fp_nokappa
        ]
        execute_command(cmd3a, f"Paso 3a: Fingerprint SIN Kappa (W={w_val}, P={p_val})")

    # Paso 3b: Fingerprint (CON Kappa)
    if fp_kappa_csv.exists():
        print(f"--- Fingerprint con κ detectado para {base_tag}. Saltando Paso 3b. ---")
    else:
        cmd3b = [
            python_cmd, str(RUN_FP_SCRIPT),
            "-i_calib", str(csv_calib), "-i_cluster", str(csv_cluster_kappa),
            "-o", str(out_dir_fp), "-l", label_fp_kappa
        ]
        execute_command(cmd3b, f"Paso 3b: Fingerprint CON Kappa (W={w_val}, P={p_val})")
    
    # Paso 4: Comparación Baseline
    if csv_cluster_kappa.exists():
        if baseline_summary.exists():
            print(f"--- Resumen baseline detectado para {base_tag}. Saltando Paso 4. ---")
        else:
            cmd4 = [
                python_cmd, str(RUN_BASELINE_SCRIPT),
                "--cluster_csv", str(csv_cluster_kappa),
                "--output", str(baseline_summary)
            ]
            execute_command(cmd4, f"Paso 4: Comparación Baseline (W={w_val}, P={p_val})")
    else:
        print("--- Resultados de cluster (κ) ausentes; saltando comparación baseline. ---")
    
    # Paso 5: Figuras para la corrida actual
    if fp_kappa_csv.exists():
        if all(path.exists() for path in figure_targets):
            print(f"--- Figuras detectadas para {base_tag}. Saltando Paso 5. ---")
        else:
            cmd5 = [
                python_cmd, str(RUN_FIGURES_SCRIPT),
                "--base_tag", base_tag,
                "--results_dir", str(output_root)
            ]
            execute_command(cmd5, f"Paso 5: Figuras (W={w_val}, P={p_val})")
    else:
        print("--- Fingerprint con κ ausente; saltando generación de figuras. ---")
    
    print(f"--- Pipeline Robustez (W={w_val}, P={p_val}) ¡COMPLETADO! ---")

# -----------------------------------------------------------------
# FUNCIÓN AUXILIAR PARA FASE 2
# -----------------------------------------------------------------
def run_calibration_only(w_val, p_val, input_file, python_cmd, output_root):
    """
    Ejecuta SOLO el Paso 1 (Calibración) para generar el 
    archivo 'doft_config.json' necesario para el jittering.
    Lo ejecuta en el directorio base para esa config.
    """
    print(f"--- Asegurando 'doft_config.json' base (W={w_val}, P={p_val}) ---")
    
    w_tag = f"w{int(w_val)}"
    p_tag = f"p{int(p_val)}"
    
    out_dir_calib = output_root / f"results_{w_tag}_{p_tag}" / "calib"
    label_calib = f"calib_{w_tag}_{p_tag}"
    out_dir_calib.mkdir(parents=True, exist_ok=True)

    cmd_calib = [
        python_cmd, str(RUN_CALIBRATION_SCRIPT),
        "-i", input_file, "-o", str(out_dir_calib), "-l", label_calib,
        "--winsor_X", str(w_val), "--prime_max_val", str(p_val)
    ]
    execute_command(cmd_calib, f"Calibración Base (W={w_val}, P={p_val})")
    print("--- 'doft_config.json' generado/actualizado. ---")


# -----------------------------------------------------------------
# FASE 2: PIPELINE DE SENSIBILIDAD (JITTER)
# -----------------------------------------------------------------
def run_sensitivity_analysis(j_val, p_val, input_file, n_runs, python_cmd, summary_dir, winsor_value):
    """
    Ejecuta el script de sensibilidad 'run_sensitivity.py'
    para una combinación de Jitter (j_val) y Prime Max (p_val).
    """
    print(f"=====================================================")
    print(f"   INICIANDO PIPELINE SENSIBILIDAD (J={j_val}%, P={p_val})   ")
    print(f"=====================================================")
    
    cmd_jitter = [
        python_cmd, str(RUN_SENSITIVITY_SCRIPT),
        "-i", input_file,
        "-n", str(n_runs),
        "-j", str(j_val),
        "--prime_max_val", str(p_val),
        "--output_dir", str(summary_dir)
    ]
    if winsor_value is not None:
        cmd_jitter.extend(["--winsor", str(winsor_value)])
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
        '-i', '--input_file', '--input_csv',
        dest='input_file',
        type=str,
        default=str(DEFAULT_INPUT_FILE),
        help="Archivo CSV de entrada (default: %(default)s)"
    )
    parser.add_argument(
        '-o', '--output_root',
        type=str,
        default=None,
        help="Directorio donde se crearán las carpetas 'results_*' (default: carpeta del script)"
    )
    parser.add_argument(
        '--python_cmd',
        type=str,
        default="python3",
        help="Comando para ejecutar Python (ej. 'python3' en Mac, 'python' en Win)\n(default: %(default)s)"
    )
    parser.add_argument(
        '--digest_fingerprint_tag',
        type=str,
        default=None,
        help="Etiqueta base (wXXX_pYYY) a usar para el digest/tablas (default: primer combo)"
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

    # Normalizar rutas
    input_file_path = Path(args.input_file).expanduser()
    if not input_file_path.is_absolute():
        input_file_path = (Path.cwd() / input_file_path).resolve()
    args.input_file = str(input_file_path)

    if args.output_root:
        output_root = Path(args.output_root).expanduser()
        if not output_root.is_absolute():
            output_root = (Path.cwd() / output_root).resolve()
    else:
        output_root = SCRIPT_DIR
    output_root.mkdir(parents=True, exist_ok=True)

    log_path = output_root / f"run_master_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    with tee_output(log_path):
        print(f"--- Log registrado en: {log_path} ---")

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
        print(f"Directorio base de resultados: {output_root}")
        print(f"Valores Winsor:   {winsor_list}")
        print(f"Valores Prime:    {prime_list}")
        print(f"Valores Jitter:   {jitter_list}")
        print(f"Jitter N-Runs:    {n_runs}")
        print("--------------------------------------------")
        fingerprint_tag = args.digest_fingerprint_tag or f"w{sanitize_tag(winsor_list[0])}_p{sanitize_tag(prime_list[0])}"
        digest_label = infer_digest_label(args.input_file)

        # --- FASE 1: ROBUSTEZ (Workflow Bifurcado) ---
        if not args.skip_robustness:
            print("\n\n--- FASE 1: INICIANDO ANÁLISIS DE ROBUSTEZ (Workflow Completo) ---")
            robustness_jobs = list(product(winsor_list, prime_list))
            print(f"Se ejecutarán {len(robustness_jobs)} pipelines de robustez (de 5 pasos c/u)...")
            
            for w_val, p_val in robustness_jobs:
                try:
                    run_robustness_pipeline(w_val, p_val, args.input_file, args.python_cmd, output_root)
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
            config_path = SCRIPT_DIR / "doft_config.json"
            
            if not config_path.exists() or args.skip_robustness:
                try:
                    print("El 'doft_config.json' no existe o se saltó la Fase 1. Generando config base...")
                    run_calibration_only(baseline_w, baseline_p, args.input_file, args.python_cmd, output_root)
                except Exception as e:
                    print("¡¡ERROR FATAL!! No se pudo generar el 'doft_config.json' base.")
                    print(f"Detalle del error: {e}")
                    sys.exit(1)
            else:
                print(f"--- 'doft_config.json' ya existe. Usando como base para Jitter. ---")
                
            sensitivity_jobs = list(product(jitter_list, prime_list))
            print(f"Se ejecutarán {len(sensitivity_jobs)} pipelines de sensibilidad...")
            sensitivity_dir = output_root / "sensitivity"
            sensitivity_dir.mkdir(parents=True, exist_ok=True)

            for j_val, p_val in sensitivity_jobs:
                j_tag = str(j_val).replace(".", "p").replace("-", "neg")
                p_tag = str(int(p_val))
                tag = f"j{j_tag}_p{p_tag}"
                summary_file = sensitivity_dir / f"sensitivity_{tag}.csv"
                if summary_file.exists():
                    if sensitivity_summary_complete(summary_file):
                        print(f"--- Resultados de sensibilidad para {tag} detectados. Saltando. ---")
                        continue
                    elif upgrade_sensitivity_summary(summary_file, baseline_w, p_val, j_val):
                        print(f"--- Resumen de sensibilidad para {tag} actualizado. ---")
                        continue
                    else:
                        print(f"--- Resumen de sensibilidad obsoleto para {tag}. Regenerando. ---")
                try:
                    run_sensitivity_analysis(j_val, p_val, args.input_file, n_runs, args.python_cmd, sensitivity_dir, baseline_w)
                except Exception as e:
                    print(f"¡¡ERROR INESPERADO!! El pipeline (J={j_val}, P={p_val}) falló.")
                    print(f"Detalle del error: {e}")
                    
        elif args.skip_sensitivity:
            print("\n--- FASE 2: Saltando Análisis de Sensibilidad (por flag). ---")
        else:
            print("\n--- FASE 2: No se especificaron --jitter_values, saltando análisis de sensibilidad. ---")

        digest_dir = output_root / "digest"
        digest_targets = [
            digest_dir / f"eta_gamma_summary_{digest_label}.csv",
            digest_dir / f"integer_fingerprint_summary_{digest_label}.csv",
            digest_dir / f"rational_q_summary_{digest_label}.csv",
            digest_dir / f"kappa_impact_{digest_label}.csv",
        ]
        result_dirs = [d for d in output_root.iterdir() if d.is_dir() and d.name.startswith("results_w")]
        if not result_dirs:
            print("\n--- No se encontraron resultados 'results_w*'. Saltando digest. ---")
        elif digest_dir.exists() and all(path.exists() for path in digest_targets):
            print("\n--- Digest detectado previamente. Saltando generación resumida. ---")
        else:
            cmd_digest = [
                args.python_cmd, str(RUN_DIGEST_SCRIPT),
                "--output_root", str(output_root),
                "--fingerprint_tag", fingerprint_tag,
                "--label", digest_label
            ]
            execute_command(cmd_digest, "Digest resumido (archivos maestros)")
            
        print("\n\n============================================")
        print("   ¡Master Pipeline finalizado!   ")
        print("============================================")


if __name__ == "__main__":
    main()
def upgrade_sensitivity_summary(path: Path, winsor_value: float, prime_val: int, jitter_val: float) -> bool:
    try:
        df = pd.read_csv(path)
    except Exception:
        return False
    changed = False
    if "prime_max_val" not in df.columns:
        df["prime_max_val"] = prime_val
        changed = True
    if winsor_value is not None and "winsor" not in df.columns:
        df["winsor"] = winsor_value
        changed = True
    if "jitter_percent" not in df.columns:
        df["jitter_percent"] = jitter_val
        changed = True
    if changed:
        df.to_csv(path, index=False)
    return "prime_max_val" in df.columns and "winsor" in df.columns
