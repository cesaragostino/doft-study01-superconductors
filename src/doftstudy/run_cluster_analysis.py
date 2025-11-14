import os, math, argparse, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Importar dependencias
try:
    from scipy.optimize import lsq_linear
except ImportError:
    print("ERROR: Faltan dependencias. Ejecuta: pip install scipy")
    exit()

# --- Constantes físicas ---
kB = 1.380649e-23; h = 6.62607015e-34; eV = 1.602176634e-19
CONFIG_FILE = "doft_config.json"

# --- Funciones de conversión ---
def f_from_K(T_K): 
    try: T = float(T_K); return (kB*T)/h if T>0 else np.nan
    except Exception: return np.nan

def f_from_eV(E_eV): 
    try: E = float(E_eV); return (E*eV)/h if E>0 else np.nan
    except Exception: return np.nan

# --- Funciones de Locking (Duplicadas) ---
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

SMALL_DENOMINATORS = np.array([1, 2, 3, 4, 5, 6, 7, 8])
PRIME_VALS = None

# --- ¡CAMBIO! La función ahora acepta max_val ---
def generate_prime_products(max_val=10000):
    global PRIME_VALS
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
    if not np.isfinite(R): return (np.nan, np.nan, lock_family)
    
    # --- Assert de coherencia DOFT (automático en el pipeline) ---
    if category in ['SC_TypeI', 'SC_TypeII'] and lock_family == 'integer':
        return best_integer_lock(R)
    if category in ['Superfluid', 'SC_Molecular'] and lock_family == 'rational':
        return best_rational_lock(R)
        
    if lock_family == 'integer':
        return best_integer_lock(R)
    if lock_family == 'rational':
        return best_rational_lock(R)
    if lock_family == 'mixed':
        if R < 1.0:
            err_r, val_r, _ = best_rational_lock(R)
            err_i, val_i, _ = best_integer_lock(R)
        else:
            err_i, val_i, _ = best_integer_lock(R)
            err_r, val_r, _ = best_rational_lock(R)
        if np.isnan(err_i) and np.isnan(err_r): return (np.nan, np.nan, "mixed")
        if np.isnan(err_r): return (err_i, val_i, "integer")
        if np.isnan(err_i): return (err_r, val_r, "rational")
        
        if (err_i <= err_r):
            if val_i not in [12, 24, 30, 4, 6, 8, 10, 14, 18, 21, 28, 42, 60, 84, 105, 210, 2, 3, 5, 7] and val_i % 10 != 0:
                pass
            return (err_i, val_i, "integer")
        else:
            return (err_r, val_r, "rational")
            
    return (np.nan, np.nan, lock_family)

# --- Función de Ajuste de Kappa (Limpia) ---
def fit_kappa_global(df_fit, C_AB):
    """Ajusta un kappa global (>=0) para un material de cluster."""
    if not np.isfinite(C_AB) or C_AB == 0 or df_fit.empty:
        return 0.0 # No se puede ajustar kappa
    
    try:
        fit_data = df_fit.dropna(subset=['R_corr_eta', 'prime_value', 'R_obs'])
        if fit_data.empty:
            return 0.0
        T_i = fit_data['R_corr_eta'].values
        pv_i = fit_data['prime_value'].values
        R_obs_i = fit_data['R_obs'].values
        
        A_i = T_i - pv_i
        B_i = C_AB * R_obs_i
        sum_A_B = np.sum(A_i * B_i)
        sum_B_B = np.sum(B_i**2)
        
        if sum_B_B == 0: kappa_unconstrained = 0.0
        else: kappa_unconstrained = sum_A_B / sum_B_B
        
        kappa_global = max(0.0, kappa_unconstrained)
        return kappa_global
        
    except Exception as e:
        print(f"ERROR: Falló el ajuste de Kappa. {e}")
        return 0.0


# --- ¡CAMBIO! Se añadieron 'prime_max_val' y 'jitter_percent' ---
def analyze_clusters(df, outdir, run_label, g_univ, e_univ, estimate_kappa, prime_max_val=10000, jitter_percent=0.0):
    results = []
    anclas = {}
    
    generate_prime_products(max_val=prime_max_val)
    
    print(f"\n--- Aplicando Corrección Universal DOFT ---")
    print(f"--- Cargando Parámetros Calibrados desde {CONFIG_FILE}:")
    print(f"--- Gamma (g) = {g_univ:.6g}")
    print(f"--- Eta (e)   = {e_univ:.6e}")
    print(f"--- Estimando Kappa (k)? = {estimate_kappa}")
    print(f"--- Jitter de Anclas   = {jitter_percent:.2f}%")
    print(f"--- Prime Max Val      = {prime_max_val}")
    print("--------------------------------------------------")

    print("\n--- 1. Construyendo Anclas 'Límpias' y Parámetro 'X' (solo clusters) ---")
    
    df_clusters = df[df['sub_network'].str.lower() != 'single'].copy()

    if jitter_percent > 0.0 and (not df_clusters.empty):
        print(f"--- APLICANDO JITTER DE +/- {jitter_percent:.1f}% A LAS ANCLAS ---")
        jitter_factor = jitter_percent / 100.0
        anchor_cols = ['Tc_K', 'Gap_meV', 'ThetaD_K', 'EF_eV']
        for col in anchor_cols:
            if col in df_clusters.columns:
                random_mult = 1.0 + np.random.uniform(-jitter_factor, jitter_factor, size=len(df_clusters))
                df_clusters.loc[:, col] = df_clusters[col] * random_mult

    for _, row in df_clusters.iterrows():
        try:
            name = str(row['name']); sub = str(row['sub_network'])
            
            f_th = f_from_K(row['Tc_K'])
            f_gap = f_from_eV(row['Gap_meV'] * 1e-3)
            f_debye = f_from_K(row['ThetaD_K'])
            f_elec = f_from_eV(row['EF_eV'])
            
            if not (np.isfinite(f_th) and pd.notnull(row['Gap_meV']) and np.isfinite(f_debye)):
                pass
                
            Tc_val = float(row['Tc_K']); ThD_val = float(row['ThetaD_K'])
            X_c = ThD_val / Tc_val if (Tc_val > 0 and ThD_val > 0) else np.nan
            sub_order = int(row['sub_order']) if pd.notnull(row['sub_order']) else np.nan
                
            anclas[(name, sub)] = {
                'f_th': f_th, 'f_gap': f_gap, 'f_debye': f_debye, 'f_elec': f_elec,
                'X': X_c, 'lock_family': str(row['lock_family']).lower(),
                'category': str(row['category']),
                'sub_order': sub_order
            }
        
        except Exception as e:
            print(f"Error procesando fila: {row}. Error: {e}")

    if jitter_percent == 0.0:
        for (name, sub), ancla in anclas.items():
            print(f"Ancla creada para: {name}-{sub} (X = {ancla['X']:.2f}, sub_order = {ancla['sub_order']})")


    print("\n--- 2. Calculando Saltos INTRA e INTER-Canal ---")
    
    all_materials_data = []
    material_names = df_clusters['name'].unique()
    
    for name in material_names:
        material_records = []
        subs = [k[1] for k in anclas.keys() if k[0] == name]
        if not subs: continue

        if jitter_percent == 0.0:
            print(f"Procesando material: {name}")

        C_AB = 0.0
        
        if name == 'LaH10':
            sub_A_order, sub_B_order = 2, 3
        else:
            sub_A_order, sub_B_order = 1, 2

        sub_A_X, sub_B_X = np.nan, np.nan
        sub_A_name, sub_B_name = None, None

        for sub in subs:
            if (name, sub) not in anclas: continue
            if anclas[(name, sub)]['sub_order'] == sub_A_order:
                sub_A_X = anclas[(name, sub)]['X']
                sub_A_name = sub
            elif anclas[(name, sub)]['sub_order'] == sub_B_order:
                sub_B_X = anclas[(name, sub)]['X']
                sub_B_name = sub
        
        if np.isfinite(sub_A_X) and np.isfinite(sub_B_X):
            C_AB = abs(sub_A_X - sub_B_X)
            if jitter_percent == 0.0:
                print(f"Parámetro de Contraste C_AB (orden {sub_A_order} vs {sub_B_order}) para {name} = {C_AB:.4f}")
        else:
            if jitter_percent == 0.0:
                print(f"ADVERTENCIA: No se encontraron subredes con orden {sub_A_order} y {sub_B_order} para {name}. C_AB será 0.")
        
        
        for sub in subs:
            if (name, sub) not in anclas: continue
            
            freqs = anclas[(name, sub)]
            X = freqs['X']; lock_family = freqs['lock_family']; cat = freqs['category']
            
            anchors_list = [
                ("thermal@Tc", freqs['f_th']), ("gap(Δ)", freqs['f_gap']),
                ("Debye", freqs['f_debye']), ("electronic(E_F)", freqs['f_elec']),
            ]
            jump_defs = [
                (0,1,"thermal@Tc→gap(Δ)", 3), (1,2,"gap(Δ)→Debye", 2),
                (2,3,"Debye→electronic(E_F)", 1),
            ]

            for i,j,label,d in jump_defs:
                fa = anchors_list[i][1]; fb = anchors_list[j][1]
                
                if not (np.isfinite(fa) and np.isfinite(fb) and fa > 0):
                    material_records.append({
                        "name": name, "sub_network": sub, "category": cat,
                        "jump_type": "intra", "jump_desc": label, "d": d, "X": X, "C_AB": C_AB,
                        "lock_family": lock_family, "R_obs": np.nan, "chosen_lock": "skip (nan)", 
                        "prime_value": np.nan, "err_before": np.nan, 
                        "R_corr_eta": np.nan, "err_after_eta": np.nan,
                        "R_corr_kappa": np.nan, "err_after_kappa": np.nan
                    })
                    continue 

                R_obs = fb/fa
                pe, pv, chosen_lock = find_best_lock(R_obs, lock_family, cat)
                
                correction_factor = max(0.2, 1 - g_univ * (X**2) - e_univ * d * X) if np.isfinite(X) else 1.0
                R_corr_eta = R_obs * correction_factor
                err_after_eta = abs(R_corr_eta - pv) / R_corr_eta if (R_corr_eta != 0 and np.isfinite(pv)) else np.nan
                
                material_records.append({
                    "name": name, "sub_network": sub, "category": cat,
                    "jump_type": "intra", "jump_desc": label, "d": d, "X": X, "C_AB": C_AB,
                    "lock_family": lock_family, "R_obs": R_obs, "chosen_lock": chosen_lock, 
                    "prime_value": pv, "err_before": pe, 
                    "R_corr_eta": R_corr_eta, "err_after_eta": err_after_eta,
                    "R_corr_kappa": np.nan, "err_after_kappa": np.nan
                })

        if sub_A_name and sub_B_name:
            if (name, sub_A_name) not in anclas or (name, sub_B_name) not in anclas:
                continue
                
            ancla_A = anclas[(name, sub_A_name)]
            ancla_B = anclas[(name, sub_B_name)]
            
            if jitter_percent == 0.0:
                print(f"Calculando saltos INTER-canal para: {name} ({sub_A_name} vs {sub_B_name})")
            
            nan_cols = {"d": np.nan, "X": np.nan, "C_AB": np.nan, "R_corr_eta": np.nan, "err_after_eta": np.nan, "R_corr_kappa": np.nan, "err_after_kappa": np.nan}
            
            R_gap_ratio = ancla_A['f_gap'] / ancla_B['f_gap']
            pe_gap, pv_gap, lock_gap = best_rational_lock(R_gap_ratio)
            material_records.append({
                "name": name, "sub_network": f"{sub_A_name}-vs-{sub_B_name}", "category": ancla_A['category'],
                "jump_type": "inter", "jump_desc": f"gap({sub_A_name}) / gap({sub_B_name})",
                "lock_family": "rational", "R_obs": R_gap_ratio, "chosen_lock": lock_gap,
                "prime_value": pv_gap, "err_before": pe_gap, **nan_cols
            })

            R_deb_ratio = ancla_A['f_debye'] / ancla_B['f_debye']
            pe_deb, pv_deb, lock_deb = best_rational_lock(R_deb_ratio)
            material_records.append({
                "name": name, "sub_network": f"{sub_A_name}-vs-{sub_B_name}", "category": ancla_A['category'],
                "jump_type": "inter", "jump_desc": f"Debye({sub_A_name}) / Debye({sub_B_name})",
                "lock_family": "rational", "R_obs": R_deb_ratio, "chosen_lock": lock_deb,
                "prime_value": pv_deb, "err_before": pe_deb, **nan_cols
            })
        else:
            if jitter_percent == 0.0:
                print(f"ADVERTENCIA: No se pudieron encontrar pares de subredes (orden {sub_A_order} y {sub_B_order}) para {name}. Omitiendo saltos INTER.")


        df_material = pd.DataFrame(material_records)
        kappa_global = 0.0 
        
        if estimate_kappa:
            if jitter_percent == 0.0:
                print(f"--- 3. Ajustando Coeficiente de Cluster (kappa >= 0) para {name} ---")
            
            df_fit = df_material[df_material['jump_type'] == 'intra']
            kappa_global = fit_kappa_global(df_fit, C_AB)
            
            print(f"--- Coeficiente de Cluster (kappa) para {name} = {kappa_global:.6g} ---")
        
        T_i = df_material['R_corr_eta']
        R_obs_i = df_material['R_obs']
        pv_i = df_material['prime_value']
        
        R_corr_kappa = T_i - kappa_global * C_AB * R_obs_i
        err_after_kappa = np.abs(R_corr_kappa - pv_i) / np.abs(R_corr_kappa)
        
        intra_mask = df_material['jump_type'] == 'intra'
        df_material.loc[intra_mask, 'R_corr_kappa'] = R_corr_kappa[intra_mask]
        df_material.loc[intra_mask, 'err_after_kappa'] = err_after_kappa[intra_mask]
        
        inter_mask = df_material['jump_type'] == 'inter'
        df_material.loc[inter_mask, 'R_corr_eta'] = df_material.loc[inter_mask, 'R_obs']
        df_material.loc[inter_mask, 'err_after_eta'] = df_material.loc[inter_mask, 'err_before']
        df_material.loc[inter_mask, 'R_corr_kappa'] = df_material.loc[inter_mask, 'R_obs']
        df_material.loc[inter_mask, 'err_after_kappa'] = df_material.loc[inter_mask, 'err_before']
        
        all_materials_data.append(df_material)

    if jitter_percent == 0.0:
        print("\n--- 4. Guardando Reporte de Diagnóstico de Cluster ---")
    
    if not all_materials_data:
        print("ADVERTENCIA: No se generaron resultados (¿CSV de entrada vacío o sin sub-redes?)")
        return

    out_df = pd.concat(all_materials_data, ignore_index=True)
    
    # --- ¡LA LÍNEA CORREGIDA! ---
    cols_order = [
        "name", "sub_network", "category", "lock_family", "jump_type", "jump_desc", "d", "X", "C_AB",
        "R_obs", "chosen_lock", "prime_value", "err_before", 
        "R_corr_eta", "err_after_eta", 
        "R_corr_kappa", "err_after_kappa"
    ]
    # --- FIN DE LA CORRECCIÓN ---
    
    final_cols = [c for c in cols_order if c in out_df.columns]
    out_df = out_df[final_cols]
    
    if outdir:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        out_path = os.path.join(outdir, f"results_{run_label}.csv")
        out_df.to_csv(out_path, index=False, float_format='%.6g')
        if jitter_percent == 0.0:
            print(f"--- ¡Reporte de Diagnóstico de Cluster guardado en: {out_path} ---")

    return out_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DOFT Study - Cluster Analyzer (v-Master 2.4)")
    parser.add_argument("-i", "--input", dest="input_csv", type=str, required=True, help="Input CSV file (e.g., materials_clusters_real_v2.csv)")
    parser.add_argument("-o", "--outdir", dest="outdir", type=str, required=True, help="Output directory (e.g., 'results_clusters')")
    parser.add_argument("-l", "--label", dest="run_label", type=str, default="cluster_analysis", help="Run label for output files")
    parser.add_argument("--estimate_kappa", action='store_true', help="Activa el ajuste de kappa (default: False)")
    
    parser.add_argument("--prime_max_val", dest="prime_max_val", type=int, default=10000, help="Valor máximo para la generación de 'prime values' (default: 10000)")
    parser.add_argument("--jitter_percent", dest="jitter_percent", type=float, default=0.0, help="Porcentaje de jitter +/- aleatorio para aplicar a las anclas (default: 0.0)")

    args = parser.parse_args()
    
    try:
        with open(CONFIG_FILE, 'r') as f:
            config = json.load(f)
        CALIBRATED_GAMMA = config.get('CALIBRATED_GAMMA', 0.0)
        CALIBRATED_ETA = config.get('CALIBRATED_ETA', 0.0)
        
        config_prime_max = config.get('METADATA', {}).get('prime_max_val', 10000)
        prime_max_val_to_use = args.prime_max_val
        if args.prime_max_val == 10000 and config_prime_max != 10000:
            print(f"--- Usando 'prime_max_val' del archivo de config: {config_prime_max} ---")
            prime_max_val_to_use = config_prime_max
        
        if CALIBRATED_ETA == 0.0 and args.jitter_percent == 0.0:
            print(f"ADVERTENCIA: {CONFIG_FILE} contiene Eta=0.0. ¿Ejecutaste 'run_calibration.py' primero?")
    except FileNotFoundError:
        print(f"ERROR: No se pudo encontrar el archivo de configuración: {CONFIG_FILE}")
        print("Por favor, ejecuta 'run_calibration.py' primero para generar este archivo.")
        exit()
    except Exception as e:
        print(f"ERROR: No se pudo leer {CONFIG_FILE}. Error: {e}"); exit()

    try:
        df = pd.read_csv(args.input_csv)
        
        analyze_clusters(df, args.outdir, args.run_label, 
                         CALIBRATED_GAMMA, CALIBRATED_ETA, 
                         args.estimate_kappa,
                         prime_max_val_to_use,
                         args.jitter_percent)
            
    except FileNotFoundError:
        print(f"ERROR: No se pudo encontrar el archivo de entrada: {args.input_csv}")
    except Exception as e:
        print(f"Ha ocurrido un error inesperado: {e}")