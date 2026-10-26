# -*- coding: utf-8 -*-
"""
Procesamiento_DWI.py

Pipeline de difusion del Lab. de Visualizacion UdeC: preprocesamiento, tractografia, segmentacion de
fasciculos (DWM/SWM), metricas DTI por fasciculo y QC. Sirve para datos investigativos (b0 de fase inversa:
topup) y clinicos (solo AP o solo PA: Synb0-DisCo); el tipo se detecta por sujeto. Ver README.md.

Los atlas, las plantillas MNI, main_index, tck2bundles.npy y read_write_bundle.py deben estar en la misma
carpeta que este script (y que qc_dwi.py).

@author: snava
"""

import argparse
import glob
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
import traceback

import nibabel as nb
import nibabel.streamlines.tck as TF
import numpy as np
from nibabel.streamlines.tractogram import Tractogram

# Atlas, plantillas y main_index: en la carpeta de este script, como en la v1
RECURSOS = os.path.dirname(os.path.abspath(__file__))
MNI_1mm_brain = os.path.join(RECURSOS, "MNI152_T1_1mm_brain.nii.gz")
ANTS_PREFIX = "ants_T12MNI_"  # salidas de antsRegistrationSyN.sh (0GenericAffine.mat, 1Warp, 1InverseWarp)
MAIN_INDEX = os.path.join(RECURSOS, "main_index")
WORKDIR_NAME = "proc_dwi_v2"
METODOS_FILE = "metodos.json"  # entradas, correccion de distorsion, parametros y versiones del sujeto
PE_FILE = "pe.json"            # fase y readout de la DWI, en los ejes del NIfTI de entrada (para MRtrix)
SYNB0_IMAGE = "leonyichencai/synb0-disco:v3.1"
ATLAS = {"SWM": ("AtlasRo", "AtlasRo_estables.txt", ""), "DWM": ("atlas_faisceaux", "atlas_faisceaux.txt", "_MNI")}
PE_VALIDAS = ("i", "i-", "j", "j-", "k", "k-")
T1_BRAIN = "T1w_acpc_dc_restore_brain.nii.gz"
T1_HEAD = "T1w_acpc_dc_restore.nii.gz"
TCK_DWI = "tractography_prob_sift_3M.tck"         # tractografia SIFT original (espacio DWI, paso ~1 mm)
TCK_21P = "tractography_prob_sift_3M_21p.tck"     # remuestreada a 21 puntos, solo para llevarla a MNI
FINAL_FILE = "tractography_prob_sift_3M_21p_MNI.tck"
N_PUNTOS = 21  # puntos por fibra que espera main_index
# Mapas que se promedian por fasciculo (calculados con -dti, en espacio DWI)
METRICAS = [("FA", "dti_FA.nii.gz"), ("AD", "dti_AD.nii.gz"), ("RD", "dti_RD.nii.gz"), ("MD", "dti_MD.nii.gz")]

CLEANUP_FILE = ".cleanup_done"

# Intermedios que borra -cleanup. Quedan dwi_preproc_unbiased, mascara, DTI, registros, 5TT, FOD
# normalizados, tractogramas, warps, DWM/SWM, metricas, el PDF de QC y los logs.
LIMPIEZA = [
    "run[0-9].nii.gz", "run[0-9].bvec", "run[0-9].bval",      # copias de las series de entrada
    "run*_den.nii.gz", "run*_den_unr.nii.gz",                # quedan incluidos en dwi_den
    "run*_noise.nii.gz",                                      # mapas de ruido de dwidenoise
    "b0_PA.nii.gz", "b0_PA_unr.nii.gz", "b0_AP.nii.gz", "b0_PA_sel.nii.gz", "b0_pair.nii.gz",
    "bzero.nii.gz",                                           # b0 medio de dwi_preproc_unbiased
    "dwi_den.nii.gz",                                         # entrada de eddy
    "dwi_preproc.nii.gz",                                     # = dwi_preproc_unbiased * bias_field
    TCK_21P,                                                  # se conservan TCK_DWI y FINAL_FILE
    "tractography_prob_sift_3M_21p_MNI.bundles",              # tractograma completo en .bundles
    "tractography_prob_sift_3M_21p_MNI.bundlesdata",          # (se conserva el .tck)
    "wmfod.nii.gz", "gm.nii.gz", "csf.nii.gz",                # reemplazados por *_norm
    "tensor.nii.gz", "RF_voxels.nii.gz", "b0_mean_hdbet.nii.gz", "b0_hdbet_bet.nii.gz",
    "diff2struct_epireg*", "diff2struct_init.mat", "diff2struct_fast*",  # auxiliares de epi_reg / BBR
    f"{ANTS_PREFIX}Warped.nii.gz", f"{ANTS_PREFIX}InverseWarped.nii.gz",  # T1w deformado (QC lo recalcula)
    "output.tck", "tmp_dwicat_*.mif", "gmwmi_seed.nii.gz",
    "synb0", "b0_pre.nii.gz", "b0_topup*",                    # Synb0 (queda b0_synb0.nii.gz para el QC)
    "eddy_in.nii", "eddy_out*", "eddy_mask.nii.gz", "eddy_acqp.txt", "eddy_index.txt", "bvecs", "bvals",
    "eddy_qc",                                                # metricas de eddy (su reporte queda en QC/)
    "qc_tmp", "qc",                                           # imagenes intermedias del QC (queda QC/)
]


def tck2bundles(tck, salida):
    """Tractografia en MNI (.tck, mm) al formato .bundles que lee main_index (transformacion tck2bundles.npy)."""
    import read_write_bundle as bt
    T = np.load(os.path.join(RECURSOS, "tck2bundles.npy"))
    sl = nb.streamlines.load(tck).streamlines
    datos = sl.get_data()
    P = np.empty_like(datos, dtype=np.float32)
    for i in range(0, len(datos), 10_000_000):  # por bloques para no duplicar la memoria
        P[i:i + 10_000_000] = datos[i:i + 10_000_000] @ T[:3, :3].T + T[:3, 3]
    bt.write_bundle(salida, np.split(P, np.cumsum(sl._lengths)[:-1]))


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

LOG_FILE = None  # pipeline_v2.log del sujeto en proceso


def log_msg(texto):
    linea = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {texto}"
    print(linea, flush=True)
    if LOG_FILE:
        with open(LOG_FILE, "a") as f:
            f.write(linea + "\n")


def run(cmd, cwd):
    """Ejecuta un comando de shell en cwd y detiene el pipeline si falla."""
    log_msg(f"$ {cmd}")
    shell_cmd = cmd
    if LOG_FILE:
        # el comando (con fecha y hora) y su salida quedan en pipeline_v2.log
        shell_cmd = f"set -o pipefail; ( {cmd} ) 2>&1 | tee -a '{LOG_FILE}'"
    try:
        subprocess.run(shell_cmd, shell=True, check=True, cwd=cwd, executable="/bin/bash")
    except subprocess.CalledProcessError:
        log_msg(f"ERROR en: {cmd}")
        raise


def step(outputs, cmds, cwd):
    """Ejecuta cmds solo si alguna de las salidas no existe (permite reanudar)."""
    outputs = [outputs] if isinstance(outputs, str) else outputs
    if all(existe(cwd, o) for o in outputs):
        print(f"  [omitido, ya existe] {', '.join(outputs)}")
        return
    for c in ([cmds] if isinstance(cmds, str) else cmds):
        run(c, cwd)


def existe(wd, archivo):
    return os.path.exists(os.path.join(wd, archivo))


def limpiar(wd):
    """Borra los intermedios de LIMPIEZA y deja la lista en .cleanup_done (el sujeto queda terminado)."""
    borrados, liberado = [], 0
    for patron in LIMPIEZA:
        for p in sorted(glob.glob(os.path.join(wd, patron))):
            if os.path.isdir(p):
                for raiz, _, archivos in os.walk(p):
                    liberado += sum(os.path.getsize(os.path.join(raiz, a)) for a in archivos)
                shutil.rmtree(p)
            else:
                liberado += os.path.getsize(p)
                os.remove(p)
            borrados.append(os.path.relpath(p, wd))
    with open(os.path.join(wd, CLEANUP_FILE), "a") as f:
        f.writelines(b + "\n" for b in borrados)
    log_msg(f"$ __cleanup__ {len(borrados)} archivos borrados, {liberado / 1e9:.2f} GB liberados")
    print(f"  Cleanup: {len(borrados)} archivos borrados, {liberado / 1e9:.2f} GB liberados")


def contar_puntos(tck):
    """Numero de fibras por cantidad de puntos (main_index requiere que todas tengan N_PUNTOS)."""
    n = {}
    for f in nb.streamlines.load(tck, lazy_load=True).streamlines:
        n[len(f)] = n.get(len(f), 0) + 1
    return n


def puntos_densos(fibras, paso):
    """Puntos de las fibras remuestreadas cada `paso` mm como maximo, con el indice de fibra de cada uno."""
    datos = np.asarray(fibras.get_data(), dtype=np.float32)
    largos = np.asarray(fibras._lengths)
    id_pt = np.repeat(np.arange(len(largos)), largos)
    fin = np.cumsum(largos) - 1                              # ultimo punto de cada fibra
    ini = np.setdiff1d(np.arange(len(datos) - 1), fin)        # segmentos dentro de una misma fibra
    seg = datos[ini + 1] - datos[ini]
    k = max(1, int(np.ceil(np.linalg.norm(seg, axis=1).max() / paso)))
    pts = [datos] + [datos[ini] + t * seg for t in np.arange(1, k) / k]
    ids = [id_pt] + [id_pt[ini]] * (k - 1)
    return np.concatenate(pts), np.concatenate(ids)


def metricas_fasciculos(wd, tipo, atlas_txt, sufijo, salida):
    """
    Para cada fasciculo del atlas y cada mapa de METRICAS (espacio DWI):
      - n_streamlines: fibras del fasciculo.
      - vol_mm3: volumen de los voxeles que atraviesan las fibras, dentro de la mascara DWI.
      - mean/std/median/IQR: se promedia el mapa a lo largo de cada fibra (interpolacion trilineal, puntos
        cada <= medio voxel) y se resumen esos promedios entre las fibras.
    """
    import pandas as pd
    from scipy.ndimage import map_coordinates
    mask_img = nb.load(os.path.join(wd, "dwi_mask.nii.gz"))
    mask = mask_img.get_fdata() > 0
    inv_aff = np.linalg.inv(mask_img.affine)
    zooms = mask_img.header.get_zooms()[:3]
    mapas = {m: nb.load(os.path.join(wd, f)).get_fdata(dtype=np.float32) for m, f in METRICAS}
    filas = []
    nombres = [l.split()[0].removesuffix(sufijo) for l in open(os.path.join(RECURSOS, atlas_txt)) if l.strip()]
    for nombre in nombres:
        tck = os.path.join(wd, tipo, "results_folder_DWI", nombre + ".tck")
        fibras = nb.streamlines.load(tck).streamlines if os.path.exists(tck) else []
        n, vol, prom = len(fibras), 0.0, {m: np.array([]) for m in mapas}
        if n:
            pts, ids = puntos_densos(fibras, 0.5 * min(zooms))
            ijk = pts @ inv_aff[:3, :3].T + inv_aff[:3, 3]
            v = np.rint(ijk).astype(np.int64)
            v = v[np.all((v >= 0) & (v < mask.shape), axis=1)]
            vox = np.unique(np.ravel_multi_index(v.T, mask.shape))
            vol = float(mask.ravel()[vox].sum() * np.prod(zooms))
            cuenta = np.bincount(ids, minlength=n)
            for m, mapa in mapas.items():
                val = map_coordinates(mapa, ijk.T, order=1, mode="nearest")
                prom[m] = np.bincount(ids, weights=val, minlength=n) / cuenta
        for m in mapas:
            x = prom[m]
            q1, med, q3 = np.percentile(x, [25, 50, 75]) if len(x) else (np.nan,) * 3
            filas.append({"bundle": nombre, "metrics": m, "space": "DWI", "n_streamlines": n, "vol_mm3": vol,
                          "mean": x.mean() if len(x) else np.nan, "std": x.std(ddof=1) if len(x) > 1 else np.nan,
                          "median": med, "IQR": q3 - q1})
    os.makedirs(os.path.dirname(salida), exist_ok=True)
    pd.DataFrame(filas).round(6).to_csv(salida, index=False)


def metodo_mascara(pedido):
    """Resuelve -mask: 'auto' usa hd-bet si esta instalado y si no dwi2mask."""
    if pedido == "auto":
        return "hdbet" if shutil.which("hd-bet") else "dwi2mask"
    if pedido == "hdbet" and not shutil.which("hd-bet"):
        raise SystemExit("Error: -mask hdbet pero hd-bet no esta instalado (pip install hd-bet).")
    return pedido


def hdbet_gpu():
    """True si el PyTorch de hd-bet puede usar la GPU (no puede, por ejemplo, con drivers de NVIDIA antiguos)."""
    if not shutil.which("nvidia-smi"):
        return False
    with open(shutil.which("hd-bet")) as f:
        primera = f.readline()
    python = shlex.split(primera[2:]) if primera.startswith("#!") else [sys.executable]
    try:
        return subprocess.run(python + ["-c", "import torch; torch.zeros(1).cuda()"],
                              capture_output=True, timeout=300,
                              env=dict(os.environ, MKL_THREADING_LAYER="GNU")).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def opposite_pe(a, b):
    return a.rstrip("-") == b.rstrip("-") and a.endswith("-") != b.endswith("-")


def find_dwi_series(suj_dir):
    """
    Busca recursivamente (excluyendo la carpeta de trabajo) las series de difusion ORIGINAL convertidas con
    dcm2niix. Las que tienen bval/bvec con b > 50 son runs DWI; las de solo b0 son candidatas a b0 de fase
    inversa. Las series derivadas del escaner (ADC, TRACEW, FA, ColFA) se descartan por ImageType.
    Sin series dcm2niix, usa diff.nii.gz/.bval/.bvec de la carpeta del sujeto (formato v1).
    Devuelve (runs_dwi, series_b0, T1w_brain, T1w_head).
    """
    runs, b0s = [], []
    for js in sorted(glob.glob(os.path.join(suj_dir, "**", "*.json"), recursive=True)):
        if f"/{WORKDIR_NAME}/" in js:
            continue
        base = js[:-5]
        nii = base + ".nii.gz"
        if not os.path.exists(nii):
            continue
        meta = json.load(open(js))
        itype = meta.get("ImageType", [])
        if "DIFFUSION" not in itype or "ORIGINAL" not in itype:
            continue
        serie = {"base": base, "nii": nii, "json": js, "meta": meta, "name": os.path.basename(base),
                 "num": meta.get("SeriesNumber", 0)}
        bval = base + ".bval"
        if os.path.exists(bval) and os.path.exists(base + ".bvec") and np.loadtxt(bval).max() > 50:
            runs.append(serie)
        else:
            b0s.append(serie)

    base = os.path.join(suj_dir, "diff")
    if not runs and all(os.path.exists(base + e) for e in (".nii.gz", ".bval", ".bvec")):
        js = base + ".json"
        runs = [{"base": base, "nii": base + ".nii.gz", "json": js if os.path.exists(js) else None,
                 "meta": json.load(open(js)) if os.path.exists(js) else {}, "name": "diff", "num": 0}]

    runs = sorted(runs, key=lambda s: s["num"])
    b0s = sorted(b0s, key=lambda s: s["num"])
    t1_brain = [p for p in glob.glob(os.path.join(suj_dir, "**", T1_BRAIN), recursive=True)
                if f"/{WORKDIR_NAME}/" not in p]
    t1_head = [p for p in glob.glob(os.path.join(suj_dir, "**", T1_HEAD), recursive=True)
               if f"/{WORKDIR_NAME}/" not in p]
    return runs, b0s, t1_brain, t1_head


def detectar_docker():
    """'docker' si esta instalado, con permisos y el daemon activo (Synb0 corre en Docker)."""
    if shutil.which("docker") and subprocess.run("docker info", shell=True, capture_output=True).returncode == 0:
        return "docker"
    return None


def buscar_licencia(arg):
    """Licencia de FreeSurfer que necesita Synb0."""
    for p in [arg, os.environ.get("FS_LICENSE"),
              os.path.join(os.environ.get("FREESURFER_HOME", ""), "license.txt"),
              os.path.expanduser("~/license.txt")]:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def synb0_cmd(args, licencia, stripped):
    """Comando de Synb0-DisCo; INPUTS/OUTPUTS son synb0/ dentro de la carpeta de trabajo (cwd)."""
    return (f"docker run --rm --user {os.getuid()}:{os.getgid()} -v \"$PWD/synb0/INPUTS\":/INPUTS "
            f"-v \"$PWD/synb0/OUTPUTS\":/OUTPUTS -v '{licencia}':/extra/freesurfer/license.txt "
            f"{args.synb0_image}{' --stripped' if stripped else ''}")


def synb0_version(args):
    r = subprocess.run(f"docker image inspect --format '{{{{index .RepoDigests 0}}}}' {args.synb0_image}",
                       shell=True, capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else args.synb0_image


def recursos_faltantes(segmentacion):
    """Archivos de la carpeta del script que faltan para registrar a MNI y segmentar."""
    req = ["MNI152_T1_1mm_brain.nii.gz", "tck2bundles.npy", "read_write_bundle.py"]
    if segmentacion:
        req.append("main_index")
        for tipo in segmentacion:
            req += ATLAS[tipo][:2]
    return [f for f in req if not os.path.exists(os.path.join(RECURSOS, f))]


def validar_sujeto(suj_dir, args, runtime, licencia):
    """
    Verifica los datos de un sujeto y decide la correccion de distorsion. Devuelve
    (dict_con_rutas, lista_de_errores, lista_de_avisos).
    """
    errores, avisos = [], []
    runs, b0s, t1_brain, t1_head = find_dwi_series(suj_dir)

    if not runs:
        errores.append("no se encontraron series de difusion con .bval/.bvec (dcm2niix) ni diff.nii.gz/.bval/.bvec")
    if len(t1_brain) != 1:
        errores.append(f"se esperaba 1 archivo {T1_BRAIN}, se encontraron {len(t1_brain)} "
                       "(¿se corrio el pipeline estructural HCP?)")
    if not runs:
        return None, errores, avisos

    # Fase y readout: .json de la serie > acq.json del sujeto > argumentos
    ref = runs[0]
    ref_img = nb.load(ref["nii"])
    acq_json = os.path.join(suj_dir, "acq.json")
    acq = json.load(open(acq_json)) if os.path.exists(acq_json) else {}
    pe, trt, fuente = None, None, "argumentos"
    for origen, meta in ((os.path.basename(ref["json"]) if ref["json"] else None, ref["meta"]),
                         ("acq.json", acq), ("argumentos", {"PhaseEncodingDirection": args.pe_dir,
                                                            "TotalReadoutTime": args.readout})):
        if pe is None and meta.get("PhaseEncodingDirection"):
            pe, fuente = meta["PhaseEncodingDirection"], origen
        if trt is None and meta.get("TotalReadoutTime"):
            trt = float(meta["TotalReadoutTime"])
    if pe is not None and pe not in PE_VALIDAS:
        errores.append(f"PhaseEncodingDirection '{pe}' no valida (usar {', '.join(PE_VALIDAS)})")

    for s in runs[1:]:
        img = nb.load(s["nii"])
        if s["meta"].get("PhaseEncodingDirection") != ref["meta"].get("PhaseEncodingDirection"):
            errores.append(f"{s['name']}: PhaseEncodingDirection distinta a {ref['name']}")
        if not np.isclose(s["meta"].get("TotalReadoutTime", -1), ref["meta"].get("TotalReadoutTime", -1)):
            errores.append(f"{s['name']}: TotalReadoutTime distinto a {ref['name']}")
        if img.shape[:3] != ref_img.shape[:3] or not np.allclose(img.affine, ref_img.affine, atol=1e-3):
            errores.append(f"{s['name']}: geometria distinta a {ref['name']}")

    # b0 de fase inversa (adquisicion investigativa): serie de solo b0 con la fase opuesta a la DWI
    rev = [s for s in b0s if pe and opposite_pe(pe, s["meta"].get("PhaseEncodingDirection") or "")]
    faltan_synb0 = [f for f, ok in (("Docker", runtime),
                                    ("licencia de FreeSurfer (-fs_license o FS_LICENSE)", licencia),
                                    ("direccion de fase (.json, acq.json o -pe_dir)", pe)) if not ok]
    sdc = args.sdc
    if sdc == "auto":
        if rev:
            sdc = "topup"
        elif not faltan_synb0:
            sdc = "synb0"
        else:
            sdc = "none"
            avisos.append("sin b0 de fase inversa y Synb0 no disponible (falta: " + ", ".join(faltan_synb0) +
                          "): solo eddy, sin corregir la distorsion por susceptibilidad")
            if pe is None and runtime and licencia:
                avisos.append("para usar Synb0 basta indicar la direccion de fase, p. ej. -pe_dir j- (o un acq.json "
                              "en la carpeta del sujeto)")

    if sdc == "topup":
        if len(rev) != 1:
            errores.append(f"topup requiere 1 serie b0 con la fase opuesta a la DWI ({pe}), se encontraron "
                           f"{len(rev)}: {[s['name'] for s in rev]}")
        else:
            for s in runs + rev:
                for k in ("PhaseEncodingDirection", "TotalReadoutTime"):
                    if s["meta"].get(k) is None:
                        errores.append(f"{s['name']}: falta {k} en el .json")
            inv = rev[0]
            inv_img = nb.load(inv["nii"])
            if trt is not None and inv["meta"].get("TotalReadoutTime") is not None \
                    and not np.isclose(inv["meta"]["TotalReadoutTime"], trt):
                errores.append(f"{inv['name']}: TotalReadoutTime={inv['meta']['TotalReadoutTime']} distinto al de "
                               f"la DWI ({trt})")
            if inv_img.shape[:3] != ref_img.shape[:3] or not np.allclose(
                    inv_img.header.get_zooms()[:3], ref_img.header.get_zooms()[:3], atol=1e-3):
                errores.append(f"{inv['name']}: geometria {inv_img.shape[:3]} / {inv_img.header.get_zooms()[:3]} "
                               f"distinta a la DWI {ref_img.shape[:3]} / {ref_img.header.get_zooms()[:3]}")
    elif sdc == "synb0" and faltan_synb0:
        errores.append("Synb0 no disponible, falta: " + ", ".join(faltan_synb0))
    if sdc != "topup" and rev:
        avisos.append(f"hay una b0 de fase inversa ({rev[0]['name']}) que no se usa con -sdc {sdc}")
    if sdc == "none" and pe is None:
        pe, fuente = "j-", "supuesta"
        avisos.append("direccion de fase desconocida, se asume j- (sin topup solo la usa eddy)")
    if trt is None:
        trt = 0.05
        avisos.append("TotalReadoutTime desconocido, se usa 0.05 s (topup y eddy usan el mismo valor: solo "
                      "escala el campo estimado)")

    if errores:
        return None, errores, avisos
    datos = {"runs": runs, "rev": rev[0] if sdc == "topup" else None, "t1_brain": t1_brain[0],
             "t1_head": t1_head[0] if len(t1_head) == 1 else None, "pe": pe, "trt": trt, "fuente_pe": fuente,
             "sdc": sdc}
    return datos, errores, avisos


def guardar_metodos(wd, **campos):
    """Agrega campos a metodos.json (entradas, correccion de distorsion, parametros y versiones)."""
    p = os.path.join(wd, METODOS_FILE)
    m = json.load(open(p)) if os.path.exists(p) else {}
    m.update(campos)
    with open(p, "w") as f:
        json.dump(m, f, indent=2, ensure_ascii=False)


def versiones():
    def primera_linea(cmd):
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return r.stdout.strip().splitlines()[0] if r.returncode == 0 and r.stdout.strip() else None
    fsl = os.path.join(os.environ.get("FSLDIR", ""), "etc", "fslversion")
    return {"mrtrix": primera_linea("mrconvert -version"),
            "fsl": open(fsl).read().strip() if os.path.exists(fsl) else None,
            "ants": primera_linea("antsRegistration --version")}


def corregir_synb0(wd, datos, args, runtime, licencia, nt):
    """
    Synb0-DisCo (Schilling et al. 2019, 2020): a partir del T1w genera una b0 sin distorsion, que topup usa
    como b0 de fase "inversa" (readout 0) para estimar el campo; eddy corrige luego con ese campo.
    Produce dwi_preproc.nii.gz/.bvec/.bval y eddy_qc/, igual que dwifslpreproc en la ruta topup.
    """
    G_DEN = "-fslgrad dwi_den.bvec dwi_den.bval"
    print("--Synb0-DisCo (b0 sintetica sin distorsion + topup) y eddy--")
    if all(existe(wd, f) for f in ("dwi_preproc.nii.gz", "dwi_preproc.bvec", "dwi_preproc.bval")):
        print("  [omitido, ya existe] dwi_preproc.nii.gz")
        return
    for d in ("INPUTS", "OUTPUTS"):
        os.makedirs(os.path.join(wd, "synb0", d), exist_ok=True)
    # T1 con craneo si existe (recomendado por Synb0); si no, el T1 sin craneo con --stripped
    t1_in, stripped = (T1_HEAD, False) if datos["t1_head"] else (T1_BRAIN, True)
    # Los vectores de fase de topup y eddy se refieren a los ejes de la imagen en disco: MRtrix los exporta
    # desde pe.json con los strides de FSL, de modo que la b0 de Synb0 y la entrada de eddy coinciden.
    step(["synb0/INPUTS/b0.nii.gz", "synb0/INPUTS/T1.nii.gz", "synb0/b0_pe.txt"], [
        f"dwiextract dwi_den.nii.gz {G_DEN} -bzero - | mrmath - mean -axis 3 b0_pre.nii.gz -force",
        f"mrconvert b0_pre.nii.gz -json_import {PE_FILE} synb0/INPUTS/b0.nii.gz -strides -1,+2,+3 "
        "-datatype float32 -export_pe_topup synb0/b0_pe.txt -force",
        f"cp {t1_in} synb0/INPUTS/T1.nii.gz",
    ], wd)
    acqp = os.path.join(wd, "synb0", "INPUTS", "acqparams.txt")
    if not os.path.exists(acqp):
        linea = open(os.path.join(wd, "synb0", "b0_pe.txt")).read().split()
        with open(acqp, "w") as f:  # b0 real y b0 sintetica (sin distorsion: readout 0)
            f.write(" ".join(linea) + "\n" + " ".join(linea[:3]) + " 0.000\n")
    step("synb0/OUTPUTS/topup_fieldcoef.nii.gz", synb0_cmd(args, licencia, stripped), wd)
    step("b0_synb0.nii.gz", "mrconvert synb0/OUTPUTS/b0_u.nii.gz b0_synb0.nii.gz -force", wd)

    for c in [f"mrconvert dwi_den.nii.gz {G_DEN} -json_import {PE_FILE} eddy_in.nii -strides -1,+2,+3,+4 "
              "-datatype float32 -export_grad_fsl bvecs bvals -export_pe_eddy eddy_acqp.txt eddy_index.txt -force",
              "mrmath synb0/OUTPUTS/b0_all_topup.nii.gz mean -axis 3 b0_topup_mean.nii.gz -force",
              "bet b0_topup_mean.nii.gz b0_topup_brain -m -f 0.2 -n",
              "maskfilter b0_topup_brain_mask.nii.gz dilate eddy_mask.nii.gz -npass 2 -force"]:
        run(c, wd)
    # eddy_acqp.txt tiene solo la linea de la b0 real (FSL >= 6.0.6 rechaza readout 0), igual a la 1a de topup
    eddy_args = ("--imain=eddy_in.nii --mask=eddy_mask.nii.gz --acqp=eddy_acqp.txt --index=eddy_index.txt "
                 "--bvecs=bvecs --bvals=bvals --topup=synb0/OUTPUTS/topup --slm=linear --repol --data_is_shelled "
                 "--out=eddy_out --verbose")
    eddy_bins = [b for b in ("eddy_cuda", "eddy_cpu") if shutil.which(b)] or ["eddy"]
    for i, b in enumerate(eddy_bins):
        try:
            run(f"{b} {eddy_args}", wd)
            break
        except subprocess.CalledProcessError:
            if i == len(eddy_bins) - 1:
                raise
            log_msg(f"AVISO: {b} fallo, se intenta con {eddy_bins[i + 1]}")
    for c in ["mrconvert eddy_out.nii.gz dwi_preproc.nii.gz -fslgrad eddy_out.eddy_rotated_bvecs bvals "
              "-export_grad_fsl dwi_preproc.bvec dwi_preproc.bval -force",
              # mismos archivos que deja dwifslpreproc -eddyqc_text, para el QC
              "rm -rf eddy_qc && mkdir eddy_qc && for f in eddy_out.eddy_*; do case \"$f\" in *.nii*) ;; "
              "*) cp \"$f\" \"eddy_qc/${f#eddy_out.}\";; esac; done",
              "eddy_quad eddy_out -idx eddy_index.txt -par eddy_acqp.txt -m eddy_mask.nii.gz -b bvals -g bvecs "
              "-o eddy_qc/quad || echo 'AVISO: eddy_quad fallo'",
              "rm -f eddy_in.nii eddy_out.nii.gz eddy_out.eddy_outlier_free_data.nii.gz"]:
        run(c, wd)


# ---------------------------------------------------------------------------
# Pipeline por sujeto
# ---------------------------------------------------------------------------

def procesar_sujeto(suj_dir, datos, args, runtime=None, licencia=None):
    global LOG_FILE
    suj = os.path.basename(os.path.normpath(suj_dir))
    wd = os.path.join(suj_dir, WORKDIR_NAME)
    os.makedirs(wd, exist_ok=True)
    if existe(wd, CLEANUP_FILE):
        print(f"  [omitido] {suj} ya se proceso y limpio (-cleanup); para rehacerlo, borrar {WORKDIR_NAME}/")
        return
    LOG_FILE = os.path.join(wd, "pipeline_v2.log")
    sdc = datos["sdc"]
    log_msg(f"INICIO procesamiento (correccion de distorsion: {sdc}, PE {datos['pe']}, readout {datos['trt']} s, "
            f"fuente {datos['fuente_pe']})")
    nt = args.nthreads

    # Registro de entradas y metodo (lo usa el QC). La fase se guarda en los ejes del NIfTI de entrada.
    rel = lambda p: os.path.relpath(p, suj_dir) if p else None
    hecho = json.load(open(os.path.join(wd, METODOS_FILE))).get("sdc") \
        if os.path.exists(os.path.join(wd, METODOS_FILE)) else None
    if hecho and hecho != sdc and existe(wd, "dwi_preproc.nii.gz"):
        raise RuntimeError(f"el sujeto ya se preproceso con -sdc {hecho}; para usar {sdc} hay que borrar "
                           f"{WORKDIR_NAME}/ (o usar -sdc {hecho})")
    with open(os.path.join(wd, PE_FILE), "w") as f:
        json.dump({"PhaseEncodingDirection": datos["pe"], "TotalReadoutTime": datos["trt"]}, f)
    guardar_metodos(wd, pipeline="Procesamiento_DWI.py", sdc=sdc,
                    entradas={"dwi": [rel(s["nii"]) for s in datos["runs"]],
                              "b0_fase_inversa": rel(datos["rev"]["nii"]) if datos["rev"] else None,
                              "t1w": rel(datos["t1_brain"]), "t1w_con_craneo": rel(datos["t1_head"])},
                    PhaseEncodingDirection=datos["pe"], TotalReadoutTime=datos["trt"],
                    fuente_fase=datos["fuente_pe"], versiones=versiones())
    if sdc == "synb0":
        guardar_metodos(wd, synb0={"imagen": synb0_version(args), "runtime": runtime,
                                   "t1": T1_HEAD if datos["t1_head"] else T1_BRAIN + " (--stripped)"})

    # Todas las imagenes se guardan en .nii.gz. Los gradientes van en archivos FSL (bvec/bval) y se pasan
    # con -fslgrad: denoise, degibbs y bias no los modifican, por lo que hay dos tablas:
    #   dwi_den.bvec/bval     (originales, concatenadas)  -> entrada de dwifslpreproc
    #   dwi_preproc.bvec/bval (rotadas por eddy)          -> todos los pasos posteriores
    G_DEN = "-fslgrad dwi_den.bvec dwi_den.bval"
    G = "-fslgrad dwi_preproc.bvec dwi_preproc.bval"

    # ---- Copia de entradas ------------------------------------------------------
    print("--Conversion de series--")
    runs = []
    for i, s in enumerate(datos["runs"]):
        r = f"run{i + 1}"
        step([f"{r}.nii.gz", f"{r}.bvec", f"{r}.bval"],
             f"mrconvert '{s['nii']}' {r}.nii.gz -fslgrad '{s['base']}.bvec' '{s['base']}.bval' "
             f"-export_grad_fsl {r}.bvec {r}.bval -force", wd)
        runs.append(r)
    if sdc == "topup":
        step("b0_PA.nii.gz", f"mrconvert '{datos['rev']['nii']}' b0_PA.nii.gz -force", wd)
    step(T1_BRAIN, f"cp '{datos['t1_brain']}' {T1_BRAIN}", wd)
    if datos["t1_head"]:
        step(T1_HEAD, f"cp '{datos['t1_head']}' {T1_HEAD}", wd)
    # Direccion de fase en los ejes internos de MRtrix (la de pe.json se refiere a los ejes del NIfTI)
    s0 = datos["runs"][0]
    pe_dir, trt = subprocess.run(
        f"mrconvert '{s0['nii']}' -json_import {PE_FILE} -coord 3 0 - -quiet | "
        "mrinfo - -property PhaseEncodingDirection -property TotalReadoutTime",
        shell=True, cwd=wd, capture_output=True, text=True, check=True).stdout.split()

    # ---- Denoise (y degibbs opcional) por run ------------------------------
    print("--DWI Denoise por run--")
    den = []
    for r in runs:
        step(f"{r}_den.nii.gz", f"dwidenoise {r}.nii.gz {r}_den.nii.gz -noise {r}_noise.nii.gz -nthreads {nt} -force", wd)
        if args.degibbs:
            step(f"{r}_den_unr.nii.gz", f"mrdegibbs {r}_den.nii.gz {r}_den_unr.nii.gz -nthreads {nt} -force", wd)
            den.append(f"{r}_den_unr.nii.gz")
        else:
            den.append(f"{r}_den.nii.gz")
    pa_in = "b0_PA.nii.gz"
    if args.degibbs and sdc == "topup":
        step("b0_PA_unr.nii.gz", f"mrdegibbs b0_PA.nii.gz b0_PA_unr.nii.gz -nthreads {nt} -force", wd)
        pa_in = "b0_PA_unr.nii.gz"

    print("--Concatenacion de runs (dwicat)--")
    if len(den) > 1:
        # dwicat no acepta -fslgrad ni escribir a una tuberia ("-"): se usan .mif temporales que se
        # borran en el mismo paso
        tmp = [f"tmp_dwicat_{k}.mif" for k in range(len(den))]
        step(["dwi_den.nii.gz", "dwi_den.bvec", "dwi_den.bval"],
             [f"mrconvert {d} -fslgrad {r}.bvec {r}.bval {t} -force" for d, r, t in zip(den, runs, tmp)] + [
              f"dwicat {' '.join(tmp)} tmp_dwicat_out.mif -nthreads {nt} -force",
              "mrconvert tmp_dwicat_out.mif dwi_den.nii.gz -export_grad_fsl dwi_den.bvec dwi_den.bval -force",
              "rm -f tmp_dwicat_*.mif"], wd)
    else:
        step(["dwi_den.nii.gz", "dwi_den.bvec", "dwi_den.bval"],
             f"cp {den[0]} dwi_den.nii.gz && cp {runs[0]}.bvec dwi_den.bvec && cp {runs[0]}.bval dwi_den.bval", wd)

    # ---- Correccion de distorsion + eddy ------------------------------------------
    # Unica etapa que distingue la adquisicion investigativa (b0 de fase inversa real: topup) de la clinica
    # (Synb0, o solo eddy si no hay correccion posible).
    preproc = ["dwi_preproc.nii.gz", "dwi_preproc.bvec", "dwi_preproc.bval"]
    if sdc == "topup":
        # Par de b0 AP/PA: n_b0 volumenes con la fase de la DWI seguidos de n_b0 con la fase opuesta
        rev = datos["rev"]
        n_pa = nb.load(rev["nii"]).shape[3] if nb.load(rev["nii"]).ndim == 4 else 1
        bvals_run1 = np.loadtxt(datos["runs"][0]["base"] + ".bval")
        n_b0_ap = int(np.sum(bvals_run1 < 50))
        n_b0 = max(1, min(n_pa, n_b0_ap))
        step("b0_pair.nii.gz", [
            f"dwiextract dwi_den.nii.gz {G_DEN} -bzero - | mrconvert - -coord 3 0:{n_b0 - 1} -axes 0,1,2,3 b0_AP.nii.gz -force",
            f"mrconvert {pa_in} -coord 3 0:{n_b0 - 1} -axes 0,1,2,3 b0_PA_sel.nii.gz -force",
            "mrcat b0_AP.nii.gz b0_PA_sel.nii.gz -axis 3 b0_pair.nii.gz -force",
        ], wd)
        print("--DWI dwifslpreproc (topup + eddy)--")
        step(preproc,
             f"dwifslpreproc dwi_den.nii.gz dwi_preproc.nii.gz {G_DEN} -export_grad_fsl dwi_preproc.bvec dwi_preproc.bval "
             f"-rpe_pair -se_epi b0_pair.nii.gz -pe_dir {pe_dir} -readout_time {trt} -align_seepi "
             "-eddy_options ' --slm=linear --repol --data_is_shelled' "
             f"-eddyqc_text eddy_qc -scratch . -nthreads {nt} -force", wd)
    elif sdc == "synb0":
        corregir_synb0(wd, datos, args, runtime, licencia, nt)
    else:
        print("--DWI dwifslpreproc sin correccion de susceptibilidad (solo eddy)--")
        step(preproc,
             f"dwifslpreproc dwi_den.nii.gz dwi_preproc.nii.gz {G_DEN} -export_grad_fsl dwi_preproc.bvec dwi_preproc.bval "
             f"-rpe_none -pe_dir {pe_dir} -readout_time {trt} "
             "-eddy_options ' --slm=linear --repol --data_is_shelled' "
             f"-eddyqc_text eddy_qc -scratch . -nthreads {nt} -force", wd)

    print("--DWI bias correction (ANTs)--")
    step("dwi_preproc_unbiased.nii.gz",
         f"dwibiascorrect ants dwi_preproc.nii.gz dwi_preproc_unbiased.nii.gz {G} -bias bias_field.nii.gz "
         f"-nthreads {nt} -force", wd)

    print(f"--DWI mask ({args.mask})--")
    if args.mask == "hdbet":
        # HD-BET solo acepta 3D: se aplica al b0 medio de los datos preprocesados. Se fuerza la capa
        # GNU de MKL porque numpy (importado aqui) deja MKL_THREADING_LAYER=INTEL, incompatible con torch.
        # Con GPU utilizable (hdbet_gpu) se usa cuda; si aun asi falla (p. ej. sin memoria) se repite con CPU.
        hdbet = "MKL_THREADING_LAYER=GNU hd-bet -i b0_mean_hdbet.nii.gz -o b0_hdbet.nii.gz --save_bet_mask --no_bet_image"
        cpu = f"{hdbet} -device cpu --disable_tta"
        gpu = f"{hdbet} -device cuda || {{ echo 'AVISO: HD-BET fallo con GPU, se repite con CPU'; {cpu}; }}"
        step("dwi_mask.nii.gz", [
            f"dwiextract dwi_preproc_unbiased.nii.gz {G} - -bzero | mrmath - mean b0_mean_hdbet.nii.gz -axis 3 -force",
            gpu if args.hdbet_gpu else cpu,
            "mrconvert b0_hdbet_bet.nii.gz dwi_mask.nii.gz -datatype uint8 -force",
        ], wd)
    else:
        step("dwi_mask.nii.gz", f"dwi2mask dwi_preproc_unbiased.nii.gz {G} dwi_mask.nii.gz -nthreads {nt} -force", wd)

    shells = [int(round(float(b))) for b in subprocess.run(
        f"mrinfo dwi_preproc_unbiased.nii.gz {G} -shell_bvalues", shell=True, cwd=wd,
        capture_output=True, text=True, check=True).stdout.split()]
    shells_b = [b for b in shells if b > 50]
    guardar_metodos(wd, shells=shells, csd_tejidos=3 if len(shells_b) >= 2 else 2)

    # ---- DTI opcional sobre la shell mas baja --------------------------------
    if args.dti:
        print("--DTI (shell b mas baja)--")
        b_low = min(shells_b)
        step(["dti_FA.nii.gz", "dti_MD.nii.gz", "dti_AD.nii.gz", "dti_RD.nii.gz"], [
            f"dwiextract dwi_preproc_unbiased.nii.gz {G} -shells 0,{b_low} - | "
            f"dwi2tensor - -mask dwi_mask.nii.gz tensor.nii.gz -nthreads {nt} -force",
            "tensor2metric tensor.nii.gz -mask dwi_mask.nii.gz -fa dti_FA.nii.gz -adc dti_MD.nii.gz "
            "-ad dti_AD.nii.gz -rd dti_RD.nii.gz -vector dti_V1.nii.gz -force",
        ], wd)

    # ---- Alineacion T1w <-> DWI ----------------------------------------------
    print("--Alineación T1w con DWI--")
    step("bzero.nii.gz",
         f"dwiextract dwi_preproc_unbiased.nii.gz {G} - -bzero | mrmath - mean bzero.nii.gz -axis 3 -force", wd)
    if datos["t1_head"]:
        # epi_reg: registro rigido con boundary-based registration (BBR)
        step("diff2struct_fsl.mat", [
            f"epi_reg --epi=bzero.nii.gz --t1={T1_HEAD} --t1brain={T1_BRAIN} --out=diff2struct_epireg",
            "cp diff2struct_epireg.mat diff2struct_fsl.mat",
        ], wd)
        reg = "epi_reg (BBR)"
    else:
        # Sin T1w con craneo no se puede usar epi_reg, pero BBR solo necesita la sustancia blanca: flirt 6 dof
        # como inicio y BBR con la sustancia blanca segmentada por FAST en el T1w sin craneo (lo que hace epi_reg)
        step("diff2struct_fsl.mat", [
            f"flirt -in bzero.nii.gz -ref {T1_BRAIN} -omat diff2struct_init.mat -dof 6",
            f"fast -t 1 -n 3 -o diff2struct_fast {T1_BRAIN}",
            "fslmaths diff2struct_fast_pve_2 -thr 0.5 -bin diff2struct_fast_wmseg",
            f"flirt -in bzero.nii.gz -ref {T1_BRAIN} -dof 6 -cost bbr -wmseg diff2struct_fast_wmseg "
            "-init diff2struct_init.mat -schedule $FSLDIR/etc/flirtsch/bbr.sch -omat diff2struct_fsl.mat",
        ], wd)
        # sujetos registrados con una version anterior (flirt 6 dof sin BBR) no se rehacen
        reg = "flirt 6 dof + BBR" if existe(wd, "diff2struct_init.mat") else "flirt 6 dof"
    guardar_metodos(wd, registro_dwi_t1=reg)
    step("T1w_acpc_dc_restore_brain_diff_space.nii.gz", [
        f"transformconvert diff2struct_fsl.mat bzero.nii.gz {T1_BRAIN} flirt_import "
        "diff2struct_mrtrix.txt -force",
        f"mrtransform {T1_BRAIN} -linear diff2struct_mrtrix.txt -inverse "
        "T1w_acpc_dc_restore_brain_diff_space.nii.gz -force",
    ], wd)

    print("--Calculo de 5TT--")
    step("5TT.nii.gz",
         f"5ttgen fsl T1w_acpc_dc_restore_brain_diff_space.nii.gz 5TT.nii.gz -premasked -nthreads {nt} -force", wd)

    # ---- CSD multitejido -------------------------------------------------------
    # Con >= 2 shells b > 0 se separan 3 tejidos (WM, GM, CSF). Con una sola shell (clinico) solo hay datos
    # para 2 (WM + CSF): GM y CSF no se distinguen sin una segunda shell.
    print(f"--CSD ({len(shells_b)} shell(s) b>0: {shells_b}); dhollander; msmt_csd; mtnormalise--")
    step(["response_wm.txt", "response_gm.txt", "response_csf.txt"],
         f"dwi2response dhollander dwi_preproc_unbiased.nii.gz {G} response_wm.txt response_gm.txt response_csf.txt "
         f"-mask dwi_mask.nii.gz -voxels RF_voxels.nii.gz -nthreads {nt} -force", wd)
    if len(shells_b) >= 2:
        step(["wmfod.nii.gz", "gm.nii.gz", "csf.nii.gz"],
             f"dwi2fod msmt_csd dwi_preproc_unbiased.nii.gz {G} response_wm.txt wmfod.nii.gz response_gm.txt gm.nii.gz "
             f"response_csf.txt csf.nii.gz -mask dwi_mask.nii.gz -nthreads {nt} -force", wd)
        step("wmfod_norm.nii.gz",
             "mtnormalise wmfod.nii.gz wmfod_norm.nii.gz gm.nii.gz gm_norm.nii.gz csf.nii.gz csf_norm.nii.gz "
             f"-mask dwi_mask.nii.gz -nthreads {nt} -force", wd)
    else:
        step(["wmfod.nii.gz", "csf.nii.gz"],
             f"dwi2fod msmt_csd dwi_preproc_unbiased.nii.gz {G} response_wm.txt wmfod.nii.gz "
             f"response_csf.txt csf.nii.gz -mask dwi_mask.nii.gz -nthreads {nt} -force", wd)
        step("wmfod_norm.nii.gz",
             "mtnormalise wmfod.nii.gz wmfod_norm.nii.gz csf.nii.gz csf_norm.nii.gz "
             f"-mask dwi_mask.nii.gz -nthreads {nt} -force", wd)

    # ---- Tractografia ---------------------------------------------------------
    print("--Inicio de calculo de tractografía Probabilistica--")
    # Semillas en la interfaz sustancia gris - sustancia blanca (ACT, Smith et al. 2012): las fibras parten
    # donde terminan anatomicamente, en vez de en cualquier punto de la mascara DWI (v1).
    # De la tractografia SIFT (TCK_DWI) salen los fasciculos en espacio DWI; la de 21 puntos (TCK_21P) se
    # lleva a MNI para segmentar.
    step([TCK_DWI, TCK_21P], [
        "5tt2gmwmi 5TT.nii.gz gmwmi_seed.nii.gz -force",
        "tckgen -algorithm iFOD2 wmfod_norm.nii.gz output.tck -act 5TT.nii.gz -backtrack -crop_at_gmwmi "
        "-seed_gmwmi gmwmi_seed.nii.gz -maxlength 250 -minlength 25 -select 10M -cutoff 0.06 "
        f"-nthreads {nt} -force",
        f"tcksift output.tck wmfod_norm.nii.gz -act 5TT.nii.gz -term_number 3M {TCK_DWI} -nthreads {nt} -force",
        f"tckresample {TCK_DWI} {TCK_21P} -num_points {N_PUNTOS} -force",
        "rm -f output.tck",
    ], wd)

    # ---- Tractografia a espacio MNI ------------------------------------------
    # ANTs SyN del T1w (en espacio DWI) a MNI152 1 mm. Para tcktransform se necesita, en cada voxel del
    # sujeto, su coordenada MNI: se muestrea una identidad de MNI con la transformacion inversa. La grilla
    # MNI se amplia 40 mm para que ningun punto quede fuera: tcktransform descarta esos puntos y main_index
    # lee las fibras con largo fijo (N_PUNTOS), por lo que una fibra mas corta corrompe la segmentacion.
    print("--Registro de tractografía al espacio MNI152 (ANTs SyN)--")
    if existe(wd, FINAL_FILE):
        print(f"  [omitido, ya existe] {FINAL_FILE}")
    else:
        t1d = "T1w_acpc_dc_restore_brain_diff_space.nii.gz"
        if not os.path.exists(os.path.join(wd, f"{ANTS_PREFIX}1InverseWarp.nii.gz")):
            run(f"antsRegistrationSyN.sh -d 3 -f '{MNI_1mm_brain}' -m {t1d} -t s -n {nt} -o {ANTS_PREFIX}", wd)
        for c in [f"mrgrid '{MNI_1mm_brain}' pad -uniform 40 tmp_mni_pad.nii.gz -force",
                  "warpinit tmp_mni_pad.nii.gz tmp_identity[].nii.gz -force"] + [
                  f"antsApplyTransforms -d 3 -e 0 -i tmp_identity{i}.nii.gz -o tmp_warp{i}.nii.gz -r {t1d} "
                  f"-t [{ANTS_PREFIX}0GenericAffine.mat,1] -t {ANTS_PREFIX}1InverseWarp.nii.gz "
                  "--default-value 2147483647" for i in range(3)] + [
                  "warpcorrect tmp_warp[].nii.gz mrtrix_warp_MNI2dwi.nii.gz -marker 2147483647 -force",
                  "rm -f tmp_mni_pad.nii.gz tmp_identity?.nii.gz tmp_warp?.nii.gz",
                  "warpinvert mrtrix_warp_MNI2dwi.nii.gz mrtrix_warp_dwi2MNI.nii.gz -force",
                  f"tcktransform {TCK_21P} mrtrix_warp_MNI2dwi.nii.gz tmp_MNI.tck -force"]:
            run(c, wd)
        puntos = contar_puntos(os.path.join(wd, "tmp_MNI.tck"))
        if list(puntos) != [N_PUNTOS]:
            raise RuntimeError(f"tcktransform dejo fibras con distinto numero de puntos: {puntos}")
        os.replace(os.path.join(wd, "tmp_MNI.tck"), os.path.join(wd, FINAL_FILE))
        log_msg(f"Tractografia en MNI: {sum(puntos.values())} fibras, todas con {N_PUNTOS} puntos")

    bundles = os.path.join(wd, "tractography_prob_sift_3M_21p_MNI.bundles")
    if not existe(wd, os.path.basename(bundles)):
        print("--Transformacion de formato de tractografia (.tck a .bundles)--")
        log_msg("$ __bundles__ tck2bundles")
        tck2bundles(os.path.join(wd, FINAL_FILE), bundles)

    # ---- Segmentacion de fasciculos -------------------------------------------
    if args.segmentacion is not None:
        print("--Inicio de proceso de segmentación de fascículos--")
        tractography = None  # TCK_DWI (mismo orden de fibras que la de 21 puntos), se carga si hace falta
        for tipo in args.segmentacion:
            atlas_dir, atlas_txt, sufijo = ATLAS[tipo]
            print(f"Calculo de Fascículos {tipo}")
            out = os.path.join(wd, tipo)
            n_atlas = sum(1 for l in open(os.path.join(RECURSOS, atlas_txt)) if l.strip())
            ind_dir = os.path.join(out, "indices_folder")
            res_dwi = os.path.join(out, "results_folder_DWI")
            if os.path.isdir(res_dwi) and len(os.listdir(res_dwi)) == n_atlas:
                print(f"  [omitido, ya existe] {tipo}/results_folder_DWI ({n_atlas} fasciculos)")
                continue
            os.makedirs(res_dwi, exist_ok=True)
            # main_index usa buffers de largo fijo para las rutas: con rutas absolutas largas corrompe
            # memoria y aborta sin escribir resultados. Se ejecuta desde wd con enlaces y rutas cortas.
            for link, destino in ((atlas_dir, os.path.join(RECURSOS, atlas_dir)),
                                  (atlas_txt, os.path.join(RECURSOS, atlas_txt))):
                if os.path.islink(os.path.join(wd, link)):
                    os.remove(os.path.join(wd, link))
                os.symlink(destino, os.path.join(wd, link))
            try:
                run(f"'{MAIN_INDEX}' 21 {os.path.basename(bundles)} subject {atlas_dir}/ {atlas_txt} "
                    f"{tipo}/results_folder {tipo}/indices_folder", wd)
            except subprocess.CalledProcessError:
                # main_index puede abortar al liberar memoria al salir ("free(): invalid next size")
                # despues de escribir todos sus resultados. Se acepta solo si estan todos los indices.
                n_ind = len(os.listdir(ind_dir)) if os.path.isdir(ind_dir) else 0
                if n_ind != n_atlas:
                    raise
                log_msg(f"AVISO: main_index ({tipo}) termino con error al salir, pero genero los "
                          f"{n_ind}/{n_atlas} archivos de indices; se continua.")
            finally:
                for link in (atlas_dir, atlas_txt):
                    if os.path.islink(os.path.join(wd, link)):
                        os.remove(os.path.join(wd, link))
            # Los indices de main_index se refieren al orden de las fibras, que tckresample conserva: los
            # fasciculos en espacio DWI se extraen de la tractografia original, con todos sus puntos.
            if tractography is None:
                tractography = nb.streamlines.load(os.path.join(wd, TCK_DWI)).streamlines
            for fasc_ind in os.listdir(ind_dir):
                nombre_fasc = fasc_ind.removesuffix(".txt").removesuffix(sufijo)
                indices = open(os.path.join(out, "indices_folder", fasc_ind)).readlines()
                indices = [int(ind.rstrip("\n")[:-7]) for ind in indices]
                bundle = Tractogram(streamlines=tractography[indices], affine_to_rasmm=np.eye(4))
                TF.TckFile(bundle, header={'timestamp': 0}).save(
                    os.path.join(out, "results_folder_DWI", nombre_fasc + ".tck"))

        # ---- Metricas por fasciculo -------------------------------------------
        for tipo in args.segmentacion:
            atlas_dir, atlas_txt, sufijo = ATLAS[tipo]
            csv = os.path.join("metricas", f"bundle_metrics_{tipo}.csv")
            if existe(wd, csv):
                print(f"  [omitido, ya existe] {csv}")
            elif not all(os.path.exists(os.path.join(wd, f)) for _, f in METRICAS):
                log_msg(f"AVISO: no hay mapas DTI (usar -dti); se omiten las metricas de {tipo}.")
            else:
                print(f"--Metricas por fasciculo ({tipo})--")
                log_msg(f"$ __metricas__ {csv}")
                metricas_fasciculos(wd, tipo, atlas_txt, sufijo, os.path.join(wd, csv))

    # ---- QC en PDF -------------------------------------------------------------
    qc_pdf = os.path.join(wd, "QC", f"QC_{suj}.pdf")
    if os.path.exists(qc_pdf):
        print(f"  [omitido, ya existe] {os.path.basename(qc_pdf)} (usar qc_dwi.py para regenerarlo)")
    elif not args.no_qc:
        print("--Generacion de QC (PDF)--")
        log_msg("$ __qc__ qc_dwi.generar_qc")
        try:
            import qc_dwi
            qc_dwi.generar_qc(wd, suj)
        except Exception as e:
            # el QC no debe invalidar un procesamiento completo
            log_msg(f"AVISO: no se pudo generar el QC: {e}")

    # ---- Limpieza de intermedios ----------------------------------------------
    if args.cleanup:
        print("--Cleanup de archivos intermedios--")
        if not args.no_qc and not os.path.exists(qc_pdf):
            # sin el PDF no se limpia, para poder revisar o regenerar el QC
            log_msg("AVISO: no existe el PDF de QC; se omite el cleanup.")
        else:
            import qc_dwi
            qc_dwi.copiar_eddy_quad(wd, suj)  # antes de borrar eddy_qc/ (tambien con -no_qc)
            limpiar(wd)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="""Pipeline de preprocesamiento, tractografia y segmentacion de fasciculos para datos
        investigativos (multishell con b0 AP/PA: topup) y clinicos (solo AP o solo PA: Synb0-DisCo), creado por
        Sebastian Navarrete sebastian.navarrete@biomedica.udec.cl basado en MRtrix3, FSL, ANTs y Synb0-DisCo.
        Cada sujeto debe tener las series de difusion convertidas con dcm2niix (.nii.gz/.json/.bval/.bvec), en
        cualquier subcarpeta, o diff.nii.gz/.bval/.bvec, y el T1w procesado con el pipeline estructural HCP
        (T1w_acpc_dc_restore_brain.nii.gz). Los resultados se guardan en <sujeto>/proc_dwi_v2/."""
    )
    parser.add_argument('-folder', type=str, required=True,
                        help="Ruta a la carpeta donde se encuentran los datos de sujetos a procesar.")
    parser.add_argument('-segmentacion', nargs='+', choices=['SWM', 'DWM'], default=None,
                        help="Tipo de segmentación de fascículos: SWM (fibras cortas), DWM (fibras largas), o ambas.")
    parser.add_argument('-degibbs', action='store_true',
                        help="Aplica mrdegibbs. Desactivado por defecto: la adquisicion usa Partial Fourier 6/8.")
    parser.add_argument('-dti', action='store_true',
                        help="Calcula FA/MD/AD/RD con la shell b mas baja (comparable con datos clinicos).")
    parser.add_argument('-mask', choices=['auto', 'hdbet', 'dwi2mask'], default='auto',
                        help="Metodo de la mascara DWI. auto (defecto): hd-bet si esta instalado, si no dwi2mask.")
    parser.add_argument('-cleanup', action='store_true',
                        help="Al terminar cada sujeto borra los archivos intermedios (se conservan resultados, "
                             "bundles y el PDF de QC).")
    parser.add_argument('-no_qc', action='store_true',
                        help="No genera el PDF de QC (QC/QC_<sujeto>.pdf) al final de cada sujeto.")
    parser.add_argument('-sdc', choices=['auto', 'topup', 'synb0', 'none'], default='auto',
                        help="Correccion de distorsion por susceptibilidad. topup: con la b0 de fase inversa "
                             "(investigativo); synb0: b0 sintetica de Synb0-DisCo (clinico, solo AP o solo PA); "
                             "none: solo eddy. auto (defecto): topup si hay b0 inversa, si no synb0 si hay "
                             "Docker y licencia de FreeSurfer, si no none.")
    parser.add_argument('-pe_dir', type=str, default=None, choices=PE_VALIDAS,
                        help="Direccion de fase de la DWI (ejes del NIfTI) si no esta en el .json ni en acq.json.")
    parser.add_argument('-readout', type=float, default=None,
                        help="TotalReadoutTime en s si no esta en el .json ni en acq.json (por defecto 0.05).")
    parser.add_argument('-synb0_image', type=str, default=SYNB0_IMAGE,
                        help=f"Imagen de Synb0-DisCo (por defecto {SYNB0_IMAGE}).")
    parser.add_argument('-fs_license', type=str, default=None,
                        help="Licencia de FreeSurfer para Synb0 (por defecto FS_LICENSE, "
                             "$FREESURFER_HOME/license.txt o ~/license.txt).")
    parser.add_argument('-sujetos', nargs='+', default=None,
                        help="Procesa solo estos sujetos (nombres de carpeta).")
    parser.add_argument('-check', action='store_true',
                        help="Solo valida los datos de cada sujeto, sin procesar.")
    parser.add_argument('-nthreads', type=int, default=os.cpu_count(),
                        help="Número de hilos (por defecto: todos).")
    args = parser.parse_args()
    args.mask = metodo_mascara(args.mask)
    args.hdbet_gpu = args.mask == "hdbet" and hdbet_gpu()
    print(f"Mascara DWI: {args.mask}" + (f" ({'GPU' if args.hdbet_gpu else 'CPU'})" if args.mask == "hdbet" else ""))

    if not os.path.isdir(args.folder):
        print(f"Error: La carpeta '{args.folder}' no existe.")
        return 1

    faltan = recursos_faltantes(args.segmentacion)
    if faltan:
        print(f"Error: faltan en {RECURSOS}: {', '.join(faltan)}.\nProcesamiento_DWI.py y qc_dwi.py deben estar "
              "en la misma carpeta que los atlas, las plantillas MNI y main_index.")
        return 1
    for exe in ("antsRegistrationSyN.sh", "antsApplyTransforms"):
        if not shutil.which(exe):
            print(f"Error: no se encontro {exe} (ANTs) en el PATH.")
            return 1

    folder = os.path.abspath(args.folder)
    sujetos = sorted(d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)))
    if args.sujetos:
        sujetos = [s for s in sujetos if s in args.sujetos]
    print("Sujetos encontrados:", sujetos)

    runtime = licencia = None
    if args.sdc in ("auto", "synb0"):
        runtime, licencia = detectar_docker(), buscar_licencia(args.fs_license)
        print(f"Synb0-DisCo: Docker {'OK' if runtime else 'no disponible'} | licencia FreeSurfer {licencia or 'no encontrada'}")

    # Validacion de todos los sujetos antes de procesar
    validos = {}
    for suj in sujetos:
        datos, errores, avisos = validar_sujeto(os.path.join(folder, suj), args, runtime, licencia)
        if errores:
            print(f"\n[X] {suj}: se omitira")
            for e in errores:
                print(f"     - {e}")
            continue
        validos[suj] = datos
        modo = {"topup": "topup con b0 de fase inversa (investigativo)",
                "synb0": "Synb0-DisCo + topup (clinico)", "none": "sin correccion de susceptibilidad"}
        print(f"\n[OK] {suj}  -> {modo[datos['sdc']]}")
        for s in datos["runs"]:
            bvals = np.atleast_1d(np.loadtxt(s["base"] + ".bval"))
            print(f"     DWI: {s['name']}  ({len(bvals)} vol, b={sorted(set(int(b) for b in np.round(bvals, -2)))})")
        if datos["rev"]:
            print(f"     b0 inversa: {datos['rev']['name']}  (PE {datos['rev']['meta']['PhaseEncodingDirection']})")
        print(f"     PE: {datos['pe']}  readout: {datos['trt']} s  (fuente: {datos['fuente_pe']})")
        print(f"     T1w: {datos['t1_brain']}")
        print(f"     T1w con craneo (epi_reg{', Synb0' if datos['sdc'] == 'synb0' else ''}): "
              f"{datos['t1_head'] or 'no encontrado'}")
        for a in avisos:
            print(f"     ! {a}")

    if args.check:
        return 0

    fallidos = []
    for suj, datos in validos.items():
        print(f"\n===== Procesando sujeto: {suj} =====")
        try:
            procesar_sujeto(os.path.join(folder, suj), datos, args, runtime, licencia)
        except Exception:
            traceback.print_exc()
            print(f"[X] Error procesando {suj}, se continua con el siguiente sujeto.")
            fallidos.append(suj)

    omitidos = [s for s in sujetos if s not in validos]
    print("\nResumen:")
    print("  Procesados:", [s for s in validos if s not in fallidos])
    print("  Con errores:", fallidos)
    print("  Omitidos (validacion):", omitidos)
    return 1 if fallidos else 0


if __name__ == "__main__":
    sys.exit(main())
