# -*- coding: utf-8 -*-
"""
qc_dwi.py

Genera un PDF de control de calidad (QC) para los resultados de Procesamiento_DWI.py, pensado para
detectar rapido los problemas de cada sujeto antes de analisis mas complejos. Todo el QC esta en el PDF.

Paginas:
    1. Resumen: series usadas y correccion de distorsion, movimiento y cortes outlier (eddy), metricas
       con valores de referencia orientativos, fasciculos segmentados y cobertura del FOV.
    2. Correccion de distorsion: b0 sin corregir, b0 de fase inversa (topup) o b0 sintetica (Synb0) y
       b0 corregido, con el borde de sustancia blanca del T1w superpuesto.
    3. Registro T1w -> espacio DWI: FA con borde de sustancia blanca del T1w y T1w con mascara DWI.
    4. Segmentacion 5TT sobre el T1w en espacio de difusion.
    5. Registro T1w -> MNI152: T1w deformado con bordes de la plantilla MNI y viceversa.
    6. Cobertura del FOV: mascara DWI sobre la plantilla MNI y % de fibras de cada fasciculo del atlas
       dentro de ella (separa fasciculos cortados por la adquisicion de fallas de tractografia).
    7. Tractografia: densidad de fibras (TDI) en espacio de difusion y en MNI.
    8-9. Fasciculos DWM / SWM (proyeccion "glass brain" sobre el T1w), con fasciculos no encontrados.

Uso:
    python3 qc_dwi.py <carpeta>          # carpeta de sujetos, de un sujeto o su proc_dwi_v2
    python3 qc_dwi.py <carpeta> -force   # regenera las imagenes intermedias

El PDF se guarda en <sujeto>/proc_dwi_v2/QC/QC_<sujeto>.pdf, junto a una copia del reporte de eddy_quad
(QC_eddy_<sujeto>.pdf). Las imagenes intermedias quedan en <sujeto>/proc_dwi_v2/qc_tmp/.

@author: snava
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
import nibabel as nb
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Plantillas y atlas compartidos (misma ubicacion que en Procesamiento_DWI.py)
RECURSOS = os.path.abspath(os.environ.get(
    "NEUROCISC_RECURSOS", os.path.join(SCRIPT_DIR, "..", "Pipeline_ProcesamientoDiff_TractySegmentacion")))
MNI_1mm_brain = os.path.join(RECURSOS, "MNI152_T1_1mm_brain.nii.gz")
# Atlas en MNI para la cobertura del FOV: (lista de fasciculos, carpeta de .tck, prefijo/sufijo a quitar)
ATLAS_TCK = {"DWM": ("atlas_faisceaux.txt", "atlas_faisceaux_tck", "atlas_", "_MNI"),
             "SWM": ("AtlasRo_estables.txt", "AtlasRo_tck", "", "")}
SDC_TITULO = {"topup": "topup con b0 de fase inversa + eddy",
              "synb0": "Synb0-DisCo (b0 sintética) + topup + eddy",
              "none": "sin corrección de susceptibilidad, solo eddy"}
WORKDIR_NAME = "proc_dwi_v2"
A4 = (11.69, 8.27)

# Valores de referencia orientativos propuestos para este pipeline (no son un estandar publicado ni
# criterios de exclusion): solo sirven para senalar sujetos que conviene revisar visualmente.
REF = {
    "dice": 0.85,           # mascara DWI vs cerebro T1w
    "corr_mni": 0.80,       # T1w deformado vs plantilla MNI
    "pts_in": 0.95,         # fraccion de puntos de streamlines dentro del cerebro
    "mot_rel": 0.5,         # movimiento relativo medio (mm)
    "rot_epireg": 15.0,     # rotacion DWI->T1w (grados), valores grandes sugieren fallo
    "cob_completo": 85.0,   # % de fibras del atlas enteras dentro del FOV: fasciculo completo
    "cob_parcial": 60.0,    # bajo esto, incompleto (entre ambos, parcial). sujetos con FOV completo: >= 92 %
}


# ---------------------------------------------------------------------------
# Preparacion de imagenes
# ---------------------------------------------------------------------------

def sh(cmd, cwd):
    r = subprocess.run(cmd, shell=True, cwd=cwd, executable="/bin/bash",
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Fallo: {cmd}\n{r.stderr[-2000:]}")


def entradas(wd):
    """
    (base de la 1a serie DWI, b0 de fase inversa o sintetica, correccion usada, lista de series) segun
    metodos.json. Sujetos sin metodos.json (versiones anteriores, siempre con topup): se buscan las series
    como lo hace el pipeline.
    """
    wd = os.path.abspath(wd)
    suj_dir = os.path.dirname(wd)
    p = os.path.join(wd, "metodos.json")
    if os.path.exists(p):
        m = json.load(open(p))
        e, sdc = m["entradas"], m["sdc"]
        dwi = [os.path.join(suj_dir, d) for d in e["dwi"]]
        inv = os.path.join(suj_dir, e["b0_fase_inversa"]) if e.get("b0_fase_inversa") else None
    else:
        from Procesamiento_DWI import find_dwi_series, opposite_pe
        runs, b0s, _, _ = find_dwi_series(suj_dir)
        if not runs:
            raise RuntimeError("No se encontraron las series DWI originales del sujeto")
        pe = runs[0]["meta"].get("PhaseEncodingDirection") or ""
        rev = [s["nii"] for s in b0s if opposite_pe(pe, s["meta"].get("PhaseEncodingDirection") or "")]
        dwi, inv, sdc = [s["nii"] for s in runs], (rev[0] if rev else None), ("topup" if rev else "none")
    if sdc == "synb0":
        inv = os.path.join(wd, "b0_synb0.nii.gz")
    series = [os.path.relpath(d, suj_dir) for d in dwi] + ([os.path.relpath(inv, suj_dir)] if inv else [])
    return dwi[0].removesuffix(".nii.gz"), inv, sdc, series


def tck_dwi(wd):
    """Tractografia en espacio DWI: la original (v2) o la de 21 puntos (sujetos de versiones anteriores)."""
    t = "tractography_prob_sift_3M.tck"
    return t if os.path.exists(os.path.join(wd, t)) else "tractography_prob_sift_3M_21p.tck"


def metodo_mni(wd):
    return "ANTs SyN" if os.path.exists(os.path.join(wd, "ants_T12MNI_1Warp.nii.gz")) else "flirt afín + fnirt no lineal"


def preparar(wd, force=False):
    """Lleva todas las imagenes a una grilla comun de 1 mm (T1w en espacio DWI) y a MNI 1 mm."""
    q = os.path.join(wd, "qc_tmp")
    os.makedirs(q, exist_ok=True)
    tpl = "qc_tmp/T1.nii.gz"
    # b0 sin corregir desde las series originales (-cleanup solo conserva estas y dwi_preproc_unbiased)
    ap, inv, _, _ = entradas(wd)
    # T1w y mascara DWI deformados a MNI con ANTs SyN (o con fnirt en sujetos de versiones anteriores; FSL
    # aplica el warp en la grilla del T1w, por eso la mascara se remuestrea antes a esa grilla)
    t1d = "T1w_acpc_dc_restore_brain_diff_space.nii.gz"
    if metodo_mni(wd) == "ANTs SyN":
        ants = f"antsApplyTransforms -d 3 -r '{MNI_1mm_brain}' -t ants_T12MNI_1Warp.nii.gz -t ants_T12MNI_0GenericAffine.mat"
        t1_mni = f"{ants} -i {t1d} -o qc_tmp/T1_MNI.nii.gz"
        mask_mni = f"{ants} -i dwi_mask.nii.gz -n NearestNeighbor -o qc_tmp/mask_MNI.nii.gz"
    else:
        fsl = f"applywarp --ref='{MNI_1mm_brain}' --warp=warps_T12MNI.nii.gz"
        t1_mni = f"{fsl} --in={t1d} --out=qc_tmp/T1_MNI.nii.gz"
        mask_mni = (f"mrgrid dwi_mask.nii.gz regrid -template {t1d} -interp nearest qc_tmp/mask_t1.nii.gz -force "
                    f"-quiet && {fsl} --in=qc_tmp/mask_t1.nii.gz --interp=nn --out=qc_tmp/mask_MNI.nii.gz")
    b0_corr = "dwiextract dwi_preproc_unbiased.nii.gz -fslgrad dwi_preproc.bvec dwi_preproc.bval -bzero - -quiet"
    tareas = [
        ("qc_tmp/T1.nii.gz", "mrgrid T1w_acpc_dc_restore_brain_diff_space.nii.gz regrid -voxel 1 qc_tmp/T1.nii.gz -force"),
        ("qc_tmp/b0_pre.nii.gz", f"dwiextract '{ap}.nii.gz' -fslgrad '{ap}.bvec' '{ap}.bval' -bzero - -quiet | "
                             f"mrmath - mean - -axis 3 -quiet | "
                             f"mrtransform - -template {tpl} -interp linear qc_tmp/b0_pre.nii.gz -force"),
        ("qc_tmp/b0_inv.nii.gz", f"mrconvert '{inv}' -axes 0,1,2,3 - -quiet | mrmath - mean - -axis 3 -quiet | "
                             f"mrtransform - -template {tpl} -interp linear qc_tmp/b0_inv.nii.gz -force"),
        ("qc_tmp/b0_post.nii.gz", f"{b0_corr} | mrmath - mean - -axis 3 -quiet | "
                              f"mrtransform - -template {tpl} -interp linear qc_tmp/b0_post.nii.gz -force"),
        ("qc_tmp/FA.nii.gz", f"mrtransform dti_FA.nii.gz -template {tpl} -interp linear qc_tmp/FA.nii.gz -force"),
        ("qc_tmp/mask.nii.gz", f"mrtransform dwi_mask.nii.gz -template {tpl} -interp nearest qc_tmp/mask.nii.gz -force"),
        ("qc_tmp/5tt.nii.gz", f"mrtransform 5TT.nii.gz -template {tpl} -interp linear qc_tmp/5tt.nii.gz -force"),
        ("qc_tmp/T1_MNI.nii.gz", t1_mni),
        ("qc_tmp/mask_MNI.nii.gz", mask_mni),
        ("qc_tmp/tdi_diff.nii.gz", f"tckmap {tck_dwi(wd)} -template {tpl} qc_tmp/tdi_diff.nii.gz -quiet -force"),
        ("qc_tmp/tdi_MNI.nii.gz", f"tckmap tractography_prob_sift_3M_21p_MNI.tck -template '{MNI_1mm_brain}' "
                              "qc_tmp/tdi_MNI.nii.gz -quiet -force"),
    ]
    # FA solo existe si el pipeline se corrio con -dti; si no, se usa el b0 corregido
    if not os.path.exists(os.path.join(wd, "dti_FA.nii.gz")):
        tareas = [t for t in tareas if t[0] != "qc_tmp/FA.nii.gz"]
    # sin correccion de susceptibilidad no hay b0 de fase inversa ni sintetica
    if inv is None or not os.path.exists(inv):
        tareas = [t for t in tareas if t[0] != "qc_tmp/b0_inv.nii.gz"]
    for salida, cmd in tareas:
        if force or not os.path.exists(os.path.join(wd, salida)):
            print(f"  QC: generando {salida}", flush=True)
            sh(cmd, wd)


def cargar(path):
    img = nb.as_closest_canonical(nb.load(path))
    return np.asanyarray(img.dataobj).astype(np.float32), img.affine


# ---------------------------------------------------------------------------
# Utilidades de visualizacion
# ---------------------------------------------------------------------------

def corte(vol, plano, i):
    if plano == "ax":
        return vol[:, :, i].T
    if plano == "cor":
        return vol[:, i, :].T
    return vol[i, :, :].T


def posiciones(mask, plano, fracs):
    eje = {"sag": 0, "cor": 1, "ax": 2}[plano]
    idx = np.where(mask.any(axis=tuple(a for a in range(3) if a != eje)))[0]
    lo, hi = idx.min(), idx.max()
    return [int(lo + f * (hi - lo)) for f in fracs]


def ventana(vol, mask=None):
    v = vol[mask > 0] if mask is not None else vol[vol > 0]
    if v.size == 0:
        return 0, 1
    return np.percentile(v, 1), np.percentile(v, 99.8)


def panel(ax, bg, plano, i, vmin, vmax, contornos=(), overlay=None, titulo=None):
    ax.imshow(corte(bg, plano, i), cmap="gray", vmin=vmin, vmax=vmax, origin="lower", interpolation="nearest")
    if overlay is not None:
        ov, cmap, ovmin, ovmax, alpha = overlay
        s = np.ma.masked_less_equal(corte(ov, plano, i), ovmin)
        ax.imshow(s, cmap=cmap, vmin=ovmin, vmax=ovmax, alpha=alpha, origin="lower", interpolation="nearest")
    for vol, nivel, color in contornos:
        c = corte(vol, plano, i)
        if c.min() < nivel < c.max():
            ax.contour(c, levels=[nivel], colors=color, linewidths=0.6)
    ax.set_xticks([]); ax.set_yticks([])
    if plano in ("ax", "cor"):
        ax.text(0.97, 0.03, "R", transform=ax.transAxes, color="yellow", fontsize=7, ha="right")
    if titulo:
        ax.set_title(titulo, fontsize=7)


def encuadre(ax, mask, plano):
    """Encuadre cuadrado centrado en el cerebro, igual para los 3 planos (misma escala en todos los paneles)."""
    ind = np.argwhere(mask > 0)
    centro = (ind.min(0) + ind.max(0)) / 2
    semi = (ind.max(0) - ind.min(0)).max() / 2 + 4
    u, v = {"ax": (0, 1), "cor": (0, 2), "sag": (1, 2)}[plano]
    ax.set_facecolor("black")
    ax.set_xlim(centro[u] - semi, centro[u] + semi); ax.set_ylim(centro[v] - semi, centro[v] + semi)


def pagina(pdf, titulo, subtitulo=""):
    fig = plt.figure(figsize=A4)
    fig.suptitle(titulo, fontsize=13, fontweight="bold", x=0.02, ha="left")
    if subtitulo:
        fig.text(0.02, 0.935, subtitulo, fontsize=8, color="0.3", va="top", linespacing=1.4)
    return fig


def fila_cortes(fig, gs, r, bg, vmin, vmax, mask, etiqueta, contornos=(), overlay=None,
                cortes=(("ax", .25), ("ax", .45), ("ax", .65), ("cor", .5), ("sag", .4))):
    for k, (plano, f) in enumerate(cortes):
        ax = fig.add_subplot(gs[r, k])
        i = posiciones(mask, plano, [f])[0]
        panel(ax, bg, plano, i, vmin, vmax, contornos, overlay)
        encuadre(ax, mask, plano)
        if k == 0:
            ax.set_ylabel(etiqueta, fontsize=8)


# ---------------------------------------------------------------------------
# Metricas
# ---------------------------------------------------------------------------

def dice(a, b):
    a, b = a > 0, b > 0
    return 2 * np.sum(a & b) / (np.sum(a) + np.sum(b) + 1e-9)


def fraccion_puntos_dentro(tck, mask_img_path, n=20000):
    """Fraccion de puntos de las primeras n streamlines que caen dentro de la mascara."""
    img = nb.load(mask_img_path)
    data = np.asanyarray(img.dataobj) > 0
    inv = np.linalg.inv(img.affine)
    pts = []
    for k, s in enumerate(nb.streamlines.load(tck, lazy_load=True).streamlines):
        if k >= n:
            break
        pts.append(s)
    if not pts:
        return np.nan
    p = np.vstack(pts)
    v = np.rint(nb.affines.apply_affine(inv, p)).astype(int)
    ok = np.all((v >= 0) & (v < data.shape), axis=1)
    dentro = np.zeros(len(v), bool)
    dentro[ok] = data[v[ok, 0], v[ok, 1], v[ok, 2]]
    return dentro.mean()


def cobertura_fov(mask_mni_path):
    """
    Para cada fasciculo de los atlas (fibras en MNI): % de fibras que quedan enteras dentro de la mascara DWI
    del sujeto llevada a MNI y % medio de sus puntos fuera. No depende de la tractografia del sujeto: indica
    que fasciculos podia capturar la adquisicion. Devuelve {tipo: [(nombre, pct_enteras, pct_fuera), ...]}.
    """
    img = nb.load(mask_mni_path)
    M = np.asanyarray(img.dataobj) > 0
    inv = np.linalg.inv(img.affine)
    res = {}
    for tipo, (txt, carpeta, prefijo, sufijo) in ATLAS_TCK.items():
        filas = []
        if not os.path.isdir(os.path.join(RECURSOS, carpeta)):
            continue
        for l in open(os.path.join(RECURSOS, txt)):
            if not l.strip():
                continue
            nombre = l.split()[0]
            tck = os.path.join(RECURSOS, carpeta, nombre + ".tck")
            if not os.path.exists(tck):
                continue
            sl = nb.streamlines.load(tck).streamlines
            largos = np.asarray(sl._lengths)
            v = np.rint(nb.affines.apply_affine(inv, np.asarray(sl.get_data()))).astype(int)
            ok = np.all((v >= 0) & (v < M.shape), axis=1)
            dentro = np.zeros(len(v), bool)
            dentro[ok] = M[tuple(v[ok].T)]
            frac = np.bincount(np.repeat(np.arange(len(largos)), largos), weights=dentro,
                               minlength=len(largos)) / largos
            filas.append((nombre.removeprefix(prefijo).removesuffix(sufijo), 100 * np.mean(frac == 1),
                          100 * (1 - frac.mean())))
        res[tipo] = filas
    return res


def estado_cobertura(pct):
    return "completo" if pct >= REF["cob_completo"] else ("parcial" if pct >= REF["cob_parcial"] else "incompleto")


def params_fsl_mat(path):
    m = np.loadtxt(path)
    R = m[:3, :3]
    U, _, Vt = np.linalg.svd(R)
    Rn = U @ Vt
    ang = np.degrees(np.arccos(np.clip((np.trace(Rn) - 1) / 2, -1, 1)))
    return np.linalg.norm(m[:3, 3]), ang


def eje_corte(nii):
    """Orientacion anatomica del eje de corte (3er eje de voxel), que es el que usa eddy."""
    aff = nb.load(nii).affine
    k = int(np.argmax(np.abs(aff[:3, 2])))
    return ("sagital", "coronal", "axial")[k]


def lista_fasciculos(atlas_txt, sufijo, carpeta):
    """Compara los fasciculos del atlas con los .tck segmentados (con al menos 1 fibra)."""
    esperados = [l.split()[0].removesuffix(sufijo) for l in open(atlas_txt) if l.strip()]
    encontrados, faltan = [], []
    for nombre in esperados:
        tck = os.path.join(carpeta, nombre + ".tck")
        n = 0
        if os.path.exists(tck):
            n = int(nb.streamlines.TckFile.load(tck, lazy_load=True).header.get("count", 0))
        (encontrados if n > 0 else faltan).append(nombre)
    return esperados, encontrados, faltan


def semaforo(valor, ref, mayor_mejor=True):
    if valor is None or np.isnan(valor):
        return "-"
    ok = valor >= ref if mayor_mejor else valor <= ref
    return "OK" if ok else "REVISAR"


# ---------------------------------------------------------------------------
# Paginas
# ---------------------------------------------------------------------------

def pagina_resumen(pdf, wd, suj, M):
    fig = pagina(pdf, f"QC registro y preprocesamiento DWI — {suj}",
                 f"Generado {time.strftime('%Y-%m-%d %H:%M')} · .../{os.path.relpath(wd, os.path.dirname(os.path.dirname(os.path.dirname(wd))))}")

    # Series de entrada y correccion de distorsion (de metodos.json)
    try:
        _, _, sdc, series = entradas(wd)
        series = series[:-1] if sdc == "synb0" else series  # la b0 sintetica no es una serie de entrada
        m = json.load(open(os.path.join(wd, "metodos.json"))) if os.path.exists(os.path.join(wd, "metodos.json")) else {}
        if m.get("entradas", {}).get("t1w"):
            series.append(m["entradas"]["t1w"])
        metodo = SDC_TITULO[sdc]
        if m.get("PhaseEncodingDirection"):
            metodo += f"  ·  PE {m['PhaseEncodingDirection']}, readout {m.get('TotalReadoutTime')} s ({m.get('fuente_fase')})"
        if m.get("shells"):
            metodo += f"  ·  shells {m['shells']}, CSD {m.get('csd_tejidos')} tejidos"
    except Exception as e:
        series, metodo, sdc = [f"(no disponible: {e})"], "n/d", None
    y = 0.86
    fig.text(0.03, y, "Series de entrada", fontsize=10, fontweight="bold"); y -= 0.03
    for s in series:
        fig.text(0.04, y, s, fontsize=8, family="monospace"); y -= 0.025
    fig.text(0.04, y, "Corrección de distorsión: " + metodo, fontsize=8,
             color="red" if sdc == "none" else "black"); y -= 0.025

    # Tabla de metricas
    y -= 0.02
    fig.text(0.03, y, "Métricas (referencias orientativas*)", fontsize=10, fontweight="bold"); y -= 0.035
    filas = [
        ("Dice máscara DWI vs cerebro T1w", M.get("dice"), f"≥ {REF['dice']}", semaforo(M.get("dice"), REF["dice"])),
        ("Rotación DWI→T1w (epi_reg/flirt) [°]", M.get("rot"), f"≤ {REF['rot_epireg']}",
         semaforo(M.get("rot"), REF["rot_epireg"], False)),
        ("Traslación DWI→T1w [mm]", M.get("tras"), "-", "-"),
        ("Correlación T1w deformado vs MNI", M.get("corr_mni"), f"≥ {REF['corr_mni']}",
         semaforo(M.get("corr_mni"), REF["corr_mni"])),
        ("Puntos de fibras en máscara DWI", M.get("pts_diff"), f"≥ {REF['pts_in']}",
         semaforo(M.get("pts_diff"), REF["pts_in"])),
        ("Puntos de fibras MNI en cerebro MNI", M.get("pts_mni"), f"≥ {REF['pts_in']}",
         semaforo(M.get("pts_mni"), REF["pts_in"])),
        ("Movimiento absoluto medio / máx [mm]", M.get("mot_abs"), "-", "-"),
        ("Movimiento relativo medio [mm]", M.get("mot_rel"), f"≤ {REF['mot_rel']}",
         semaforo(M.get("mot_rel"), REF["mot_rel"], False)),
        ("Cortes outlier eddy (--repol) [%]", M.get("outliers"), "-", "-"),
    ]
    for tipo in ("DWM", "SWM"):
        if f"fasc_{tipo}" in M:
            n_ok, n_tot = M[f"fasc_{tipo}"]
            filas.append((f"Fascículos {tipo} segmentados", f"{n_ok}/{n_tot}", f"{n_tot}/{n_tot}",
                          "OK" if n_ok == n_tot else "REVISAR"))
    for tipo in ("DWM", "SWM"):
        if f"cob_{tipo}" in M:
            n_ok, n_tot = M[f"cob_{tipo}"]
            filas.append((f"Fascículos {tipo} dentro del FOV (pág. 6)", f"{n_ok}/{n_tot}", f"{n_tot}/{n_tot}",
                          "OK" if n_ok == n_tot else "REVISAR"))
    for nombre, val, ref, st in filas:
        if isinstance(val, str):
            txt = val
        elif isinstance(val, tuple):
            txt = " / ".join(f"{v:.2f}" for v in val)
        elif val is None or (isinstance(val, float) and np.isnan(val)):
            txt = "n/d"
        else:
            txt = f"{val:.3f}"
        color = {"OK": "green", "REVISAR": "red"}.get(st, "0.3")
        fig.text(0.04, y, nombre, fontsize=8.5)
        fig.text(0.33, y, txt, fontsize=8.5, family="monospace")
        fig.text(0.42, y, ref, fontsize=8.5, color="0.4")
        fig.text(0.48, y, st, fontsize=8.5, color=color, fontweight="bold")
        y -= 0.03

    y -= 0.01
    fig.text(0.03, y, "* Las referencias no corresponden a un estándar publicado: son valores orientativos propuestos "
                      "para este pipeline,\n   útiles para identificar sujetos que conviene revisar con más detalle en "
                      "las páginas siguientes.\n   No deben usarse como umbral estricto ni como criterio de exclusión; "
                      "se recomienda ajustarlos\n   según la distribución observada en la cohorte.",
             fontsize=7, color="0.3", va="top", linespacing=1.4)
    y -= 0.1
    if os.path.exists(os.path.join(wd, "QC", f"QC_eddy_{suj}.pdf")):
        fig.text(0.03, y, f"Reporte detallado de eddy (eddy_quad): proc_dwi_v2/QC/QC_eddy_{suj}.pdf",
                 fontsize=7, color="0.3")

    # Grafico de movimiento
    rms = os.path.join(wd, "eddy_qc", "eddy_movement_rms")
    if os.path.exists(rms):
        m = np.loadtxt(rms)
        ax = fig.add_axes([0.6, 0.55, 0.37, 0.3])
        ax.plot(m[:, 0], label="absoluto", lw=1)
        ax.plot(m[:, 1], label="relativo", lw=1)
        ax.set_xlabel("volumen", fontsize=8); ax.set_ylabel("mm", fontsize=8)
        ax.set_title("Movimiento estimado por eddy", fontsize=9)
        ax.tick_params(labelsize=7); ax.legend(fontsize=7)
    outl = os.path.join(wd, "eddy_qc", "eddy_outlier_map")
    if os.path.exists(outl):
        o = np.loadtxt(outl, skiprows=1)
        plano = eje_corte(os.path.join(wd, "dwi_preproc_unbiased.nii.gz"))
        ax = fig.add_axes([0.6, 0.17, 0.37, 0.28])
        ax.imshow(o.T, aspect="auto", cmap="Reds", origin="lower", interpolation="nearest",
                  extent=(-0.5, o.shape[0] - 0.5, -0.5, o.shape[1] - 0.5))
        # Grilla cada 1 (fina) y cada 10 (con etiqueta) para ubicar el volumen y corte de cada outlier
        ax.set_xticks(np.arange(0, o.shape[0], 10)); ax.set_yticks(np.arange(0, o.shape[1], 10))
        ax.set_xticks(np.arange(-0.5, o.shape[0], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, o.shape[1], 1), minor=True)
        ax.grid(which="minor", color="0.85", lw=0.2)
        ax.grid(which="major", color="0.55", lw=0.5)
        ax.tick_params(which="minor", length=0)
        ax.set_xlabel("volumen (índice desde 0)", fontsize=8)
        ax.set_ylabel(f"corte {plano} (índice desde 0)", fontsize=8)
        ax.set_title(f"Cortes {plano}es reemplazados por eddy (outliers)", fontsize=9)
        ax.tick_params(labelsize=7)
        vols, cortes = np.nonzero(o)
        try:
            bvals = [int(round(b)) for b in np.loadtxt(os.path.join(wd, "dwi_den.bval"))]
        except Exception:
            bvals = None
        pares = [f"v{v}/c{c}" + (f" (b={bvals[v]})" if bvals and v < len(bvals) else "")
                 for v, c in zip(vols, cortes)]
        if not pares:
            txt = "Sin cortes reemplazados."
        else:
            n_max = 12
            txt = "Volumen/corte reemplazados: " + ", ".join(pares[:n_max])
            if len(pares) > n_max:
                txt += f" ... (+{len(pares) - n_max})"
        import textwrap
        fig.text(0.6, 0.09, textwrap.fill(txt, 75), fontsize=7, va="top", family="monospace")
    pdf.savefig(fig); plt.close(fig)


def pagina_distorsion(pdf, I):
    sdc = I["sdc"]
    sub = ("Contorno cian: borde de sustancia blanca del T1w (5TT). Tras la corrección, el b0 debe "
           "ajustarse al contorno, sobre todo en lóbulos frontales y temporales.")
    if sdc == "none":
        sub += ("\nSin b0 de fase inversa ni Synb0: la distorsión por susceptibilidad NO se corrigió (eddy solo "
                "corrige corrientes inducidas y movimiento).")
    fig = pagina(pdf, f"2. Corrección de distorsión ({SDC_TITULO[sdc]})", sub)
    filas = [("b0_pre", "b0 sin corregir")]
    if "b0_inv" in I:
        filas.append(("b0_inv", "b0 fase inversa\nsin corregir" if sdc == "topup" else "b0 sintética\n(Synb0)"))
    filas.append(("b0_post", "b0 corregido"))
    gs = fig.add_gridspec(3, 5, top=0.9, bottom=0.03, left=0.04, right=0.99, hspace=0.08, wspace=0.03)
    cont = [(I["wm"], 0.5, "cyan")]
    for r, (key, etiqueta) in enumerate(filas):
        vmin, vmax = ventana(I[key], I["brain"])
        fila_cortes(fig, gs, r, I[key], vmin, vmax, I["brain"], etiqueta, cont)
    pdf.savefig(fig); plt.close(fig)


def pagina_dwi_t1(pdf, I):
    fig = pagina(pdf, "3. Registro T1w → espacio DWI",
                 "La transformación rígida se estima de DWI a T1w (epi_reg sobre el b0 corregido) y se aplica invertida para "
                 "llevar el T1w (y luego el 5TT) al espacio de difusión,\ndonde se hace la tractografía. "
                 "Fila 1: FA (calculada dentro de la máscara DWI) con el borde de sustancia blanca del T1w (cian). "
                 "Fila 2: T1w en espacio DWI con el contorno\nde la máscara DWI (amarillo) y FA = 0.3 (rojo). "
                 "Si el registro es correcto, los haces con FA alta coinciden con la sustancia blanca del T1w.")
    gs = fig.add_gridspec(2, 5, top=0.86, bottom=0.03, left=0.04, right=0.99, hspace=0.08, wspace=0.03)
    bg = I.get("FA", I["b0_post"])
    vmin, vmax = (0, 0.8) if "FA" in I else ventana(bg, I["brain"])
    fila_cortes(fig, gs, 0, bg, vmin, vmax, I["brain"], "FA" if "FA" in I else "b0", [(I["wm"], 0.5, "cyan")])
    t1min, t1max = ventana(I["T1"], I["brain"])
    cont2 = [(I["mask"], 0.5, "yellow")] + ([(I["FA"], 0.3, "red")] if "FA" in I else [])
    fila_cortes(fig, gs, 1, I["T1"], t1min, t1max, I["brain"], "T1w en espacio DWI", cont2)
    pdf.savefig(fig); plt.close(fig)


TEJIDOS_5TT = [("Sustancia gris cortical", (1.0, 0.15, 0.15)), ("Sustancia gris subcortical", (1.0, 0.0, 1.0)),
               ("Sustancia blanca", (0.0, 0.55, 1.0)), ("LCR", (1.0, 0.9, 0.0))]


def pagina_5tt(pdf, I):
    fig = pagina(pdf, "4. Segmentación 5TT (5ttgen fsl) en espacio DWI",
                 "Cada vóxel se colorea con el tejido de mayor fracción (≥ 0.5). Los bordes deben seguir la anatomía del "
                 "T1w: la sustancia gris cortical sobre la corteza, la subcortical en ganglios basales y tálamo.")
    gs = fig.add_gridspec(2, 5, top=0.88, bottom=0.08, left=0.04, right=0.99, hspace=0.08, wspace=0.03)
    tt = I["5tt"][..., :4]
    etiqueta = np.argmax(tt, axis=-1)
    valido = tt.max(axis=-1) >= 0.5
    rgb = np.zeros(tt.shape[:3] + (3,), np.float32)
    for v, (_, c) in enumerate(TEJIDOS_5TT):
        rgb[(etiqueta == v) & valido] = c
    t1min, t1max = ventana(I["T1"], I["brain"])
    cortes = [("ax", .2), ("ax", .35), ("ax", .5), ("ax", .65), ("ax", .8),
              ("cor", .35), ("cor", .5), ("cor", .65), ("sag", .3), ("sag", .5)]
    for k, (plano, f) in enumerate(cortes):
        ax = fig.add_subplot(gs[k // 5, k % 5])
        i = posiciones(I["brain"], plano, [f])[0]
        ax.imshow(corte(I["T1"], plano, i), cmap="gray", vmin=t1min, vmax=t1max, origin="lower")
        c = np.stack([corte(rgb[..., j], plano, i) for j in range(3)], -1)
        alpha = (c.sum(-1) > 0) * 0.45
        ax.imshow(np.dstack([c, alpha]), origin="lower", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        encuadre(ax, I["brain"], plano)
        if plano in ("ax", "cor"):
            ax.text(0.97, 0.03, "R", transform=ax.transAxes, color="white", fontsize=7, ha="right")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=n) for n, c in TEJIDOS_5TT], loc="lower center", ncol=4,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.01))
    pdf.savefig(fig); plt.close(fig)


def pagina_mni(pdf, I):
    fig = pagina(pdf, f"5. Registro T1w → plantilla MNI152 ({I['metodo_mni']})",
                 "Los atlas de fascículos (DWM y SWM) están definidos en el espacio MNI152. Para segmentar, la tractografía "
                 "se lleva a ese espacio con la deformación calculada\nentre el T1w del sujeto y la plantilla MNI. "
                 "Esta página verifica esa deformación.")
    gs = fig.add_gridspec(2, 5, top=0.86, bottom=0.2, left=0.04, right=0.99, hspace=0.08, wspace=0.03)
    mni, t1m, mm = I["MNI"], I["T1_MNI"], I["MNI"] > 0
    nivel_mni = np.percentile(mni[mm], 60)
    nivel_t1 = np.percentile(t1m[t1m > 0], 60) if np.any(t1m > 0) else 0
    vmin, vmax = ventana(t1m, mm)
    fila_cortes(fig, gs, 0, t1m, vmin, vmax, mm, "T1w del sujeto\ndeformado a MNI",
                [(mni, nivel_mni, "red"), (mm.astype(float), .5, "yellow")])
    vmin, vmax = ventana(mni, mm)
    fila_cortes(fig, gs, 1, mni, vmin, vmax, mm, "plantilla MNI152", [(t1m, nivel_t1, "lime")])
    fig.text(0.04, 0.17, "Cómo leerla", fontsize=9, fontweight="bold")
    fig.text(0.04, 0.155,
             "• Fila 1: imagen = T1w del sujeto después de deformarlo a MNI. Líneas = bordes de la plantilla MNI: "
             "rojo, límite entre sustancia gris y blanca; amarillo, borde del cerebro.\n"
             "• Fila 2: lo mismo al revés. Imagen = plantilla MNI. Línea verde = límite entre sustancia gris y blanca "
             "del T1w deformado del sujeto.\n"
             "• Registro correcto: las líneas siguen los mismos surcos, ventrículos y contorno del cerebro de la imagen "
             "de fondo, en todos los cortes.\n"
             "• Registro incorrecto: desplazamientos o rotaciones sistemáticos, ventrículos que no coinciden, o "
             "cerebro del sujeto que sobresale del borde amarillo.",
             fontsize=8, va="top", linespacing=1.6)
    pdf.savefig(fig); plt.close(fig)


def pagina_cobertura(pdf, I, cob):
    fig = pagina(pdf, "6. Cobertura del FOV de la adquisición",
                 "Fila superior: plantilla MNI con la máscara DWI del sujeto llevada a MNI; en rojo, el cerebro que "
                 "quedó fuera de la adquisición (p. ej. vértex cortado). Abajo: % de fibras de cada fascículo\n"
                 "del atlas que quedan enteras dentro de la máscara. No depende de la tractografía: un fascículo "
                 f"incompleto aquí sale truncado por la adquisición (sus métricas no son comparables), no por una "
                 f"falla\ndel procesamiento. Completo ≥ {REF['cob_completo']:.0f} %, parcial ≥ "
                 f"{REF['cob_parcial']:.0f} %, incompleto < {REF['cob_parcial']:.0f} % (orientativo).")
    gs = fig.add_gridspec(1, 5, top=0.86, bottom=0.6, left=0.04, right=0.99, wspace=0.03)
    mni, mm = I["MNI"], I["MNI"] > 0
    fuera = (mm & (I["mask_MNI"] <= 0)).astype(np.float32)
    vmin, vmax = ventana(mni, mm)
    fila_cortes(fig, gs, 0, mni, vmin, vmax, mm, "máscara DWI\nen MNI", [(I["mask_MNI"], 0.5, "yellow")],
                overlay=(fuera, "autumn_r", 0.5, 1, 0.55),
                cortes=(("ax", .35), ("ax", .65), ("ax", .88), ("cor", .5), ("sag", .5)))
    color = {"completo": "tab:green", "parcial": "tab:orange", "incompleto": "tab:red"}
    if "DWM" in cob:
        filas = cob["DWM"]
        ax = fig.add_axes([0.1, 0.04, 0.33, 0.5])
        ax.barh(range(len(filas)), [f[1] for f in filas], color=[color[estado_cobertura(f[1])] for f in filas])
        ax.set_yticks(range(len(filas)), [f[0] for f in filas], fontsize=5.5)
        ax.invert_yaxis()
        for ref in (REF["cob_completo"], REF["cob_parcial"]):
            ax.axvline(ref, color="0.4", lw=0.6, ls="--")
        ax.set_xlim(0, 100); ax.tick_params(axis="x", labelsize=7)
        ax.set_xlabel("% de fibras del atlas dentro del FOV", fontsize=7)
        ax.set_title("Fascículos DWM", fontsize=9)
    if "SWM" in cob:
        filas = sorted(cob["SWM"], key=lambda f: f[1])
        ax = fig.add_axes([0.53, 0.33, 0.44, 0.21])
        ax.bar(range(len(filas)), [f[1] for f in filas], width=1.0,
               color=[color[estado_cobertura(f[1])] for f in filas])
        for ref in (REF["cob_completo"], REF["cob_parcial"]):
            ax.axhline(ref, color="0.4", lw=0.6, ls="--")
        ax.set_ylim(0, 100); ax.set_xlim(-0.5, len(filas) - 0.5); ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("% fibras dentro del FOV", fontsize=7)
        n = {e: sum(estado_cobertura(f[1]) == e for f in filas) for e in color}
        ax.set_title(f"Fascículos SWM ordenados ({n['completo']} completos, {n['parcial']} parciales, "
                     f"{n['incompleto']} incompletos)", fontsize=9)
        import textwrap
        y = 0.27
        for e in ("incompleto", "parcial"):
            nombres = [f[0] for f in filas if estado_cobertura(f[1]) == e]
            if nombres:
                txt = textwrap.fill(f"{e.capitalize()}s ({len(nombres)}): " + ", ".join(nombres), 110)
                lineas = txt.split("\n")
                if len(lineas) > 6:
                    txt = "\n".join(lineas[:6]) + " ..."
                fig.text(0.53, y, txt, fontsize=6, color=color[e], va="top", family="monospace")
                y -= 0.013 * min(len(lineas), 6) + 0.02
    pdf.savefig(fig); plt.close(fig)


def pagina_tractografia(pdf, I):
    fig = pagina(pdf, "7. Tractografía (densidad de fibras, TDI)",
                 "Fila 1: TDI en espacio DWI sobre T1w. Fila 2: TDI de la tractografía transformada a MNI sobre la "
                 "plantilla: debe quedar dentro del cerebro y siguiendo la sustancia blanca.")
    gs = fig.add_gridspec(2, 5, top=0.9, bottom=0.03, left=0.04, right=0.99, hspace=0.08, wspace=0.03)
    for r, (bgk, tdik, maskk, etiqueta) in enumerate([("T1", "tdi_diff", "brain", "espacio DWI"),
                                                     ("MNI", "tdi_MNI", "MNI_mask", "espacio MNI")]):
        bg, tdi, m = I[bgk], I[tdik], I[maskk]
        vmin, vmax = ventana(bg, m)
        tmax = np.percentile(tdi[tdi > 0], 99) if np.any(tdi > 0) else 1
        fila_cortes(fig, gs, r, bg, vmin, vmax, m, etiqueta, overlay=(tdi, "hot", tmax * 0.02, tmax, 0.7))
    pdf.savefig(fig); plt.close(fig)


def pagina_fasciculos(pdf, I, carpeta, titulo, max_fibras, leyenda, atlas):
    esperados, encontrados, faltan = atlas
    if not esperados:
        return
    tcks = [os.path.join(carpeta, n + ".tck") for n in encontrados]
    sub = (f"{len(encontrados)}/{len(esperados)} fascículos del atlas segmentados (máx. {max_fibras} fibras "
           "dibujadas por fascículo). Proyección sobre el T1w en espacio DWI: deben ubicarse en la anatomía esperada.")
    fig = pagina(pdf, titulo, sub)
    if faltan:
        import textwrap
        fig.text(0.02, 0.9, textwrap.fill(f"No segmentados ({len(faltan)}): " + ", ".join(faltan), 160),
                 fontsize=7, color="red", va="top")
    inv = np.linalg.inv(I["aff"])
    cmap = plt.get_cmap("tab20" if len(tcks) <= 20 else "hsv")
    ejes = [("Axial", (0, 1), 2), ("Coronal", (0, 2), 1), ("Sagital", (1, 2), 0)]
    # Paneles cuadrados del mismo tamaño y escala (1 vóxel = 1 mm) centrados en el cerebro
    ind = np.argwhere(I["brain"] > 0)
    centro = (ind.min(0) + ind.max(0)) / 2
    semi = (ind.max(0) - ind.min(0)).max() / 2 + 4
    ancho = 0.25 if leyenda else 0.31
    alto = ancho * A4[0] / A4[1]
    axs = [fig.add_axes([0.02 + k * (ancho + 0.01), 0.47 - alto / 2, ancho, alto]) for k in range(3)]
    t1 = I["T1"] * (I["brain"] > 0)
    for ax, (nombre, (u, v), eje) in zip(axs, ejes):
        ax.set_facecolor("black")
        ax.imshow(t1.mean(axis=eje).T, cmap="gray", origin="lower")
        ax.set_xlim(centro[u] - semi, centro[u] + semi); ax.set_ylim(centro[v] - semi, centro[v] + semi)
        ax.set_aspect("equal")
        ax.set_title(nombre, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        if eje != 0:
            ax.text(0.97, 0.03, "R", transform=ax.transAxes, color="yellow", fontsize=8, ha="right")
    rng = np.random.default_rng(0)
    handles = []
    for n, t in enumerate(tcks):
        sl = list(nb.streamlines.load(t).streamlines)
        if len(sl) > max_fibras:
            sl = [sl[i] for i in rng.choice(len(sl), max_fibras, replace=False)]
        vox = [nb.affines.apply_affine(inv, s) for s in sl]
        color = cmap(n / max(1, len(tcks) - 1))
        for ax, (_, (u, v), _) in zip(axs, ejes):
            ax.add_collection(LineCollection([s[:, [u, v]] for s in vox], colors=[color], linewidths=0.3, alpha=0.6))
        handles.append(plt.Line2D([], [], color=color, lw=2,
                                  label=os.path.basename(t)[:-4].removeprefix("atlas_")))
    if leyenda:
        fig.legend(handles=handles, loc="center right", fontsize=5.5, ncol=1, frameon=False,
                   bbox_to_anchor=(1.0, 0.47))
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# Principal
# ---------------------------------------------------------------------------

def copiar_eddy_quad(wd, suj):
    """Copia el reporte de eddy_quad a QC/ (eddy_qc/ se borra con -cleanup)."""
    quad = os.path.join(wd, "eddy_qc", "quad", "qc.pdf")
    destino = os.path.join(wd, "QC", f"QC_eddy_{suj}.pdf")
    if os.path.exists(quad) and not os.path.exists(destino):
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        shutil.copy(quad, destino)


def generar_qc(wd, suj=None, force=False):
    wd = os.path.abspath(wd)
    suj = suj or os.path.basename(os.path.dirname(wd))
    preparar(wd, force)
    q = lambda f: os.path.join(wd, "qc_tmp", f)

    I = {}
    I["T1"], I["aff"] = cargar(q("T1.nii.gz"))
    I["brain"] = (I["T1"] > 0).astype(np.float32)
    for k in ("b0_pre", "b0_post", "mask", "tdi_diff", "T1_MNI", "tdi_MNI", "mask_MNI"):
        I[k], _ = cargar(q(f"{k}.nii.gz"))
    I["sdc"] = entradas(wd)[2]
    if I["sdc"] != "none" and os.path.exists(q("b0_inv.nii.gz")):
        I["b0_inv"], _ = cargar(q("b0_inv.nii.gz"))
    if os.path.exists(q("FA.nii.gz")):
        I["FA"], _ = cargar(q("FA.nii.gz"))
    I["5tt"], _ = cargar(q("5tt.nii.gz"))
    I["wm"] = I["5tt"][..., 2]
    I["MNI"], _ = cargar(MNI_1mm_brain)
    I["metodo_mni"] = metodo_mni(wd)
    I["MNI_mask"] = (I["MNI"] > 0).astype(np.float32)

    # Metricas
    M = {"dice": dice(I["mask"], I["brain"])}
    mat = os.path.join(wd, "diff2struct_fsl.mat")
    if os.path.exists(mat):
        M["tras"], M["rot"] = params_fsl_mat(mat)
    mm = I["MNI_mask"] > 0
    M["corr_mni"] = float(np.corrcoef(I["T1_MNI"][mm], I["MNI"][mm])[0, 1])
    M["pts_diff"] = fraccion_puntos_dentro(os.path.join(wd, tck_dwi(wd)), q("mask.nii.gz"))
    M["pts_mni"] = fraccion_puntos_dentro(os.path.join(wd, "tractography_prob_sift_3M_21p_MNI.tck"), MNI_1mm_brain)
    rms = os.path.join(wd, "eddy_qc", "eddy_movement_rms")
    if os.path.exists(rms):
        m = np.loadtxt(rms)
        M["mot_abs"] = (m[:, 0].mean(), m[:, 0].max())
        M["mot_rel"] = m[1:, 1].mean()
    outl = os.path.join(wd, "eddy_qc", "eddy_outlier_map")
    if os.path.exists(outl):
        M["outliers"] = 100 * np.loadtxt(outl, skiprows=1).mean()
    atlas = {}
    for tipo, txt, sufijo in (("DWM", "atlas_faisceaux.txt", "_MNI"), ("SWM", "AtlasRo_estables.txt", "")):
        carpeta = os.path.join(wd, tipo, "results_folder_DWI")
        if os.path.isdir(carpeta) and os.path.exists(os.path.join(RECURSOS, txt)):
            atlas[tipo] = lista_fasciculos(os.path.join(RECURSOS, txt), sufijo, carpeta)
            M[f"fasc_{tipo}"] = (len(atlas[tipo][1]), len(atlas[tipo][0]))
    cob = cobertura_fov(q("mask_MNI.nii.gz"))
    for tipo, filas in cob.items():
        M[f"cob_{tipo}"] = (sum(estado_cobertura(f[1]) == "completo" for f in filas), len(filas))

    os.makedirs(os.path.join(wd, "QC"), exist_ok=True)
    copiar_eddy_quad(wd, suj)
    out = os.path.join(wd, "QC", f"QC_{suj}.pdf")
    with PdfPages(out) as pdf:
        pagina_resumen(pdf, wd, suj, M)
        pagina_distorsion(pdf, I)
        pagina_dwi_t1(pdf, I)
        pagina_5tt(pdf, I)
        pagina_mni(pdf, I)
        pagina_cobertura(pdf, I, cob)
        pagina_tractografia(pdf, I)
        if "DWM" in atlas:
            pagina_fasciculos(pdf, I, os.path.join(wd, "DWM", "results_folder_DWI"),
                              "8. Fascículos de fibras largas (DWM)", 150, True, atlas["DWM"])
        if "SWM" in atlas:
            pagina_fasciculos(pdf, I, os.path.join(wd, "SWM", "results_folder_DWI"),
                              "9. Fascículos de fibras cortas (SWM)", 40, False, atlas["SWM"])
        info = pdf.infodict()
        info["Title"] = f"QC pipeline DWI - {suj}"
    print(f"  QC guardado en {out}")
    return out


def buscar_sujetos(ruta):
    ruta = os.path.abspath(ruta)
    if os.path.basename(ruta) == WORKDIR_NAME:
        return [(os.path.basename(os.path.dirname(ruta)), ruta)]
    if os.path.isdir(os.path.join(ruta, WORKDIR_NAME)):
        return [(os.path.basename(ruta), os.path.join(ruta, WORKDIR_NAME))]
    return [(d, os.path.join(ruta, d, WORKDIR_NAME)) for d in sorted(os.listdir(ruta))
            if os.path.isdir(os.path.join(ruta, d, WORKDIR_NAME))]


def main():
    parser = argparse.ArgumentParser(description="Genera el PDF de QC de Procesamiento_DWI.py")
    parser.add_argument("carpeta", help="Carpeta de sujetos, de un sujeto, o su proc_dwi_v2.")
    parser.add_argument("-force", action="store_true", help="Regenera las imagenes intermedias de qc_tmp/.")
    args = parser.parse_args()
    for suj, wd in buscar_sujetos(args.carpeta):
        if not os.path.exists(os.path.join(wd, "tractography_prob_sift_3M_21p_MNI.tck")):
            print(f"[X] {suj}: pipeline incompleto, se omite el QC")
            continue
        print(f"== QC {suj}")
        try:
            generar_qc(wd, suj, args.force)
        except Exception as e:
            print(f"[X] {suj}: error generando QC: {e}")


if __name__ == "__main__":
    main()
