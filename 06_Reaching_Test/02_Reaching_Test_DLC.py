#!/usr/bin/env python3
"""
Batch processor for DeepLabCut CSVs with:
- Palm Advance/Grasp/Retract phases, grasp accuracy/distance metrics
- mouth bouts and metrics
- Two figures per mouse: (1) traces+ellipses, (2) grasps-only (dots + ellipse)
- Summary CSV + per-bout palm phase CSV

Install:
    pip install pandas numpy matplotlib
"""

from pathlib import Path
from typing import List, Tuple, Optional, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse

# ========= USER SETTINGS =========
DATA_DIR = Path(r"C:\Users\ioann\OneDrive\Desktop\DLC_reach_mito_day1")
CSV_PATTERN = "*.csv"
OUT_DIR = DATA_DIR / "reaching_mito_d1_results"
SHOW_PLOTS = False  # True to show plots, False to just save PDFs

NOSE_NAME   = "snout"
PALM_NAME   = "palm"
SPOUT_NAME  = "spout"
MOUTH_NAME = "mouth"   # <— NEW (add aliases below if needed)

LIKELIHOOD_THR = 0.8
MAX_DX_FROM_SPOUT = 200        # keep points with |x - spout_x| ≤ this (pixels)
MIN_SEG_LEN = 6                # contiguous frames to form a bout
FRAMERATE_FPS = 30             # video frame rate

# Palm/mouth bout exclusion (by horizontal span)
MIN_PALM_BOUT_XRANGE_PX   = 10
MIN_MOUTH_BOUT_XRANGE_PX = 10  # <— NEW

# Grasp detection thresholds (palm)
GRASP_RADIUS_X_PX   = 15
GRASP_SPEED_THR_PX_S = 5.0
GRASP_MIN_FRAMES     = 2

# Grasp accuracy radius (Euclidean)
ACCURACY_RADIUS_PX = 0.045

# Plot style
PALM_COLOR   = "blue"
SNOUT_COLOR  = "red"
MOUTH_COLOR = "darkorange"  # <— NEW
PALM_ALPHA  = 0.6
SNOUT_ALPHA = 0.35
MOUTH_ALPHA = 0.45          # <— NEW
LINE_WIDTH  = 1.5

# Mean-trajectory style (palm)
MEAN_TRAJ_POINTS = 200
MEAN_TRAJ_COLOR  = "black"
MEAN_TRAJ_ALPHA  = 0.9
MEAN_TRAJ_WIDTH  = 3.0

# Covariance ellipse style
ELLIPSE_CONFIDENCE = 0.95      # 0.90, 0.95, 0.99 supported
ELLIPSE_LW = 2.0
ELLIPSE_ALPHA = 0.9

# Grasps-only plot style
GRASP_COLOR = "seagreen"
GRASP_ALPHA = 0.7
GRASP_DOT_SIZE = 10
# =================================


def load_dlc_multiindex(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=[0, 1, 2])
    if not isinstance(df.columns, pd.MultiIndex):
        raise RuntimeError("Expected DLC CSV with 3-row MultiIndex header.")
    return df


def get_main_scorer(df: pd.DataFrame) -> str:
    lvl0 = list(df.columns.get_level_values(0).unique())
    candidates = [v for v in lvl0 if str(v).lower() != "scorer"]
    return candidates[0] if candidates else lvl0[0]


def pick_bodypart(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    bodyparts = list(df.columns.get_level_values(1).unique())
    m = {str(bp).lower(): bp for bp in bodyparts}
    for c in candidates:
        if c.lower() in m:
            return m[c.lower()]
    return None


def extract_xyz(df: pd.DataFrame, scorer: str, bodypart: str):
    sub = df.xs((scorer, bodypart), axis=1, level=[0, 1])
    cols = {c.lower(): c for c in sub.columns}
    x = sub[cols["x"]].astype(float).to_numpy()
    y = sub[cols["y"]].astype(float).to_numpy()
    l = sub[cols["likelihood"]].astype(float).to_numpy()
    return x, y, l


def find_contiguous_segments(mask: np.ndarray, min_len: int = 3) -> List[Tuple[int, int]]:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    segs, start, prev = [], idx[0], idx[0]
    for i in idx[1:]:
        if i != prev + 1:
            segs.append((start, prev))
            start = i
        prev = i
    segs.append((start, prev))
    return [s for s in segs if (s[1] - s[0] + 1) >= min_len]


def total_path_length(x: np.ndarray, y: np.ndarray, mask: np.ndarray, min_len: int) -> float:
    total = 0.0
    for s, e in find_contiguous_segments(mask, min_len=min_len):
        xs, ys = x[s:e+1], y[s:e+1]
        total += float(np.nansum(np.hypot(np.diff(xs), np.diff(ys))))
    return total


def bout_metrics(x: np.ndarray, y: np.ndarray, mask: np.ndarray, min_len: int, fps: float):
    """Basic bout metrics: count, avg path length, avg speed."""
    segs = find_contiguous_segments(mask, min_len=min_len)
    lengths, speeds = [], []
    for s, e in segs:
        xs, ys = x[s:e+1], y[s:e+1]
        L = float(np.nansum(np.hypot(np.diff(xs), np.diff(ys))))
        dur_s = max((e - s) / float(fps), 0.0)
        if dur_s > 0:
            lengths.append(L)
            speeds.append(L / dur_s)
        else:
            lengths.append(L)
    n = len(segs)
    avg_len = float(np.nan) if not lengths else float(np.nanmean(lengths))
    avg_spd = float(np.nan) if not speeds else float(np.nanmean(speeds))
    return n, avg_len, avg_spd


def bout_tortuosity(x: np.ndarray, y: np.ndarray, mask: np.ndarray, min_len: int) -> List[float]:
    vals = []
    for s, e in find_contiguous_segments(mask, min_len=min_len):
        xs, ys = x[s:e+1], y[s:e+1]
        if xs.size < 2:
            continue
        path = float(np.nansum(np.hypot(np.diff(xs), np.diff(ys))))
        disp = float(np.hypot(xs[-1] - xs[0], ys[-1] - ys[0]))
        if np.isfinite(path) and disp > 0:
            vals.append(path / disp)
    return vals


def trajectory_variance_px2(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> float:
    xs = x[mask]; ys = y[mask]
    valid = np.isfinite(xs) & np.isfinite(ys)
    if not np.any(valid):
        return float("nan")
    vx = float(np.nanvar(xs[valid])); vy = float(np.nanvar(ys[valid]))
    return vx + vy


def mean_trajectory_over_bouts(x: np.ndarray, y: np.ndarray, mask: np.ndarray, min_len: int, n_points: int = 100):
    segs = find_contiguous_segments(mask, min_len=min_len)
    curves_x, curves_y = [], []
    for s, e in segs:
        xs, ys = x[s:e+1], y[s:e+1]
        valid = np.isfinite(xs) & np.isfinite(ys)
        xs, ys = xs[valid], ys[valid]
        if xs.size < 2:
            continue
        d = np.hypot(np.diff(xs), np.diff(ys))
        cum = np.concatenate(([0.0], np.cumsum(d)))
        if not np.isfinite(cum[-1]) or cum[-1] == 0:
            continue
        t = cum / cum[-1]
        grid = np.linspace(0.0, 1.0, n_points)
        xi = np.interp(grid, t, xs)
        yi = np.interp(grid, t, ys)
        curves_x.append(xi); curves_y.append(yi)
    if len(curves_x) == 0:
        return None, None
    mean_x = np.nanmean(np.vstack(curves_x), axis=0)
    mean_y = np.nanmean(np.vstack(curves_y), axis=0)
    return mean_x, mean_y


def add_covariance_ellipse(ax, x: np.ndarray, y: np.ndarray, mask: np.ndarray,
                           confidence: float, color: str, lw: float, alpha: float,
                           zorder: int):
    xs = x[mask]; ys = y[mask]
    valid = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[valid]; ys = ys[valid]
    if xs.size < 2:
        return None
    mx, my = float(np.nanmean(xs)), float(np.nanmean(ys))
    cov = np.cov(np.vstack([xs, ys]), bias=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]; eigvecs = eigvecs[:, order]
    chi2_map = {0.50: 1.3863, 0.68: 2.279, 0.90: 4.605, 0.95: 5.991, 0.99: 9.210}
    chi2_val = chi2_map.get(round(confidence, 2), 5.991)
    scale = np.sqrt(chi2_val)
    width  = 2.0 * scale * np.sqrt(max(eigvals[0], 0.0))
    height = 2.0 * scale * np.sqrt(max(eigvals[1], 0.0))
    angle = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    ell = Ellipse((mx, my), width=width, height=height, angle=angle,
                  edgecolor=color, facecolor='none', linewidth=lw, alpha=alpha,
                  zorder=zorder)
    ax.add_patch(ell)
    return ell


def detect_palm_phases_for_bout(xs: np.ndarray, ys: np.ndarray, s: int, e: int,
                                sx: float, sy: float, fps: float,
                                radius_x: float, speed_thr: float, min_frames: int):
    """Detect advance / grasp / retract for a palm bout and compute grasp distances."""
    idx = np.arange(s, e + 1)
    x = xs[s:e+1]; y = ys[s:e+1]
    if x.size == 0:
        return None

    d_x = np.abs(x - sx)
    i_rel_min = int(np.nanargmin(d_x))
    i_min = s + i_rel_min

    step_speed = np.hypot(np.diff(x), np.diff(y)) * float(fps)
    speed = np.empty_like(x)
    if step_speed.size:
        speed[0] = step_speed[0]
        speed[1:] = step_speed
    else:
        speed[:] = 0.0

    near = d_x <= float(radius_x)
    slow = speed <= float(speed_thr)
    grasp_cand = near & slow

    def contiguous_runs(mask_bool):
        ii = np.where(mask_bool)[0]
        if ii.size == 0: return []
        runs, a, prev = [], ii[0], ii[0]
        for k in ii[1:]:
            if k != prev + 1:
                runs.append((a, prev))
                a = k
            prev = k
        runs.append((a, prev))
        return runs

    runs = contiguous_runs(grasp_cand)
    if runs:
        containing = [r for r in runs if r[0] <= i_rel_min <= r[1]]
        if containing:
            g0_rel, g1_rel = containing[0]
        else:
            centers = [(((r[0] + r[1]) // 2), r) for r in runs]
            g0_rel, g1_rel = min(centers, key=lambda t: abs(t[0] - i_rel_min))[1]
    else:
        g0_rel = g1_rel = i_rel_min

    if (g1_rel - g0_rel + 1) < int(min_frames):
        need = int(min_frames) - (g1_rel - g0_rel + 1)
        left = need // 2
        right = need - left
        g0_rel = max(0, g0_rel - left)
        g1_rel = min(len(idx) - 1, g1_rel + right)

    grasp_start = s + g0_rel
    grasp_end   = s + g1_rel

    adv_start, adv_end = s, max(s, grasp_start - 1)
    ret_start, ret_end = min(e, grasp_end + 1), e

    def seg_len_dur(a: int, b: int):
        if b - a + 1 < 2:
            return 0.0, 0.0
        xx = xs[a:b+1]; yy = ys[a:b+1]
        L = float(np.nansum(np.hypot(np.diff(xx), np.diff(yy))))
        T = (b - a) / float(fps)
        return L, T

    adv_len, adv_dur = seg_len_dur(adv_start, adv_end)
    gra_len, gra_dur = seg_len_dur(grasp_start, grasp_end)
    ret_len, ret_dur = seg_len_dur(ret_start, ret_end)
    bout_len, bout_dur = seg_len_dur(s, e)

    gx = xs[grasp_start:grasp_end+1]
    gy = ys[grasp_start:grasp_end+1]
    gdist = np.hypot(gx - sx, gy - sy)
    grasp_mean_dist = float(np.nanmean(gdist)) if gdist.size else float("nan")
    grasp_min_dist  = float(np.nanmin(gdist)) if gdist.size else float("nan")

    return {
        "start_frame": s, "end_frame": e, "closest_frame": i_min,
        "grasp_start_frame": grasp_start, "grasp_end_frame": grasp_end,
        "advance_length_px": adv_len, "advance_duration_s": adv_dur,
        "grasp_length_px": gra_len,   "grasp_duration_s": gra_dur,
        "retract_length_px": ret_len, "retract_duration_s": ret_dur,
        "bout_length_px": bout_len,   "bout_duration_s": bout_dur,
        "grasp_mean_dist_to_spout_px": grasp_mean_dist,
        "grasp_min_dist_to_spout_px":  grasp_min_dist,
    }


def process_file(csv_path: Path) -> Dict[str, Optional[float]]:
    df = load_dlc_multiindex(csv_path)
    scorer = get_main_scorer(df)

    nose_bp   = pick_bodypart(df, [NOSE_NAME, "snout", "nose"])
    palm_bp   = pick_bodypart(df, [PALM_NAME, "palm", "forepaw", "wrist"])
    spout_bp  = pick_bodypart(df, [SPOUT_NAME, "spout", "lickspout", "nozzle", "water"])
    mouth_bp = pick_bodypart(df, [MOUTH_NAME, "mouth", "mouth_tip", "mouthtip"])  # <— NEW
    if nose_bp is None or palm_bp is None or spout_bp is None:
        missing = []
        if nose_bp is None: missing.append("snout/nose")
        if palm_bp is None: missing.append("palm/forepaw/wrist")
        if spout_bp is None: missing.append("spout")
        raise RuntimeError(f"Missing markers in {csv_path.name}: {', '.join(missing)}")

    nx, ny, nl = extract_xyz(df, scorer, nose_bp)
    px, py, pl = extract_xyz(df, scorer, palm_bp)
    sx_series, sy_series, sl = extract_xyz(df, scorer, spout_bp)
    if mouth_bp is not None:
        tx, ty, tl = extract_xyz(df, scorer, mouth_bp)
    else:
        tx = ty = tl = None

    # spout median coordinate (prefer confident frames)
    if np.sum(sl >= LIKELIHOOD_THR) > 0:
        sx = float(np.nanmedian(sx_series[sl >= LIKELIHOOD_THR]))
        sy = float(np.nanmedian(sy_series[sl >= LIKELIHOOD_THR]))
    else:
        sx = float(np.nanmedian(sx_series))
        sy = float(np.nanmedian(sy_series))

    # masks: likelihood + horizontal proximity to spout
    nose_mask = (nl >= LIKELIHOOD_THR) & (np.abs(nx - sx) <= float(MAX_DX_FROM_SPOUT))
    palm_mask = (pl >= LIKELIHOOD_THR) & (np.abs(px - sx) <= float(MAX_DX_FROM_SPOUT))
    if mouth_bp is not None:
        mouth_mask = (tl >= LIKELIHOOD_THR) & (np.abs(tx - sx) <= float(MAX_DX_FROM_SPOUT))
    else:
        mouth_mask = None

    # --- Exclude PALM bouts with x-range < MIN_PALM_BOUT_XRANGE_PX ---
    palm_mask_filt = palm_mask.copy()
    for s, e in find_contiguous_segments(palm_mask_filt, MIN_SEG_LEN):
        xs = px[s:e+1]; xs = xs[np.isfinite(xs)]
        xrng = (np.nanmax(xs) - np.nanmin(xs)) if xs.size else 0.0
        if xrng < float(MIN_PALM_BOUT_XRANGE_PX):
            palm_mask_filt[s:e+1] = False

    # --- Exclude MOUTH bouts with x-range < MIN_MOUTH_BOUT_XRANGE_PX (if present) ---
    if mouth_mask is not None:
        mouth_mask_filt = mouth_mask.copy()
        for s, e in find_contiguous_segments(mouth_mask_filt, MIN_SEG_LEN):
            xs = tx[s:e+1]; xs = xs[np.isfinite(xs)]
            xrng = (np.nanmax(xs) - np.nanmin(xs)) if xs.size else 0.0
            if xrng < float(MIN_MOUTH_BOUT_XRANGE_PX):
                mouth_mask_filt[s:e+1] = False
    else:
        mouth_mask_filt = None

    # totals (pixels)
    nose_total  = total_path_length(nx, ny, nose_mask, MIN_SEG_LEN)
    palm_total  = total_path_length(px, py, palm_mask_filt, MIN_SEG_LEN)
    if mouth_mask_filt is not None:
        mouth_total = total_path_length(tx, ty, mouth_mask_filt, MIN_SEG_LEN)
    else:
        mouth_total = float("nan")

    # basic bout metrics
    n_palm_bouts,  avg_palm_len,  avg_palm_spd  = bout_metrics(px, py, palm_mask_filt, MIN_SEG_LEN, FRAMERATE_FPS)
    n_snout_bouts, avg_snout_len, avg_snout_spd = bout_metrics(nx, ny, nose_mask,      MIN_SEG_LEN, FRAMERATE_FPS)
    if mouth_mask_filt is not None:
        n_mouth_bouts, avg_mouth_len, avg_mouth_spd = bout_metrics(tx, ty, mouth_mask_filt, MIN_SEG_LEN, FRAMERATE_FPS)
    else:
        n_mouth_bouts = avg_mouth_len = avg_mouth_spd = float("nan")

    # palm smoothness / variance
    palm_torts = bout_tortuosity(px, py, palm_mask_filt, MIN_SEG_LEN)
    avg_palm_tortuosity = float(np.nan) if not palm_torts else float(np.nanmean(palm_torts))
    smooth_list = [1.0 / t for t in palm_torts if np.isfinite(t) and t > 0]
    avg_palm_smooth_idx = float(np.nan) if not smooth_list else float(np.nanmean(smooth_list))
    palm_traj_var = trajectory_variance_px2(px, py, palm_mask_filt)

    # mouth variance
    if mouth_mask_filt is not None:
        mouth_traj_var = trajectory_variance_px2(tx, ty, mouth_mask_filt)
    else:
        mouth_traj_var = float("nan")

    # Mean palm trajectory across bouts (after filtering)
    mt_x, mt_y = mean_trajectory_over_bouts(px, py, palm_mask_filt, MIN_SEG_LEN, n_points=MEAN_TRAJ_POINTS)

    # ====== PALM PHASES & GRASPS ======
    phase_rows = []
    for s, e in find_contiguous_segments(palm_mask_filt, MIN_SEG_LEN):
        row = detect_palm_phases_for_bout(
            px, py, s, e, sx, sy, FRAMERATE_FPS,
            radius_x=GRASP_RADIUS_X_PX,
            speed_thr=GRASP_SPEED_THR_PX_S,
            min_frames=GRASP_MIN_FRAMES
        )
        if row is not None:
            row["bout_id"] = len(phase_rows) + 1
            phase_rows.append(row)

    grasp_mask = np.zeros_like(palm_mask_filt, dtype=bool)
    for r in phase_rows:
        gs, ge = int(r["grasp_start_frame"]), int(r["grasp_end_frame"])
        grasp_mask[gs:ge+1] = True
    grasp_mask &= palm_mask_filt

    if phase_rows:
        df_phase = pd.DataFrame(phase_rows)
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        phase_csv_path = OUT_DIR / f"{csv_path.stem}_palm_bout_phases.csv"
        df_phase.to_csv(phase_csv_path, index=False)

        avg_adv_len = float(df_phase["advance_length_px"].mean())
        avg_adv_dur = float(df_phase["advance_duration_s"].mean())
        avg_gra_dur = float(df_phase["grasp_duration_s"].mean())
        avg_ret_len = float(df_phase["retract_length_px"].mean())
        avg_ret_dur = float(df_phase["retract_duration_s"].mean())
        n_grasps = int((df_phase["grasp_duration_s"] > 0).sum())
        avg_grasp_mean_dist = float(df_phase["grasp_mean_dist_to_spout_px"].mean())
        acc_ratio = float((df_phase["grasp_min_dist_to_spout_px"] <= ACCURACY_RADIUS_PX).mean())
    else:
        df_phase = None
        phase_csv_path = None
        avg_adv_len = avg_adv_dur = avg_gra_dur = avg_ret_len = avg_ret_dur = float("nan")
        n_grasps = 0
        avg_grasp_mean_dist = float("nan")
        acc_ratio = float("nan")

    # grasp positions variance (px^2)
    grasp_positions_variance = trajectory_variance_px2(px, py, grasp_mask) if np.any(grasp_mask) else float("nan")

    # ================== PLOTS ==================
    # (1) Main traces figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # Snout behind
    for s, e in find_contiguous_segments(nose_mask, MIN_SEG_LEN):
        xs, ys = nx[s:e+1], ny[s:e+1]
        ax.plot(xs, ys, linewidth=LINE_WIDTH, color=SNOUT_COLOR, alpha=SNOUT_ALPHA, zorder=1)

    # mouth (middle layer) — only if present
    if mouth_mask_filt is not None:
        for s, e in find_contiguous_segments(mouth_mask_filt, MIN_SEG_LEN):
            xs, ys = tx[s:e+1], ty[s:e+1]
            ax.plot(xs, ys, linewidth=LINE_WIDTH, color=MOUTH_COLOR, alpha=MOUTH_ALPHA, zorder=2)

    # Palm in front
    for s, e in find_contiguous_segments(palm_mask_filt, MIN_SEG_LEN):
        xs, ys = px[s:e+1], py[s:e+1]
        ax.plot(xs, ys, linewidth=LINE_WIDTH, color=PALM_COLOR, alpha=PALM_ALPHA, zorder=3)

    # Mean palm trajectory
    if mt_x is not None and mt_y is not None:
        ax.plot(mt_x, mt_y, color=MEAN_TRAJ_COLOR, linewidth=MEAN_TRAJ_WIDTH,
                alpha=MEAN_TRAJ_ALPHA, zorder=4)

    # Covariance ellipses (95%)
    add_covariance_ellipse(ax, px, py, palm_mask_filt,   ELLIPSE_CONFIDENCE, PALM_COLOR,   ELLIPSE_LW, ELLIPSE_ALPHA, 5)
    add_covariance_ellipse(ax, nx, ny, nose_mask,        ELLIPSE_CONFIDENCE, SNOUT_COLOR,  ELLIPSE_LW, ELLIPSE_ALPHA, 5)
    if mouth_mask_filt is not None:
        add_covariance_ellipse(ax, tx, ty, mouth_mask_filt, ELLIPSE_CONFIDENCE, MOUTH_COLOR, ELLIPSE_LW, ELLIPSE_ALPHA, 5)

    # Spout
    ax.scatter([sx], [sy], s=80, marker="X", color="black", zorder=6)

    title_main = (
        f"{csv_path.stem}\n"
        f"L≥{LIKELIHOOD_THR}, |x-spout_x|≤{MAX_DX_FROM_SPOUT}px | "
        f"Excluded palm x-range<{MIN_PALM_BOUT_XRANGE_PX}px"
        + (f", mouth x-range<{MIN_MOUTH_BOUT_XRANGE_PX}px" if mouth_mask_filt is not None else "")
    )
    ax.set_title(title_main)
    ax.set_xlabel("x (pixels)"); ax.set_ylabel("y (pixels)")
    ax.invert_yaxis(); ax.axis("equal"); plt.tight_layout()

    legend_elements = [
        Line2D([0],[0], color=SNOUT_COLOR, alpha=SNOUT_ALPHA, lw=2, label="Snout"),
        Line2D([0],[0], color=MOUTH_COLOR, alpha=MOUTH_ALPHA, lw=2, label="Mouth") if mouth_mask_filt is not None else None,
        Line2D([0],[0], color=PALM_COLOR, alpha=PALM_ALPHA, lw=2, label="Palm"),
        Line2D([0],[0], color=MEAN_TRAJ_COLOR, alpha=MEAN_TRAJ_ALPHA, lw=MEAN_TRAJ_WIDTH, label="Palm mean trajectory"),
        Line2D([0],[0], marker="X", color="black", linestyle="None", label="Spout (median)")
    ]
    legend_elements = [h for h in legend_elements if h is not None]
    ax.legend(handles=legend_elements, loc="best")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = OUT_DIR / f"{csv_path.stem}_palm_snout_mouth.pdf"
    plt.savefig(fig_path, bbox_inches="tight")
    if SHOW_PLOTS: plt.show()
    plt.close(fig)

    # (2) Grasps-only figure (dots + covariance ellipse for grasps)
    fig2, ax2 = plt.subplots(figsize=(7, 7))
    if np.any(grasp_mask):
        ax2.scatter(px[grasp_mask], py[grasp_mask], s=GRASP_DOT_SIZE, color=GRASP_COLOR, alpha=GRASP_ALPHA, zorder=1)
        add_covariance_ellipse(ax2, px, py, grasp_mask, ELLIPSE_CONFIDENCE, GRASP_COLOR, ELLIPSE_LW, 1.0, 2)
        n_pts = int(np.sum(np.isfinite(px[grasp_mask]) & np.isfinite(py[grasp_mask])))
        subtitle = f"{n_pts} grasp points; Var={grasp_positions_variance:.2f} px²"
    else:
        subtitle = "No grasps detected"
    ax2.scatter([sx], [sy], s=80, marker="X", color="black", zorder=3)
    ax2.set_title(f"{csv_path.stem} — Grasps (dots) + covariance ellipse\n{subtitle}")
    ax2.set_xlabel("x (pixels)"); ax2.set_ylabel("y (pixels)")
    ax2.invert_yaxis(); ax2.axis("equal"); plt.tight_layout()

    grasps_fig_path = OUT_DIR / f"{csv_path.stem}_grasps_only.pdf"
    plt.savefig(grasps_fig_path, bbox_inches="tight")
    if SHOW_PLOTS: plt.show()
    plt.close(fig2)
    # =========================================

    return {
        "mouse_id": csv_path.stem,
        # totals & basic bout metrics
        "total_distance_snout_px": nose_total,
        "total_distance_palm_px":  palm_total,
        "total_distance_mouth_px": mouth_total,
        "n_snout_bouts":  n_snout_bouts,
        "avg_snout_bout_length_px": avg_snout_len,
        "avg_snout_bout_speed_px_per_s": avg_snout_spd,
        "n_palm_bouts":  n_palm_bouts,
        "avg_palm_bout_length_px":  avg_palm_len,
        "avg_palm_bout_speed_px_per_s": avg_palm_spd,
        "n_mouth_bouts": n_mouth_bouts,
        "avg_mouth_bout_length_px": avg_mouth_len,
        "avg_mouth_bout_speed_px_per_s": avg_mouth_spd,
        # smoothness & variance (palm/mouth)
        "avg_palm_bout_tortuosity": avg_palm_tortuosity,
        "avg_palm_bout_smoothness_index": avg_palm_smooth_idx,
        "palm_trajectory_variance_px2":   palm_traj_var,
        "mouth_trajectory_variance_px2": mouth_traj_var,
        # palm phase & grasp metrics
        "n_palm_grasps_detected": n_grasps,
        "avg_advance_length_px":  avg_adv_len,
        "avg_advance_duration_s": avg_adv_dur,
        "avg_grasp_duration_s":   avg_gra_dur,
        "avg_retract_length_px":  avg_ret_len,
        "avg_retract_duration_s": avg_ret_dur,
        "avg_grasp_mean_dist_to_spout_px": avg_grasp_mean_dist,
        "grasp_accuracy_radius_px": ACCURACY_RADIUS_PX,
        "grasp_accuracy_ratio":    acc_ratio,
        # outputs
        "figure_path": str(fig_path),
        "grasps_figure_path": str(grasps_fig_path),
        "phase_table_path": str(phase_csv_path) if phase_csv_path else "",
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    files = sorted(DATA_DIR.glob(CSV_PATTERN))
    if not files:
        print(f"No files found in {DATA_DIR} matching '{CSV_PATTERN}'.")
        return

    rows, errors = [], []
    for f in files:
        try:
            print(f"Processing: {f.name}")
            rows.append(process_file(f))
        except Exception as e:
            print(f"  ⚠️ Skipping {f.name}: {e}")
            errors.append({"file": f.name, "error": str(e)})

    if rows:
        summary = pd.DataFrame(rows).sort_values("mouse_id")
        csv_path = OUT_DIR / "trace_distance_summary.csv"
        summary.to_csv(csv_path, index=False)
        if errors:
            pd.DataFrame(errors).to_csv(OUT_DIR / "errors.csv", index=False)
        print(f"\nSummary written to: {csv_path}")
    else:
        print("\nNo results to write.")

    if errors:
        print("\nFiles with issues:")
        for e in errors:
            print(f" - {e['file']}: {e['error']}")


if __name__ == "__main__":
    main()
