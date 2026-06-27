"""
06: 干渉バンドパスフィルター（特定波長のみを透過させる材料）の再現
====================================================================

「特定波長のみを透過させる材料」として、最も古典的・かつ実用化されているのが
**全誘電体 Fabry-Pérot 型干渉バンドパスフィルター**。

構造（Macleod, "Thin-Film Optical Filters" 5th ed., Ch.7）:

    Air | (L H)^N | 2L (cavity) | (H L)^N | Glass

  H : 高屈折率層 (TiO₂, n ≈ 2.35-2.6)   厚さ λ₀/(4 n_H)
  L : 低屈折率層 (SiO₂, n ≈ 1.46)        厚さ λ₀/(4 n_L)
  2L: 半波長スペーサ（=共振器）          厚さ λ₀/(2 n_L)

中央付近が "...L H | 2L | H L..." となり、両脇の H 層がキャビティを挟む
部分鏡として働く Fabry-Pérot 共振器を形成。両側の Bragg ミラーが λ₀ で
高反射する一方、共振器が λ₀ でのみ位相条件を満たすため、その波長だけ透過する。

ベンチマーク（市販品例）:
    Thorlabs FB550-10  (CWL = 550 nm, FWHM ≈ 10 nm)
    Edmund Optics 65-687 (550 nm bandpass)
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import japanize_matplotlib 

from src.tmm import transmittance_unpolarized, reflectance_transmittance
from src.materials import n_tio2, n_sio2, n_glass_BK7

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = HERE

# ============================================================================
# 設計
# ============================================================================
LAMBDA0_NM = 550.0
N_PAIRS    = 5

n_H = float(n_tio2(np.array([LAMBDA0_NM]))[0])
n_L = float(n_sio2(np.array([LAMBDA0_NM]))[0])

d_H = LAMBDA0_NM / (4 * n_H)
d_L = LAMBDA0_NM / (4 * n_L)
d_spacer = LAMBDA0_NM / (2 * n_L)


def build_stack(n_pairs: int = N_PAIRS):
    """Air | (L H)^n_pairs | spacer(2L) | (H L)^n_pairs | Glass の順で返す。"""
    mats, ds = [], []
    for _ in range(n_pairs):
        mats += ["SiO2", "TiO2"]
        ds += [d_L, d_H]
    mats.append("SiO2")
    ds.append(d_spacer)
    for _ in range(n_pairs):
        mats += ["TiO2", "SiO2"]
        ds += [d_H, d_L]
    return mats, ds


def peak_and_fwhm(wl, T, wl_lo=500.0, wl_hi=600.0):
    """指定波長域の局所ピークから (lam_peak, T_peak, FWHM) を返す。"""
    mask = (wl >= wl_lo) & (wl <= wl_hi)
    i_local = int(np.argmax(T[mask]))
    i_peak = int(np.where(mask)[0][i_local])
    T_peak = float(T[i_peak])
    half = T_peak / 2.0
    li = i_peak
    while li > 0 and T[li] >= half:
        li -= 1
    ri = i_peak
    while ri < len(wl) - 1 and T[ri] >= half:
        ri += 1
    fwhm = float(wl[ri] - wl[li])
    return float(wl[i_peak]), T_peak, fwhm


# ============================================================================
def main():
    materials, thicknesses = build_stack()
    L = len(materials)
    print(f"Layer count: {L}")
    print(f"Design wavelength: {LAMBDA0_NM} nm")
    print(f"n_H (TiO2) = {n_H:.3f}, n_L (SiO2) = {n_L:.3f}")
    print(f"d_H = {d_H:.2f} nm, d_L = {d_L:.2f} nm, spacer = {d_spacer:.2f} nm")

    wl = np.linspace(450, 700, 4001)
    n_table = {"TiO2": n_tio2(wl), "SiO2": n_sio2(wl)}
    layer_idx = np.stack([n_table[m] for m in materials], axis=1)
    d = np.array(thicknesses)

    T = transmittance_unpolarized(wl, layer_idx, d, n_substrate=n_glass_BK7(wl))

    lam_peak, T_peak, fwhm = peak_and_fwhm(wl, T)
    print(f"Peak: lam = {lam_peak:.2f} nm, T = {T_peak*100:.2f} %, FWHM = {fwhm:.2f} nm")

    R_s, T_s = reflectance_transmittance(wl, layer_idx, d,
                                         n_substrate=n_glass_BK7(wl), polarization="s")
    print(f"Energy conservation check  max|1-R-T| = {np.max(np.abs(1 - R_s - T_s)):.2e}")

    # ========================================================================
    # 図1: 層構造の断面図
    # ========================================================================
    fig1, ax = plt.subplots(figsize=(7, 3.5))
    color_map = {"TiO2": "#d62728", "SiO2": "#1f77b4"}
    sub_h = 80
    x = 0.0
    ax.add_patch(Rectangle((-sub_h, 0), sub_h, 1, facecolor="#cccccc",
                            edgecolor="black", linewidth=0.5))
    ax.text(-sub_h / 2, -0.1, "ガラス基板", ha="center", va="top", fontsize=10)

    spacer_index = N_PAIRS * 2
    for i, (m, dj) in enumerate(zip(materials, thicknesses)):
        is_spacer = (i == spacer_index)
        lw = 1.2 if is_spacer else 0.4
        ax.add_patch(Rectangle((x, 0), dj, 1, facecolor=color_map[m],
                                edgecolor="black", linewidth=lw, alpha=0.85))
        if is_spacer:
            ax.annotate(
                f"(SiO2, {dj:.1f} nm)",
                xy=(x + dj / 2, 1), xytext=(x + dj / 2, 1.28),
                ha="center", va="bottom", fontsize=9, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
            )
        x += dj

    air_h = 140
    ax.add_patch(Rectangle((x, 0), air_h, 1, facecolor="white",
                            edgecolor="black", linewidth=0.5))
    ax.text(x + air_h / 2, -0.1, "空気", ha="center", va="top", fontsize=10)

    ax.annotate("", xy=(x + 4, 0.5), xytext=(x + air_h - 4, 0.5),
                arrowprops=dict(arrowstyle="->", lw=2, color="orange"))
    ax.text(x + air_h / 2, 0.40, "白色光",
            color="orange", fontsize=9, fontweight="bold", ha="center", va="top")

    legend_elements = [
        Patch(facecolor=color_map["TiO2"], edgecolor="black",
              label=f"H = TiO2  ({d_H:.1f} nm,  n={n_H:.2f})"),
        Patch(facecolor=color_map["SiO2"], edgecolor="black",
              label=f"L = SiO2 ({d_L:.1f} nm,  n={n_L:.2f})"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9)

    ax.set_xlim(-sub_h - 10, x + air_h + 10)
    ax.set_ylim(-0.25, 1.55)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("ガラス面からの累積膜厚 [nm]")
    ax.set_title(f"Fabry-Pérot バンドパスフィルター  (LH)$^{{{N_PAIRS}}}$ · 2L · (HL)$^{{{N_PAIRS}}}$",
                 fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig1.tight_layout()
    out1 = os.path.join(FIGDIR, "geometry.png")
    fig1.savefig(out1, dpi=140, bbox_inches="tight")
    print(f"Saved: {out1}")

    # ========================================================================
    # 図2: 透過スペクトル
    # ========================================================================
    fig2, ax2 = plt.subplots(figsize=(5.0, 4.8))
    ax2.plot(wl, T * 100, lw=1.8, color="C3")
    ax2.fill_between(wl, 0, T * 100, color="C3", alpha=0.15)
    ax2.axvline(LAMBDA0_NM, color="k", ls=":", lw=0.8,
                label=f"設計波長 = {LAMBDA0_NM:.0f} nm")
    ax2.annotate(
        f"ピーク\nλ = {lam_peak:.1f} nm\nT = {T_peak*100:.1f}%",
        xy=(lam_peak, T_peak * 100),
        xytext=(lam_peak + 60, T_peak * 100 - 22),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="C3"),
        arrowprops=dict(arrowstyle="->", color="C3"),
    )
    ax2.set_xlabel("波長 λ [nm]")
    ax2.set_ylabel("透過率 T [%]（非偏光・垂直入射）")
    ax2.set_title("550 nm 付近を透過するバンドパスフィルター")
    ax2.set_xlim(wl.min(), wl.max())
    ax2.set_ylim(0, 105)
    ax2.legend(loc="upper right")
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    out2 = os.path.join(FIGDIR, "transmittance.png")
    fig2.savefig(out2, dpi=140, bbox_inches="tight")
    print(f"Saved: {out2}")

    # ========================================================================
    # 図3: 入射角依存性
    # ========================================================================
    angles_deg = [0, 15, 30, 45, 60]
    cmap = plt.get_cmap("viridis")
    n_rest = len(angles_deg) - 1
    colors = ["C3"] + [cmap(i / (n_rest - 1)) for i in range(n_rest)]

    fig3, ax3 = plt.subplots(figsize=(5.0, 4.8))
    for ang, col in zip(angles_deg, colors):
        T_ang = transmittance_unpolarized(
            wl, layer_idx, d,
            n_substrate=n_glass_BK7(wl),
            theta0=np.deg2rad(ang),
        )
        ax3.plot(wl, T_ang * 100, lw=1.8, color=col, label=f"θ = {ang}°")

    ax3.set_xlabel("波長 λ [nm]")
    ax3.set_ylabel("透過率 T [%]（非偏光）")
    ax3.set_title("バンドパスフィルターの入射角依存性")
    ax3.set_xlim(450, 700)
    ax3.set_ylim(0, 105)
    ax3.legend(loc="upper right")
    ax3.grid(alpha=0.3)
    fig3.tight_layout()
    out3 = os.path.join(FIGDIR, "angle_dependence.png")
    fig3.savefig(out3, dpi=140, bbox_inches="tight")
    print(f"Saved: {out3}")


if __name__ == "__main__":
    main()
