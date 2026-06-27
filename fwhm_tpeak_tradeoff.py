"""
07: FWHM と ピーク透過率のトレードオフ解析
====================================================================

干渉バンドパスフィルター（Fabry-Perot 型）で、ミラー対数 N を増やすと
帯域幅 FWHM は狭くなるが、現実材料には僅かな吸収（消衰係数 k > 0）が
あるため、共振器内の往復回数増加とともにロスが増幅され、
ピーク透過率 T_peak が下がる。

これを示すために:
  - N_PAIRS を 2..12 で変化
  - 高屈折率層 TiO2 の k_H を 4 通り (0, 1e-4, 5e-4, 2e-3)
の組み合わせで FWHM と T_peak を計算し、パレートフロントを描く。

Macleod 規約: 吸収体は N = n - i*k (k > 0)。
k = 0 では T_peak ~ 100% を保ったまま FWHM が単調に狭くなる。
k > 0 では「狭く・明るく」を同時に追求できないトレードオフが出る。
これが「材料開発 = 多次元最適化問題」の典型的な構図。
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import matplotlib.pyplot as plt
import japanize_matplotlib  # 日本語ラベルを有効化

from src.tmm import transmittance_unpolarized
from src.materials import n_tio2, n_sio2, n_glass_BK7

from bandpass_filter import build_stack, peak_and_fwhm

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = HERE


def complex_layer_idx(materials, wl, k_H=0.0, k_L=0.0):
    """Macleod 規約 N = n - i*k (k > 0) で吸収を入れた屈折率配列。"""
    n_table = {
        "TiO2": n_tio2(wl) - 1j * k_H,
        "SiO2": n_sio2(wl) - 1j * k_L,
    }
    return np.stack([n_table[m] for m in materials], axis=1)


def main():
    # 550 nm 近傍を細かくサンプリング（高フィネス時の鋭いピーク用）。
    wl_coarse = np.linspace(450, 700, 2501)
    wl_fine   = np.linspace(540, 560, 20001)  # 0.001 nm 刻み
    wl = np.unique(np.concatenate([wl_coarse, wl_fine]))

    N_list = list(range(2, 7))
    k_levels = [0.0, 1e-4, 5e-4, 2e-3]
    labels = {
        0.0:  "k = 0 (無損失・理想)",
        1e-4: "k = 1e-4 (高品質 TiO$_2$)",
        5e-4: "k = 5e-4 (標準品質 TiO$_2$)",
        2e-3: "k = 2e-3 (低品質 TiO$_2$)",
    }
    colors  = ["#444444", "#1f77b4", "#ff7f0e", "#d62728"]
    markers = ["o", "s", "^", "D"]

    results = {}
    for k_H in k_levels:
        fwhms, tpeaks = [], []
        for N in N_list:
            mats, ds = build_stack(N)
            idx = complex_layer_idx(mats, wl, k_H=k_H, k_L=0.0)
            T = transmittance_unpolarized(
                wl, idx, np.array(ds), n_substrate=n_glass_BK7(wl)
            )
            _, T_peak, fwhm = peak_and_fwhm(wl, T)
            fwhms.append(fwhm)
            tpeaks.append(T_peak * 100.0)
        results[k_H] = (fwhms, tpeaks)
        print(f"k_H = {k_H:.0e}:")
        for N, f, t in zip(N_list, fwhms, tpeaks):
            print(f"   N={N:2d}  FWHM={f:8.3f} nm   T_peak={t:6.2f} %")

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for (k_H, (fwhms, tpeaks)), c, m in zip(results.items(), colors, markers):
        ax.plot(fwhms, tpeaks, marker=m, color=c, lw=1.6, ms=7,
                label=labels[k_H])

    # 最も損失のあるカーブに N のラベル（赤線は左上→右下方向のため各点の左上に配置）
    fwhms_lossy, tpeaks_lossy = results[k_levels[-1]]
    # N=2 は最左端なので右上、他は左上にオフセット
    N_offsets = {2: (10, 6), 3: (-8, 8), 4: (-8, 8), 5: (-8, 8), 6: (-8, 8)}
    bbox_kw = dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85)
    for N, f, t in zip(N_list, fwhms_lossy, tpeaks_lossy):
        dx, dy = N_offsets.get(N, (-8, 8))
        ax.annotate(f"N={N}", xy=(f, t), xytext=(dx, dy),
                    textcoords="offset points", fontsize=10, color="#d62728",
                    bbox=bbox_kw)

    # 「右上が理想方向」を対角矢印＋星マーカーで一目で示す
    # 軸は invert_xaxis により右=狭い(FWHM小)、上=明るい(T_peak高) → 右上が理想
    ax.annotate("",
                xy=(0.94, 0.93), xycoords="axes fraction",
                xytext=(0.28, 0.18), textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#2ca02c",
                                lw=2.0, mutation_scale=20, alpha=0.5))
    ax.text(0.45, 0.45, "理想の方向", transform=ax.transAxes,
            ha="center", va="center", fontsize=14, color="#2ca02c",
            fontweight="bold", rotation=40,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="none", alpha=0.85))


    ax.set_xscale("log")
    ax.set_xlabel("半値全幅 FWHM [nm]  (狭くなる →)", fontsize=12)
    ax.set_ylabel("ピーク透過率 [%]", fontsize=12)
    ax.set_title("半値全幅 vs ピーク透過率,　(ミラー対数 N = 2〜7 を変化)", fontsize=13)
    ax.set_ylim(0, 105)
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower left", fontsize=11, title="材料の吸収レベル")
    ax.invert_xaxis()
    fig.tight_layout()

    out = os.path.join(FIGDIR, "fwhm_vs_tpeak_tradeoff.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()