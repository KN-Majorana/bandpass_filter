# 干渉バンドパスフィルターの再現 — 転送行列法によるミニ検証

設計対象は**全誘電体 Fabry-Pérot 型干渉バンドパスフィルター**で、市販品
Thorlabs FB550-10（中心波長 550 nm, FWHM ≈ 10 nm）相当の性能を再現するとともに、
**「FWHM を狭くするほどピーク透過率が落ちる」というトレードオフ**が
材料の吸収（消衰係数 k）にどう支配されるかを定量的に示した。

## 設計の概要

```
Air | (L H)^N | 2L (cavity) | (H L)^N | Glass
```

- H : TiO₂（高屈折率, n ≈ 2.65 @ 550 nm）厚さ λ₀/(4 n_H) ≈ 52 nm
- L : SiO₂（低屈折率, n ≈ 1.46 @ 550 nm）厚さ λ₀/(4 n_L) ≈ 94 nm
- 2L : 半波長スペーサ（共振器）厚さ λ₀/(2 n_L) ≈ 188 nm

両側の (LH)^N / (HL)^N が分布ブラッグ反射器として高反射ミラーを形成し、
中央の 2L 層がそのキャビティを構成する Fabry-Pérot 共振器。

## 計算結果

### 1. 層構造の断面図（geometry.png）

設計波長 550 nm に対して、各層厚を λ/4 で量子化した多層膜の断面。

### 2. 透過スペクトル（transmittance.png）

550 nm でピーク透過率 ≈ 96 %、FWHM ≈ 1 nm のバンドパス特性を確認。
エネルギー保存則 |1 − R − T| < 1e−10 で TMM の数値精度を担保。

### 3. 入射角依存性（angle_dependence.png）

斜入射ではピーク波長が短波長側にシフトする（角度シフト）。
干渉条件 2 n d cos θ = m λ から定量的に予測される挙動と一致。

### 4. FWHM と ピーク透過率のトレードオフ（fwhm_vs_tpeak_tradeoff.png）

**本検証の主結果。**ミラー対数 N を 2〜12 で振り、材料の消衰係数 k_H を 4 通り
（0, 1e-4, 5e-4, 2e-3）で変えて FWHM と T_peak を測定し、パレートフロントを描いた。

- **k = 0（無損失・理想）**: T_peak が一定のまま FWHM だけが単調に狭くなる。
- **k > 0（現実材料）**: 共振器内の往復回数（フィネス）が上がるほど、わずかな
  吸収が指数的に効き、T_peak が急落する。
- 同じ N でも、成膜プロセスの質（k 値）で性能が大きく変わることが示される。

これが「**狭く・明るく**は同時に追求できない」という、材料開発における典型的な
多次元最適化問題の構造である。

## ファイル

| ファイル | 内容 |
|---|---|
| `bandpass_filter.py` | フィルター設計と図1〜3（断面・透過スペクトル・角度依存性）の生成 |
| `fwhm_tpeak_tradeoff.py` | 図4（FWHM vs T_peak トレードオフ）の生成 |
| `geometry.png` | 多層膜の断面構造 |
| `transmittance.png` | 透過スペクトル（中心 550 nm, FWHM ≈ 1 nm） |
| `angle_dependence.png` | 入射角を変えたときの透過スペクトル変化 |
| `fwhm_vs_tpeak_tradeoff.png` | FWHM と ピーク透過率のパレートフロント |

## 依存関係

このスクリプト群は、親ディレクトリの `src/` モジュールに依存する。

- `src/tmm.py` — 転送行列法のコア実装（Macleod 規約に従う）
- `src/materials.py` — Sellmeier 分散モデルによる n(λ) — TiO₂, SiO₂, BK7

### 必要パッケージ

```
numpy >= 1.24
matplotlib >= 3.7
japanize-matplotlib   # 日本語ラベル用
```

## 実行方法

```bash
cd fujifilm/ミニ検証
python bandpass_filter.py        # 図1〜3 を生成
python fwhm_tpeak_tradeoff.py    # 図4 を生成
```

## 参考文献

- H. A. Macleod, *Thin-Film Optical Filters*, 5th ed., CRC Press (2018), Ch. 2 & 7.
- Thorlabs FB550-10 データシート（ベンチマーク）
