# explore/：探索用スクリプト

浮動小数による見積もりと実験。ここの出力は証明ではない。厳密な版は `composite/` と `B72/`。

| ファイル | 内容 |
|---|---|
| `feas.py`, `feas2.py`, `feas3.py` | Miller の補題に必要な範囲 $X$ の見積もり（素数冪層、GRH 版との比較、零点のない帯の仮定ごとの表） |
| `refined.py` | $\zeta_K$ の零点を差し引く改良の評価。得にならないことの確認 |
| `compos2.py`, `comp_run.py` | 合成層の上界（浮動小数版）と証明書探索。`composite/` の元になったもの |
| `genus_scan.py` | Chevalley の公式による $p\mid h$ の判定を $p<400$ で走査（$p=31,73$ を再発見） |
| `unit_barrier.py` | 単数を使う方法の閾値（Gauss 発見的評価、Blichfeldt、Minkowski） |
| `class65_test.py` | $\ell\equiv65 \pmod{128}$ の素数ごとの格子判定 |
| `weber_wieferich.py` | 必要条件 $\ell\mid B_{1,\psi\omega^{-1}}$ が成り立つ $(\ell,\psi)$ の探索 |
| `norm_yield.py` | 小さい元 $c+\sum a_i\eta_i$（$\eta_i$ は Gauss 周期）のノルムの大きさと、小さい分解素数だけで割り切れる割合の測定。使い方：`python3 norm_yield.py "[(7,2,40000),(3,4,30000)]"` |
| `pari_probe.py` | PARI の `bnfinit` が $k_p$ で通るかの試験。$k_{37}$ は 2 コアの環境で 9 分半では終わらなかった |
| `relations_k49/` | 小さい元の関係式から素イデアルの単項性を導く試作（$k_{7^2}$）。詳細はその中の README |
| `fig8_tower.py` | 8 の字結び目のオービフォールド塔の計算 |

`comp_run.py` と `genus_scan.py` は `compos2.py` と `feas.py` を import するので、このディレクトリで実行する。
