# 再監査報告：`relations/` の修正（前回の監査 `docs/audit_relations.md` §7 への対応）

日付：2026-11-24。監査者のコードと記録はすべて `/home/claude/rel_audit2/`（リポジトリ外）にある。リポジトリ内で作ったのはこのファイルだけである。既存ファイルの変更、commit、push、Actions の起動、release の変更はしていない（`relations/__pycache__/` に Python のキャッシュができたが、`.gitignore` の対象である）。

用語：**行**は証明書ファイルの 1 行（ある完全分解素数 $q$ について、ノルムが $q$ と小さい素数の積になる元を書いたもの）。**要約**は `verify.py --range` が合格したときに出力する JSON。**検証器**は行を検査して合格・不合格を出すプログラム（`verify.py`、`audit/relations_recheck.py`）。

## 0. 結論

- **修正後のコードに、証明を壊す誤りは見つからなかった。** `up`／`down` は整数演算だけで安全な向きに丸めている。非有限の値では例外で止まる（安全な向き）。`aggregate.py` は和と個数を要約から独立に計算し直しており、前回の 3 つの改竄はすべて落ちる。`subfields.py` の論理と数値も正しい。`relations_recheck.py` は `verify.py --range` と同じ主張を検査し、改竄 34 種類のうち通してはいけないものはすべて落ちた。
- `python3 -O aggregate.py k3_4 2` の出力は `relations/k3_4/aggregate_output.txt` と完全に一致した（2 分 14 秒）。$h<99.8767<163$。
- **直すべき点は 7 つ（すべて軽微〜中程度。重大度の順）**。主なものは次の 3 つ。
  1. （中）`recheck.yml` の `collect` は、一部の区間の再検証が失敗・欠落しても部分的な結果を commit し、失敗として扱わない。結果の JSON を検査するプログラムもない。行への信頼を取り除くのがこのワークフローの目的なので、「82 個すべてが両方の検証器で合格」を機械的に確かめる手段が要る。
  2. （低〜中）`relations.yml` の `collect` は `aggregate.py ... || true` で、集計の失敗を隠す。`aggregate.py` 自身も、上界が閾値を下回らないとき（`VERDICT: bound not below`）に終了コード 0 を返す。
  3. （低）`aggregate.py` は体を要約から取る。体 $(3,4)$ と終点 $3\cdot10^{10}$ を引数で固定していない。偽の結論を出す経路ではない（出力に体が表示される）が、「k_{3^4} の証明」としては固定すべきである。
- 行の中身を全行独立に検証した範囲ファイルは、前回の 2 個に今回の 1 個を加えて 3 個になった。残り 79 個は、引き続き Actions 上での `verify.py` の合格を信頼している（§4）。

## 1. コードの監査

### 1.1 `verify.py` の `_frac`、`_dec`、`up`、`down`

- **仕様の確認（python-flint 0.9.0、手元の docstring）。** `arb.upper()` は「上界（$+\infty$ 向き）。現在の精度に**切り上げた**正確な浮動小数を持つ arb」を返す。`lower()` はその逆である。`man_exp()` は「正確かつ有限な値に対して、整数の仮数と指数の組を返す。そうでなければ `ValueError`」。
- **論理。** `_frac` は $m\cdot2^e$ を有理数 $a/b$（$b>0$）として正確に返す。`up` は $\lfloor a\cdot10^k/b\rfloor+1>a\cdot10^k/b$、`down` は $\lceil a\cdot10^k/b\rceil-1<a\cdot10^k/b$ を使う（Python の `//` は負の数でも床関数なので、`-((-x)//b)` は天井関数になる）。よって $\mathrm{up}\cdot10^k>\mathrm{upper}\ge$ 球のすべての点、$\mathrm{down}$ も同様で、どちらも**狭義**の不等式になる。`_dec` は符号と 0 埋めを正しく扱う（$-5\to$ `-0.005`、$0\to$ `0.000`）。
- **試験**（`t_updown.py`）。$0,\pm1,\pm0.1,\pm10^{-40},\pm2^{-1000},10^{50}$、半径 $10^{-3}$ の球、$0\pm10^{10}$、$-7\pm100$ について、$k=1,4,7,10,13$ で「$\mathrm{up}>$ 上端、$\mathrm{down}<$ 下端、ずれは $10^{-k}$ 以下」を有理数で確かめた。すべて成り立った。float、int、文字列の入力でも正しかった。
- **非有限の値。** `inf`、`-inf`、`nan`、および $1/(0\pm1)$（結果は nan）では、`man_exp` が `ValueError` を出して止まる。出力も合格表示も出ないので、**安全な向きに失敗する**。
- **欠陥（軽微）。** $k=0$ のとき `_dec` は誤った文字列を返す（`s[:-0]` が空になる。`up(5.3,0)` は `.6`、`down(5.3,0)` は `.5` で、どちらも値として誤り。`up` の方は**危険な向き**）。現在のコードは $k\in\{4,7,10,13\}$ でしか呼ばないので影響はない。`need(k>=1)` を入れるべきである。

### 1.2 `aggregate.py`

- **恒等式。** $\cosh(\tfrac{m}{2}\log q)=\tfrac12(q^{m/2}+q^{-m/2})$ だから $q^{-m/2}/\cosh(\tfrac m2\log q)=2/(q^m+1)$。よって第 $m$ 項 $2\log q\cdot q^{-m/2}F(m\log q)$ は $4\log q\,e^{-(m\log q/c)^2}/(q^m+1)$ に等しい。**正しい。** コードの `w=4*L/(qm+1)`、`(-(m*m)*L2*ic[c]).exp()` もこれと一致する。
- **打ち切り。** 各項は正なので、`m==8 or qm>10**16` で止めても下界である。項 $m$ を加えたあとで判定するので、加えた項に漏れもない。最小の完全分解素数 $487$ でも、落とすのは $m=7,8$ の項（$<10^{-17}$）だけである。
- **整数への変換と合計。** `range_sums` は arb の和（精度 128 ビット）の下端を `_frac` で正確な有理数にし、$\lfloor\cdot10^{30}\rfloor$ を返す。整数の和は正確で、各床関数は下端以下なので、合計 $\le$ 真の和 $\times10^{30}$。`arb(total)/SCALE` は精度 200 ビットで計算され（$10^{30}$ による割り算は丸めを含む球になる）、その `.lower()` を使う。**厳密な下界になっている。**
- **要約の検査**（§3.1 の表）。bool（`range` に `false`、`split_primes` に `true`）、float の範囲、長さ 3 の範囲、$A=B$、負の `rows`、JSON でないファイル、JSON の配列、重複キー、同じ範囲の重複、架空の区間、個数 +1、和 ×100、和が `NaN`／文字列、はすべて失敗（または例外）になった。抜けとして見つかったのは次のとおりで、どれも偽の結論にはつながらない。
  - 体は要約どうしの一致しか見ていない。`[3.0,4]` は `[3,4]` と等しいとみなされて通る（無害）。全要約の体を書き換えると別の体の計算になる（今回の要約では個数が合わず失敗した）。
  - 巨大な整数（終点 $10^{400}$）は型の検査を通り、`split_primes` が終わらない（DoS のみ。合格は出ない）。
  - `sums` の値が数でない文字列だと、`need` ではなく `arb()` の例外で止まる（終了コード 1。安全な向き）。
  - 終点が予定の値か、`run`／`verify_sha256` が記録されているかは見ていない。既存の 83 個の要約は修正前のワークフローで作られたもので、これらの項目を持たない。
- **要約に残る信頼。** 結論に効く要約の値は、`range`（どの区間の素数を数えるか）、`field`、および「その区間が `verify.py` に合格した」という事実だけである。`split_primes` は照合に、`sums` は照合だけに、`rows` は表示にだけ使われる。`README` の説明は正しい。ただし `field` も結論に効く値に入れて書くべきである。なお、基礎部分 $[0,10^7)$ の要約は `--complete` で作られたもので、手元で再実行して同じ値になった（§2.3）。
- **multiprocessing。** `imap_unordered` の結果は、ワーカーが入力の $(A,B)$ をそのまま返し、それをキーにした辞書に入れている。同じ $(A,B)$ が 2 個ある要約は先に連結の検査で落ちるので、キーの衝突はない。**対応づけは正しい。** ワーカーは fork のあとで `K(p,n)` を作り、`ctx.prec` もワーカー内で設定している。
- **終了コード（軽微）。** `VERDICT: bound not below ... yet` でも終了コード 0 で終わる。

### 1.3 `subfields.py`

- **論理は正しい。** $j=1,\dots,n-1$ について、Miller の補題を $S=\emptyset$ で使うと $h(k_{p^j})<$ 上界になる。$\ell\ne p$ が $h$ を割るなら、階数定理（Masley Cor. 2.15、van der Linden Thm 6, 8）により $\ell^f\mid h$ となる（$f$ は $\ell \bmod p^j$ の位数）。この定理は中間体 $k_{p^i}$（$i<j$）の類数が $\ell$ で割れないことを仮定するが、それは前の段で示されている。$p\mid h$ は岩澤の定理で否定される。`hmin` の探索は、$\ell\ge$ 暫定の最小値になったところで止まる。$\ell^f\ge\ell$ なので、これで正しい。
- プログラムは $j$ の順序を強制していない（全段が合格したときだけ PASS を出す）。帰納法としては、全段が合格すれば足りるので問題ない。
- $k_{p^n}$ の真部分体は $k_{p^j}$（$j<n$）と $\mathbb{Q}$ だけである（巡回拡大）。`verify.py` の事実 (3) の仮定はこれで満たされる。

### 1.4 `audit/relations_recheck.py`

- **同じ主張を検査している。** 行ごとに、係数の形式（添字が $0\le i<N$、重複なし、係数 $\ne0$）、素因数が異なること、各 $r$ が素数で $e\ge1$、ノルム、$q$ が完全分解する素数で指数 1、$A\le q<B$、他の完全分解素因数 $<q$ を検査する。ファイル全体では、$[A,B)$ の完全分解素数と、正しい行の $q$ の集合が一致すること（欠落も余分もない）を検査する。格子ブロックの行は読めずに不合格になる（`verify.py` の範囲の形式と同じく禁止）。
- **終結式の向きと指数。** $\Phi_M$ はモニックなので $\mathrm{Res}(\Phi_M,A)=\pm\prod_{\zeta}A(\zeta)=\pm N_{\mathbb{Q}(\zeta)/\mathbb{Q}}(\alpha)=\pm N_{K/\mathbb{Q}}(\alpha)^{[\mathbb{Q}(\zeta):K]}$ で、$[\mathbb{Q}(\zeta_M):K]=\varphi(M)/N=p-1$。$x\mapsto x^{p-1}$ は非負の整数の上で単射で、$\prod r^e\ge1$ だから、$|\mathrm{Res}|=(\prod r^e)^{p-1}$ は $|N(\alpha)|=\prod r^e$ と同値である。PARI の `Pol(v[::-1])` は先頭が最高次の係数なので、$A(x)=\sum v_e x^e$ になる。**正しい。**
- **生成元と剰余類。** `relations_recheck.py` の $g$（位数が $\varphi(M)$ の最小の元）と `verify.py` の $g$ は、$p=3$ で 2、$p=7$ で 3 となって一致し、剰余類の表 `E` も `F.exps` と完全に一致した。そもそも健全性は添字の取り方に依存しない。どの $g$ でも $\alpha$ は $O_K$ の元で、主張は「そういう元が存在する」ことだけだからである。
- **素数判定。** PARI の `isprime` は、BPSW 判定が保証される範囲（$10^{13}$ ないし $10^{15}$ 未満）ではそれで決定し、それを超えると Pocklington–Lehmer 型か APRCL で証明する（[PARI 開発者の 2009 年の説明](https://pari.math.u-bordeaux.fr/archives/pari-dev-0908/msg00017.html)。手元の PARI は 2.19.1 で、この版の文書は確認していない）。ここで現れる数は $3\cdot10^{10}$ 程度で、どちらの場合も証明になる。
- **cypari2 の初期化。** 親プロセスでは初期化せず、ワーカーの中で初めて `cypari2.Pari()` を呼ぶ。指示どおりである。
- **差異（無害）。** (i) ファイル先頭の `# field` 行を読まず、引数の $p,n$ を使う（ヘッダを `# field 7 2` にしても合格する。中身は $k_{3^4}$ の行として正しい）。(ii) `str.splitlines()` は改ページ文字などでも行を分ける。`verify.py` は分けない。2 本の正しい行を改ページ文字でつないだファイルは、`relations_recheck.py` では合格し、`verify.py` では例外で落ちる。どちらの検証器も単独で健全なので、結論には影響しない。(iii) 先頭に空白のある `#` 行は行として読まれ、不合格になる（安全な向き）。

### 1.5 ワークフロー（読むだけ。起動していない）

- **`relations.yml` の `gen`。** `shell: bash` と `set -euo pipefail` があり、`verify.py | tee` の失敗も `grep` の不一致も、その段で止まる。`read A B < <(...)` は、Python が失敗すると `read` が EOF で失敗するので止まる。要約のアップロードは成功時だけ（`if: always()` を外した）で、`if-no-files-found: error`。失敗時の行は `failed-range-*` という名前でアップロードされ、`collect` の `pattern: range-*` には一致しない（パターンは名前全体に対する glob）。**ここは正しく直っている。**
- **`relations.yml` の `collect`（未修正の点）。** `shell` の指定がないので `bash -e` で、`pipefail` はない。`cp ... || true` と `(cd relations && python3 aggregate.py ...) || true` があり、**集計が失敗しても run は成功と表示される**。ただし `aggregate.py` の出力は commit されず、README の数値は手元での実行から来ているので、偽の結論が入る経路ではない。
- **`recheck.yml`。**
  - `check` は `set -euo pipefail`、チェックサムの比較、2 つの検証器の `PASS` 行の存在（`test -s`）を確かめる。`test "$(cut …)" = "$(sha256sum …)"` は、両辺が空になる状況（ファイルがない）が直前の `gunzip` の成功で排除されているので健全である。
  - **release の取得。** タグを API の順に試し、**最初にダウンロードできたもので止まる**。同じ名前のファイルが複数の release にあり（同じ区間を再実行した場合）、最初のものが現在の sha256 と違えば、チェックサムで失敗する。偽の合格にはならないが、誤った失敗になる。チェックサムが一致するまで次のタグを試すべきである。
  - **`collect`（中程度）。** `if: always()` で、届いた結果だけを commit し、件数が足りなくても失敗しない（件数は表示されるだけ）。commit される JSON を検査するプログラムもない（両方が PASS か、`rows_sha256` が commit 済みの sha256 と一致するか、`verify_sha256` がレビューした版と一致するか、82 個そろっているか）。
  - `plan` の `echo "list=$(ls … | … )"` は、コマンド置換の中の失敗では止まらない。ただ、空のリストは matrix の段で失敗になる。
  - 基礎証明書 $[0,10^7)$ は再検証の対象に入っていない（sha256 ファイルがない）。基礎証明書はリポジトリにあり、前回の監査で全行を独立に検証済み、今回も `verify.py --complete` を再実行して合格した。
  - cypari2（したがって PARI の版）は固定されていない。

## 2. 独立に再計算した数値と照合

### 2.1 `aggregate.py` の `range_sums` と自作の実装（`indep_sums.py`、`cmp_sums.py`）

自作の実装は `verify.py` も `aggregate.py` も import しない。完全分解素数は、前回の監査の Eratosthenes の篩の結果（`split_primes.npy`）を使い、さらに $243k\pm1$ の上で決定的な Miller–Rabin（底 2〜37、$3.3\cdot10^{24}$ 未満で正しい）で一覧を作り直して一致を確かめた。和は元の形 $2\log q\cdot q^{-m/2}e^{-(m\log q/c)^2}/\cosh(\tfrac m2\log q)$ で、$m=1..8$ を打ち切りなしで、arb（精度 160 ビット）で計算した。比べたのは $\lfloor\cdot10^{30}\rfloor$ の整数である。

| 区間 | 個数（自作／`range_sums`／要約） | 5 つの $c$ で「`range_sums` の下界 ≤ 自作の上端」 | 自作の下端 − `range_sums`（$\times10^{-30}$、最大） |
|---|---|---|---|
| $[0,10^7)$ | 8,160／8,160／8,160 | すべて成り立つ | $6.3\times10^{12}$（$=6.3\times10^{-18}$） |
| $[10^7,1.2\cdot10^7)$ | 1,528／1,528／1,528 | すべて成り立つ | $6.3\times10^{6}$ |
| $[1.2\cdot10^7,1.4\cdot10^7)$ | 1,532／1,532／1,532 | すべて成り立つ | $3.2\times10^{6}$ |
| $[1.4\cdot10^7,388825000)$ | 244,804（3 つとも） | すべて成り立つ | $2.1\times10^{12}$ |
| $[15007000000,15381825000)$ | 197,612（3 つとも） | すべて成り立つ | $4.7\times10^{6}$ |
| $[29625175000,3\cdot10^{10})$ | 191,725（3 つとも） | すべて成り立つ | $6.1\times10^{5}$ |

ずれはすべて非負で $10^{-17}$ 未満である（`range_sums` の $m$ の打ち切りと、精度 128 ビットの球の半径の分）。たとえば $c=100$、$[0,10^7)$ では、自作が $[0.493955675793424875770022144582,\ \ldots583]$、`range_sums` が $0.493955675793424870387945666318$ だった。

### 2.2 `aggregate.py` の全体の実行

`python3 -O aggregate.py k3_4 2`（2 分 14 秒）の出力は、`relations/k3_4/aggregate_output.txt` の注釈行を除いた全行と**完全に一致**した：83 区間、$[0,3\cdot10^{10})$、完全分解素数 16,049,777、$c=100$ で $h<99.8767$、`VERDICT: PASS (h < 99.8767 < 163)`。

### 2.3 `verify.py` の再実行

| 実行 | 結果 | 前回・commit 済みの値との比較 |
|---|---|---|
| `verify.py cert_k7_2.txt` | 29,834 行、格子ブロック（19）、`VERDICT: PASS (h < 152.0639 < 197)` | 一致（$c$ ごとの 5 行も一致） |
| `verify.py k3_4/base_cert.txt --complete 10000000` | 8,699 行、split_primes 8,160、BASE PASS | `summary_0_10000000.json` と和が一致（新しい形式は文字列で、末尾の 0 が付くだけ） |
| `--range 10000000 12000000` | RANGE PASS、1,528 行 | 要約と一致 |
| `--range 15007000000 15381825000` | RANGE PASS、197,612 行（3 分 24 秒） | 要約と一致 |

### 2.4 `subfields.py` の数値（`sub_indep.py`。`verify.py` は使わない）

積分は前回の監査の Simpson 法（`/home/claude/rel_audit/code/integral.py`：幅 $\min(0.004,x/40)$、剰余項は arb_series による 4 階導関数の上界）、判別式は導手判別式公式と PARI の `nfdisc(polsubcyclo(·))`、閾値は sympy の `n_order` で計算した。

| 体 | $c$ | 判別式（nfdisc と一致） | 閾値（自作／主張） | $h<$（自作） | 主張 |
|---|---|---|---|---|---|
| $k_3$ | 2 | $3^4$ | 4／4 | 1.345417757 | 1.3455 |
| $k_{3^2}$ | 3 | $3^{22}$ | 19／19 | 1.331330435 | 1.3314 |
| $k_{3^3}$ | 8 | $3^{94}$ | 109／109 | 5.960014879 | 5.9601 |
| $k_7$ | 5 | $7^{12}$ | 8／8 | 4.802110966 | 4.8022 |
| $k_3$（参考） | 3 | | | 1.681453171 | 前回 1.682 |

主張はすべて自作の値以上（安全な向き）で、どれも閾値より小さい。

## 3. 試した改竄と結果

### 3.1 `aggregate.py`（`mk_aggtests.py`。基礎 + 2 区間の小さな要約の組に改竄を入れた）

| 改竄 | 結果 |
|---|---|
| 改竄なし（基礎 + 2 区間） | `VERDICT: bound not below 163 yet`（終了コード 0） |
| 基礎の範囲を `[false, 10000000]` | FAILED: the range must be [A,B) |
| 2 番目の範囲の始点を `false` | FAILED（同上） |
| `split_primes` を `true` | FAILED: number of rows |
| `field` を `[true,4]` | FAILED: field |
| `field` を `[3.0,4]`（1 個の要約だけ） | 通る（`[3.0,4]==[3,4]`。無害） |
| 全要約の `field` を `[3,3]` | FAILED: 24627 split primes found, the summary states 8160 |
| 範囲を float、長さ 3、$A=B$ | FAILED |
| `rows` を負 | FAILED |
| 個数 +1（`rows` も +1） | FAILED: 1532 split primes found, the summary states 1533 |
| 和 ×100 | FAILED: the sum in the summary disagrees |
| 和を `NaN` | FAILED: disagrees |
| 和を `"abc"` | `ValueError`（終了コード 1） |
| 和を `"0.009"`（誤り） | FAILED: disagrees |
| `sums` を削除 | 通る（和は計算し直した値を使うので無害） |
| 架空の区間 $[1.4\cdot10^7,1.5\cdot10^7)$、個数 0 | FAILED: 744 split primes found |
| 同じ範囲の要約を 2 個 | FAILED: not contiguous |
| 重複キー（`range` を 2 回） | 後の値が使われ、個数が合わず FAILED |
| JSON でない、JSON の配列 | 例外／FAILED |
| 終点 $10^{400}$ | 60 秒で打ち切り（終わらない。合格は出ない） |

### 3.2 行の改竄（範囲 $[10^7,1.2\cdot10^7)$。`verify.py --range` と `relations_recheck.py`）

前回の 15 種類（`/home/claude/rel_audit/adv/t*.txt`）と、今回作った 19 種類（`/home/claude/rel_audit2/adv2/`。うち 3 種類は無害な変形）。

| 改竄 | `verify.py` | `relations_recheck.py` |
|---|---|---|
| a01 係数の符号を 1 つ反転 | FAILED: norm | BAD: norm |
| a02 添字 $=N=81$ | FAILED: coefficients | BAD: coefficients |
| a03 負の添字 $-58$（Python では 23 と同じ要素を指す） | FAILED: coefficients | BAD: coefficients |
| a04 素因数 $-3$ | FAILED: not prime | BAD: not prime |
| a05 素因数 $1^1$ を追加 | FAILED: not prime | BAD: not prime |
| a06 欄を 1 つ追加 | 例外で停止 | BAD: cannot read |
| a07 $q=B$（範囲の外、合成数） | FAILED: norm | BAD（3 項目） |
| a08 $q^1$ を 2 回書く | FAILED: norm（重複の検査と同じ `need`） | BAD: repeated prime |
| a09 別の行の $q$ を付け替え | FAILED: q must be a split prime… | BAD（2 項目） |
| a10 正しい 2 行を改ページ文字でつなぐ | 例外で停止 | PASS（2 行とも正しい。無害） |
| a11 行を 1 本削除 | FAILED: 1 split primes without a row | FAILED: 1 without a row |
| a12 ヘッダを `# field 7 2` | FAILED: norm | PASS（ヘッダを読まない。行は正しい） |
| a13 格子ブロックを追加 | FAILED: not allowed | BAD: cannot read |
| a14 $c$ を 1 ずらす | FAILED: norm | BAD: norm |
| a15 $c=10^{40}$ | FAILED: not enough auxiliary primes | BAD: norm |
| a16 素因数の欄を空に | FAILED: norm | BAD（2 項目） |
| b01 素因数の順序を入れ替え（無害） | PASS | PASS |
| b02 改行を CRLF に（無害） | PASS | PASS |
| b03 `+6:+1`、`6_803`（無害） | PASS | PASS |
| t1〜t4、t6〜t14（前回と同じ） | すべて FAILED（前回の監査 §4.1） | すべて BAD／FAILED |
| t5 行の重複 | PASS（前回と同じ。無害） | PASS（rows=1529。無害） |
| t15 ヘッダを `# field 3 3` | FAILED: norm | PASS（ヘッダを読まない） |

通してはいけない入力で通ったものはない。`relations_recheck.py` が通した a10、a12、t15、t5 は、行そのものが $k_{3^4}$ の正しい行である。

### 3.3 全行の再検証（追加分）

`relations_recheck.py` で `rows_22128675000_22503500000.txt`（前回は 1,000 行の標本だけ）の**全 194,209 行**を検査し、RECHECK PASS だった（2 プロセスで 6 分 23 秒）。sha256 は `relations/k3_4/rows_22128675000_22503500000.sha256` と一致した。小さな区間 $[10^7,1.2\cdot10^7)$ も RECHECK PASS。

## 4. 残っている信頼

1. **79 個の範囲ファイルの行。** 全行を独立に検証したのは 3 個（$[10^7,1.2\cdot10^7)$、$[15007000000,\ldots)$、$[22128675000,\ldots)$）で、標本だけのものが 4 個ある。残り 75 個と標本の残りは、Actions 上での `verify.py` の合格（ログと要約の一致は前回の監査で確認）を信頼している。`recheck.yml` を実行して 82 個すべてが合格すれば、この信頼は「`relations_recheck.py` と `verify.py` の両方が誤る」ことへの信頼にまで減る。ただし `relations_recheck.py` は作業側が書いたもので、コードは `verify.py` から独立しているが、作者は独立していない（今回の監査で読んだ）。
2. **既存の 83 個の要約** は修正前のワークフローで作られ、`run`、`commit`、`verify_sha256` を持たない。由来の確認は、前回の監査で行った Actions のログとの照合に頼っている。
3. **要約のうち結論に効く値**：`field`、`range`、「合格した」という事実（§1.2）。
4. **文献とソフトウェア**：前回の §5 と同じ（Miller の補題、階数定理、岩澤の定理、FLINT/Arb、PARI）。PARI の `isprime` の証明の方式は、2009 年の開発者の説明で確かめた。手元の版の文書は確かめていない。
5. **確かめていないこと**：ワークフローの実際の挙動（起動していない。§1.5 は YAML と GitHub Actions の一般的な仕様に基づく推論）。`g*.txt`（$k_{7^2}$ の格子ブロックの改竄）を `relations_recheck.py` にかけること（このスクリプトは格子ブロックを扱わないので対象外とした）。

## 5. 直すべき点（重大度の順）

1. **（中）`recheck.yml` の結果の検査。** `collect` で、件数が `rows_*.sha256` の個数に足りなければ失敗させる。あわせて、`recheck/*.json` を読んで、全区間について両方が PASS、`recheck.rows_sha256` が commit 済みの sha256 と一致、`verify.field==[p,n]`、`verify_sha256`／`recheck_sha256` が監査した版と一致、区間が `aggregate.py` の区間と同じ、を確かめる小さなプログラムを用意する。
2. **（低〜中）`relations.yml` の `collect`。** `aggregate.py` の `|| true` を外す（または結果を明示して失敗させる）。`shell: bash` と `set -euo pipefail` を付ける。`aggregate.py` は、PASS でないとき 0 以外の終了コードで終わる。
3. **（低）`aggregate.py` で体と終点を固定する。** 引数（またはディレクトリ名）から $(p,n)$ を取り、全要約の `field` と型まで一致させる（`[3.0,4]` を拒否）。終点も予定の値と比べる。
4. **（低）`recheck.yml` の release の取得。** チェックサムが一致するまで次のタグを試す。
5. **（低）`_dec` の $k=0$。** `need(k>=1)` を入れる。
6. **（低）`relations_recheck.py` の読み方を `verify.py` に合わせる。** 行の区切りを `\n` だけにし、ヘッダの `# field p n` を引数と照合する。
7. **（低）版の固定。** `recheck.yml` と `relations.yml` で cypari2 の版も固定する。基礎証明書を再検証の対象に入れるか、対象外であることと、前回の監査で全行を検証済みであることを README に書く。

作業用のファイル：`/home/claude/rel_audit2/`（`t_updown.py`、`indep_sums.py`、`cmp_sums.py`、`sub_indep.py`、`mk_aggtests.py`、`mk_adv2.py`、`aggtests/`、`adv2/`、各 `*.log`）。

## 対応の記録（監査者ではなく、作業側が 2026-11-24 に追記）

| 指摘 | 対応 |
|---|---|
| 1. `recheck.yml` の `collect` が欠落を見逃す | `relations/check_recheck.py` を追加。チェックサムのある範囲ファイルすべてに結果があり、両方の検査が合格し、チェックサムと完全分解素数の個数が一致することを確かめる。ワークフローの最後でも実行する |
| 2. `relations.yml` の `|| true`、終了コード | 集計を独立した段階にして `|| true` を外した。`aggregate.py` は、上界が閾値を下回らないとき終了コード 1 で終わる。`collect` に `pipefail` を付けた |
| 3. 体と終点が固定されていない | `aggregate.py DIR procs p n END` と書くと、体と終点が一致しないとき失敗する |
| 4. 最初の release で止まる | チェックサムが一致するファイルが見つかるまで、全 release を順に試す |
| 5. `_dec` の k=0 | `need(k>=1)` を追加 |
| 6. 行の読み方 | `relations_recheck.py` を `verify.py` と同じ読み方（改行だけで区切る）にした |
| 7. 版の固定、基礎証明書 | `recheck.yml` で cypari2 を 2.2.4（実行時に入っていた版）に固定。基礎証明書は `recheck.yml` の対象外で、1 回目の監査で全 8,699 行が別実装で検証されている（`audit_relations.md` §3.2） |

**全行の再検証の結果（Actions の実行 45380165270、commit 0fdb1cc）。** 82 個の範囲ファイルすべてで、`relations_recheck.py`（PARI）と `verify.py`（修正後）の両方が合格した。`python3 relations/check_recheck.py relations/k3_4` の出力は「82 range files, 16041617 rows, 16041617 split primes; both checks passed for every file」。基礎部分の 8,160 個と合わせて 16,049,777 個になる。

- この実行の後に変えたのは、`verify.py` の `need(k>=1)` の 1 行と、`relations_recheck.py` の行の読み方の 1 行だけ。
- 修正後のワークフロー 2 本（`relations.yml`、`recheck.yml`）は、修正後には実行していない。
