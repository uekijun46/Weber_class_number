# 文献照合（2026-11-22）

> 記号について：この文書は以前の記号で書かれている。$B_{p,n}$ は現在の $k_{p^n}$、$\mathbb{Q}(m)$ は現在の $k_m$、「層」「合成層」は「layer」「合成数次の layer」に当たる（README の冒頭を参照）。

到達できた文献の範囲での結果。「読んで確認」と「推論」を区別して書く。

## 合成層 $\mathbb{Q}(m)$ の類数 1

| $m$（導手） | 状況 | 根拠 |
|---|---|---|
| 6 (72) | 既知・無条件 | Masley 1978 Thm 4.5、van der Linden 1982 Thm 3、Miller 2014 Thm 1.1 から降下 |
| 12 (144) | 既知・無条件 | van der Linden Thm 3、Miller Thm 1.1 から降下 |
| 10 (200) | 既知・無条件 | Miller 2014 §7 / Table 1（$h^+_{200}\le1$）から降下。van der Linden は GRH 下 |
| 18 (216) | 既知・無条件 | Miller §7 / Table 1 から降下 |
| 24 (288) | 既知・無条件 | Miller Prop. 7.11 から降下 |
| 15 (225), 14 (392), 20 (400), 21 (441), 22 (968), 26 (1352), 28 (784), 36 (432) | GRH 下のみ | LMFDB に「(GRH)」付きの記載。Miller 2014 の範囲（$\varphi(f)\le116$）の外 |

降下：$\mathbb{Q}(\zeta_f)^+/\mathbb{Q}(m)$ に不分岐な部分拡大がないので、$h(\mathbb{Q}(\zeta_f)^+)=1$ なら $h(\mathbb{Q}(m))=1$。13 体すべてで計算により確認した。

## 類数が 1 でない既知の例

Gras（arXiv:2009.05278）の脚注と Fukuda のスライド（2018）による一覧。

- Horie 2001：$31\mid h(\mathbb{Q}(2\cdot31))$、$73\mid h(\mathbb{Q}(3\cdot73))$。
- Fukuda–Komatsu 2011：上の 2 例と $1546463\mid h(\mathbb{Q}(2\cdot1546463))$。
- Fukuda–Komatsu–Morisawa 2011：$18433\mid h(\mathbb{Q}(2^8\cdot18433))$、$114689\mid h(\mathbb{Q}(2^{10}\cdot114689))$、$487\mid h(\mathbb{Q}(3^4\cdot487))$、$238627\mid h(\mathbb{Q}(3^4\cdot238627))$、$2251\mid h(\mathbb{Q}(5^2\cdot2251))$。
- Aoki–Fukuda：$107\mid h(\mathbb{Q}(2\cdot53))$。$p\nmid m$ 型の唯一の既知例。

既知の最小次数は 62、次が 106、219。$m<62$ の例は見つかっていない。

## 背景

- Coates の問い：J. Coates, "The enigmatic Tate–Shafarevich group", Fifth ICCM, AMS/IP Stud. Adv. Math. 51 (2012), §3。$\mathbb{Q}(N)$ の類群が自明か、少なくとも極めて稀か。
- Fukuda–Komatsu, "Weber's class number problem in the cyclotomic $\mathbb{Z}_2$-extension of $\mathbb{Q}$, III", Int. J. Number Theory 7 (2011)：$\ell<10^9$ は $h(B_{2,n})$ を割らない（吉崎サーベイ arXiv:2211.15201 と Miller の引用で確認）。
- Buhler–Pomerance–Robertson, Fields Inst. Commun. 41 (2004)：有限個の $(p,n)$ を除き $h(\mathbb{Q}(\zeta_{p^{n+1}})^+)=h(\mathbb{Q}(\zeta_{p^n})^+)$ という発見的予想（Miller arXiv:1410.2921 の記述による）。

## 主な文献

- J. M. Masley, Compositio Math. 37 (1978) 297–319。
- F. J. van der Linden, Math. Comp. 39 (1982)。
- J. C. Miller, "Class numbers of real cyclotomic fields of composite conductor", LMS J. Comput. Math. 17 (2014) 404–417, doi:10.1112/S1461157014000382。
- J. C. Miller, "Class numbers in cyclotomic $\mathbb{Z}_p$-extensions", arXiv:1410.2921。
- K. Horie, "A note on the $\mathbb{Z}_p\times\mathbb{Z}_q$-extension over $\mathbb{Q}$", Proc. Japan Acad. 77A (2001) 84–86, doi:10.3792/pjaa.77.84。
- G. Gras, "Weber's class number problem and p-rationality in the cyclotomic $\hat{\mathbb{Z}}$-extension of $\mathbb{Q}$", arXiv:2009.05278。
- M. Aoki, T. Fukuda, LNCS 4076 (ANTS VII, 2006) 56–71, doi:10.1007/11792086_5。
- T. Fukuda, Bordeaux seminar slides (2018), https://lfant.math.u-bordeaux.fr/seminar/slides/2018-03-06T10:00--Takashi_Fukuda.pdf

## 未確認

- Horie 2001、Coates 2012、Buhler–Pomerance–Robertson、Fukuda–Komatsu III の原論文は読めていない（bot チェックなど）。
- Masley の導手 72 と van der Linden の範囲（$\varphi(f)\le72$）は、スキャンが不鮮明なため一部推論。
- LMFDB で「(GRH)」の印がないことを無条件と解釈している（$m=6,10,12$）。
- Horie、Morisawa の $\mathbb{Z}_p\times\mathbb{Z}_q$ 塔の論文は未読。8 体について無条件証明が他にないとは言い切れない。

## 追記（2026-11-24）：$k_{7^2}$ のノートのための確認

VPN なしで到達できる範囲（arXiv、ar5iv、出版社の公開ページ、J-STAGE、numdam）での確認。記号は現在のもの（$k_{p^n}$）。

### 既知の結果の出典（Miller, arXiv:1410.2921 の Table 1、Theorem 1.8、Table 2 による）

| 体 | 仮定 | 出典 |
|---|---|---|
| $k_{2^3}$ | なし | Cohn (1960) |
| $k_{2^4}$、$k_{3^2}$、$k_5$、$k_7$ | なし | Bauer (1969) |
| $k_{2^5}$ | なし | van der Linden (1982) |
| $k_{3^3}$ | なし | Miller の表では Masley (1978)。ただし下の 2026-11-24 夜の追記を参照 |
| $k_{2^6}$ | なし | Miller（Acta Arith. 2014） |
| $k_{5^2}$、$k_{11}$ | なし | Miller（LMS J. Comput. Math. 2014） |
| $k_{13}$、$k_{17}$、$k_{19}$ | なし | Miller（J. Number Theory 2015, Theorem 1.8） |
| $k_{2^7}$ | GRH | Miller |
| $k_{3^4}$ | GRH | van der Linden (1982) |
| $k_{23}$ | GRH | Coates–Liang–Mihăilescu (2012) |
| $k_{29}$、$k_{31}$ | GRH | Miller（J. Number Theory 2015, Table 2） |

- $k_{2^3}$ については、Weber が示したとする文献（Morisawa–Okazaki の序文、吉崎サーベイ）もある。
- $k_{23}$ の無条件証明を主張する GitHub リポジトリがある（`johnchakkour/b23-class-number`、上界 $h<37.717103$、未査読）。

### $k_{7^2}$ の先行結果

- 見つからなかった。Miller の表にあるのは $k_7$ だけ。
- ある論文の引用文に「$(p,n)=(7,2)$ で $\mathbb{Q}(\zeta_{p^n})^+$ の類数が 1」とあるが、これは導手 $7^2=49$ の実円分体のことで、含まれるのは $k_7$。$k_{7^2}$ の導手は $7^3=343$。
- Horie–Horie、Fukuda–Komatsu、Morisawa–Okazaki の結果は、素数 $\ell$ を固定したときの「$\ell$ が $h(k_{p^n})$ を割らない」型で、類数そのものの決定ではない。
- 検索した範囲で、2015 年以降に $k_{p^n}$ の類数 1 を新たに証明した査読論文は見つからなかった。

### GRH を仮定した $k_{37}$、$k_{41}$、$k_{43}$、$k_{47}$

- 検索した範囲では、結果は見つからなかった。Miller の GRH 下の結果は $k_{29}$、$k_{31}$ まで。

### 書誌情報

確認できたもの：
- Miller, Acta Arith. 164 (2014), 381–397（出版社のページ。AMS の引用では 381–398）。
- Miller, Math. Comp. 84 (2015), 2459–2469。
- Miller, LMS J. Comput. Math. 17 (2014), Special issue A, 404–417。
- Bauer, J. Number Theory 1 (1969), 161–162。
- Cohn, Numer. Math. 2 (1960), 347–362。
- van der Linden, Math. Comp. 39 (1982), 693–707。
- Fukuda–Komatsu III, Int. J. Number Theory 7 (2011), 1627–1635。
- Horie–Horie, Tohoku Math. J. 61 (2009), 551–570。Morisawa–Okazaki, Tohoku Math. J. 65 (2013), 253–272。

確認できなかったもの：
- Miller, J. Number Theory 150 (2015) のページ番号（出版社のページに到達できず）。
- Johansson の ICMS 2018 の論文のページ番号（arXiv:1802.07942 は確認）。

### まだ確認できていないこと

- MathSciNet と zbMATH での被引用文献の確認（VPN が必要）。
- Hakkarainen の表の原本の目視。
- Masley、Bauer、Cohn の原論文の本文。

## 追記（2026-11-24 夕方）：MathSciNet での確認

VPN 経由で MathSciNet を使った確認。記号は現在のもの。

- **重要な言い換え**：$k_{3^4}=\mathbb{Q}(\zeta_{243})^+$ である（$(\mathbb{Z}/243)^\times/\{\pm1\}$ が位数 81 の巡回群なので）。したがって $h(k_{3^4})=1$ は「導手 243 の実円分体の類数が 1」と同じ主張。同様に $k_{3^n}=\mathbb{Q}(\zeta_{3^{n+1}})^+$。一方 $k_{7^2}$ は $\mathbb{Q}(\zeta_{343})^+$（147 次）の真部分体。
- Miller, J. Number Theory 150 (2015), 47–73（MR3304606）。ページ番号を確認。レビューにある既知の結果の一覧は、上の表と一致。$k_{3^4}$（$\mathbb{B}_{3,4}$）は「GRH 下で van der Linden」とされている。
- 被引用文献の確認：
  - MR3304606（Miller, JNT 2015）を引用する 5 件
  - MR3244941（Miller, Acta Arith. 2014）を引用する 21 件
  - MR3240817（Miller, LMS JCM 2014）または MR3356035（Miller, Math. Comp. 2015）を引用する 25 件
  - MR0669662（van der Linden 1982）を引用する 2015 年以降の 15 件
  
  題名と、関係しそうなもののレビューを見た範囲で、$k_{7^2}$ または $k_{3^4}$（導手 343、243）の類数 1 を無条件で証明したものは見つからなかった。
- キーワード検索（"class number" cyclotomic 243、同 343、"Z_3-extension"）でも該当なし。
- 関連する文献（文章化のときに引用の候補）：
  - Fukuda–Komatsu–Morisawa, "Weber's class number one problem", Iwasawa theory 2012, 221–226, Contrib. Math. Comput. Sci. 7, Springer, 2014（MR3586814）：すべての $n$ で、$10^9$ 未満の素数は $h(k_{3^n})$ を割らない（Miller のレビューによる）。
  - Kogoshi–Morisawa, Exp. Math. 29 (2020), 470–474：$k_{29}$、$k_{31}$ の類数について。
  - Bernstein, "Fast norm computation in smooth-degree Abelian number fields", Res. Number Theory 9 (2023), Paper No. 82：次数が滑らかなアーベル体でのノルムの高速計算。
  - Gras, J. Number Theory 228 (2021), 219–252（arXiv:2009.05278 の出版版。題名は "Tate–Shafarevich groups in the cyclotomic $\hat{\mathbb{Z}}$-extension and Weber's class number problem"）。
- 確認できなかったもの：Johansson の ICMS 2018 の論文は MathSciNet の検索に出ず、ページ番号は未確認のまま。Hakkarainen の表の原本の目視も未実施。

## 追記（2026-11-24 夜）：$k_{3^3}$ の典拠

- 監査（`docs/audit_relations.md` §1.3）の指摘。Miller の Table 1 は $h(k_{3^3})=1$（導手 81）の典拠を Masley (1978) としているが、不正確。
- MathSciNet のレビュー（MR0511747）で確認：Masley が $h=1$ を示したのは、素数冪導手 **70 未満**の実巡回体。100 まで延ばせることは予想として書かれている。
- van der Linden (1982) のレビュー（MR0669662）：Masley の結果と方法を延ばして、導手 200 以下のほとんどの実アーベル体と、導手 256、243 の場合を扱う。一部は GRH を仮定。導手 81 が無条件の側に入っているかは、本文を読めていないので未確認（AMS のページは自動閲覧の確認画面で止まる）。
- 対処：文献に頼らない証明を用意した。`python3 relations/subfields.py 3 4` は、素数を使わない Miller の補題で $h(k_3)<1.3455$、$h(k_{3^2})<1.3314$、$h(k_{3^3})<5.9601$ を出す。階数定理の閾値は 4、19、109 で、岩澤の定理と合わせて 3 個とも類数 1。$k_7$ も同じ方法で $h<4.8022<8$（`subfields.py 7 2`）。

## 追記（2026-11-24 夜、その 2）：原論文の PDF で確認したこと

利用者から van der Linden (1982) と Miller (J. Number Theory 2015) の PDF を受け取り、本文を読んだ。

**van der Linden, Math. Comp. 39 (1982), 693–707**
- Theorem 1（無条件）：導手 $q$ が素数冪で $\varphi(q)\le66$ なら $h(K)=1$。
  - $k_{3^3}$（導手 81、$\varphi=54$）はこれに入る。**$h(k_{3^3})=1$ の正しい典拠は van der Linden の Theorem 1。**
  - $k_3$、$k_{3^2}$、$k_7$（導手 49、$\varphi=42$）、$k_{2^5}$（導手 128、$\varphi=64$）も入る。
- Theorem 2（GRH を仮定）：$\varphi(q)\le162$ で、導手 163 の場合（$h=4$）以外は $h(K)=1$。$k_{3^4}$（導手 243、$\varphi=162$）はここ。GRH を仮定しない場合については「主張した値で割り切れ、商に現れる素数冪に下界がある」とだけ書かれている。
- Theorem 8（Rank）：巡回拡大 $L/K$（次数 $n$）と $p\nmid n$ に対し、$\#\mathrm{Cl}_p^*(L/K)$ は $p^f$ の冪（$f$ は $p \bmod n$ の位数）。証明つき。「Masley の (2.15) をより精密にしたもの」とある。Theorem 6 と合わせると、こちらで使っている一般形の階数定理になる。**階数定理の引用先として、Masley の Cor. 2.15 よりこちらが適切。**
- Corollary 12（Pushing Down）：分岐する素点が高々 1 個の $p$ 拡大 $L/K$ で、$p\mid h(L)$ なら $p\mid h(K)$。岩澤の定理と同じ内容。

**Miller, J. Number Theory 150 (2015), 47–73**
- 出版版の番号は arXiv 版と違う。上界の補題は **Lemma 1**（arXiv 版では Lemma 2.2）、$k_{13}$、$k_{17}$、$k_{19}$ の定理は **Theorem 7**（arXiv 版では Theorem 1.8）。補題の出典は Acta Arith. の論文とされている。
- Lemma 1 の文面を確認。$S$ は「完全分解して単項素イデアルになる素数の集合の部分集合」で、こちらの使い方と合う。
- Rank Theorem は素数次の形でしか述べられていない（Masley の Cor. 2.15 を引用）。
- Table 1 は $\mathbb{B}_{3,3}$ の典拠を Masley (1978) としているが、同じ論文の Section 6 の本文には「van der Linden が $\mathbb{B}_{3,3}$ の類数 1 と、GRH 下で $\mathbb{B}_{3,4}$ の類数 1 を示した」とある。表の側の誤記と考えられる。
- Table 1 は $\mathbb{B}_{3,4}$ を「1 (under GRH)、van der Linden (1982)」としている。2015 年の時点で $k_{3^4}$ の無条件の結果は知られていない。
- Theorem 5（Fukuda–Komatsu–Morisawa）：すべての $n$ で、$10^9$ 未満の素数は $h(k_{3^n})$ を割らない。こちらの $k_{3^4}$ の証明はこの結果を使っていない（閾値 163 だけを使う）。
- Hakkarainen (Math. Comp. 78 (2009), 555–573) は、導手 2000 未満のアーベル体の類数の「約数」を計算したもの、と説明されている（$k_{37}$、$k_{41}$、$k_{43}$ で 10,000 未満の奇素数が類数を割らない）。類数そのものを決める表ではない。
- 文献表では Acta Arith. の論文のページが 381–398。MathSciNet も 381–398。以前に出版社のページとして 381–397 と記録したものと食い違う。ノートでは 381–398 に合わせた。

**ノート `paper/k49.tex` の修正**：引用の番号を出版版に合わせた（Lemma 1、Theorem 7）。階数定理の引用を van der Linden の Theorems 6, 8 に変更。$h(k_7)=1$ に van der Linden の Theorem 1 を追加。本文で確認していない Math. Comp. 2015 の論文（Lemma 3.1）への引用を削除。

## 追記（2026-11-24 夜、その 3）：Miller, Acta Arith. の PDF で確認したこと

- 書誌：Acta Arith. 164 (2014), no. 4, **381–398**（PDF は 18 ページ。最後のページ 398 は文献表と著者の住所）。EuDML の「381–397」は誤り。
- 無条件の上界は **Lemma 5.1**。文面は J. Number Theory の Lemma 1 と同じで、$S$ は「完全分解して単項素イデアルになる素数の集合の部分集合」。GRH 下の版は Lemma 5.2（$F$ は Gauss 関数、上界に因子 $e^{(c/4)^2}$）。
- (4.1)：$h=\dfrac{4\mathcal{H}(F)/n}{C-\mathcal{G}(F)-\log\mathrm{rd}(K)+\frac1{hn}\sum_\rho\Phi(\rho)+\frac2{hn}\sum_{\mathfrak{P}}\sum_m\frac{\log N\mathfrak{P}}{N\mathfrak{P}^{m/2}}F(m\log N\mathfrak{P})}$。和は Hilbert 類体のすべての素イデアルにわたる。
- $F=f/\cosh(x/2)$ で $f$ が非負かつ Fourier 変換が非負なら、$\Phi$ は臨界帯全体で非負（Odlyzko の概説を引用）。
- Section 6（導手 256、$c=210$）：完全分解素数の和 $>0.7023$ に、完全分岐する 2 の項 $>0.0331$ を足して $B>0.0021$、$h<5539$。最後は Fukuda–Komatsu の結果（小さい素数の非可除性）で $h=1$。
- 「van der Linden は導手 128 で 2 の上の分岐する素イデアルの寄与を使った」と Section 5 にある。
- ノートの修正：補題の引用を Lemma 5.1 に、ページを 381–398 に。分岐する素数の項についての注意を追加。

## 追記（2026-11-24 夜、その 4）：Johansson の 2 論文を原典で確認

利用者から受け取った PDF による。
- F. Johansson, Numerical integration in arbitrary-precision ball arithmetic, in: Mathematical Software – ICMS 2018 (J. H. Davenport et al., eds.), Lecture Notes in Comput. Sci. 10931, Springer, Cham, 2018, **pp. 255–263**（章の脚注と目次で確認。次の章は 264 ページから）。積分の方法は Petras のアルゴリズム（適応的な二分と Gauss–Legendre 求積）で、誤差の上界は厳密。
- F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017), no. 8, **1281–1292**（12 ページ、DOI 10.1109/TC.2017.2690633）。
- これで、ノート 2 本の文献表に未確認の印（CHECK）は残っていない。

## 追記（2026-11-25 夜）：GRH を仮定する $k_{37}$、$k_{41}$、$k_{43}$ についての MathSciNet での確認

VPN 経由。対象は「GRH 下で $h(k_{37})=h(k_{41})=h(k_{43})=1$」と、その方法（`docs/grh.md` の §6）。

- 検索 `any:(Weber "class number" cyclotomic)`、2014〜2026 年：20 件。題名を見た範囲で、次数 37、41、43 の layer の類数を扱うものはない。
- 検索 `any:("class number" ("Z_p-extension" …) (GRH OR "Riemann hypothesis"))`、2014〜2026 年：1 件（Miller, J. Number Theory 2015）。
- 検索 `any:("class number" cyclotomic ("37,1" OR "41,1" OR "43,1" OR "Z_{37}" OR "Z_{41}" OR "Z_{43}"))`：0 件。
- 被引用文献：MR3304606（Miller, JNT 2015）を引用する 5 件、MR3244941（Miller, Acta Arith. 2014）を引用する 21 件（2026 年の 3 件を含む）。題名を見た範囲で該当なし。
- Miller, JNT 2015 のレビュー（MR3304606）による GRH 下の既知の結果：$k_{2^7}$（Miller）、$k_{3^4}$、$k_{13}$（van der Linden）、$k_{17}$、$k_{19}$、$k_{23}$（Coates の ICCM 2010 の報告、MR2908059 を典拠としている）、$k_{29}$、$k_{31}$（Miller のこの論文）。$k_{37}$ 以上はない。
- **この範囲では、GRH 下の $h(k_{37})=h(k_{41})=h(k_{43})=1$ は見つからなかった。** レビューと題名による確認で、本文の全文検索ではない。arXiv の未出版の原稿は別に確認が要る。

**方法の位置づけ（新規性を主張しすぎないための記録）。**

- 「関係式の行列の行列式は類群の位数の倍数である」は、類群計算の標準的な事実。Hafner–McCurley, "A rigorous subexponential algorithm for computation of class groups", J. Amer. Math. Soc. 2 (1989), no. 4, 837–850（MR1002631）のレビューに、関係式の疎な行列の行列式 $h$ を小さい素数を法として計算して中国剰余定理で復元し、それが類群の位数の倍数である、とある。
- 一般の数体への拡張は Buchmann, "A subexponential algorithm for the determination of class groups and regulators of algebraic number fields", Séminaire de Théorie des Nombres, Paris 1988–1989, 27–41, Progr. Math. 91, Birkhäuser Boston, 1990（MR1104698）。レビューによると、新しい要素は単数群の扱いで、関係式の格子に単数の対数の格子を加えて類群と単数基準を同時に求める。GRH と滑らかさの仮定の下で準指数時間。
- こちらの議論はこの標準的な事実の特別な場合である：倍数が 1 になれば、単数も単数基準も要らない。違いは 3 点。(1) 類群が因子基底で生成されることを Bach の上界ではなく Miller の補題（GRH）に任せる。(2) Galois 群の作用を使って行列を $\mathbb{Z}[\zeta_p]$ 上の小さい行列にする。(3) PARI を使わず、整数の厳密な計算だけの検証器で確かめる。
- Bach, "Explicit bounds for primality testing and related problems", Math. Comp. 55 (1990), no. 191, 355–380（MR1023756）：GRH 下の明示的な上界（レビューで書誌を確認。類群の生成元の上界 $12(\log|d|)^2$ はこの論文によるとされるが、本文では未確認）。$k_{37}$ では $12(\log|d|)^2\approx8.1\cdot10^5$ で、それ未満の完全分解素数をすべて $T$ に入れる必要がある。今の $T$ は $10^6$ 未満の 2,135 個のうち 2,049 個なので、Bach の上界ではなく Miller の補題を使っている。
- 関連：Biasse–Espitau–Fouque–Gélin–Kirchner, "Computing generator in cyclotomic integer rings", EUROCRYPT 2017, Part I, 60–88, LNCS 10210（MR3652099）。円分整数環での単項イデアルの生成元の計算。Miller の Acta Arith. の論文を引用している。

## 追記（2026-11-26）：arXiv などの検索（GRH 下の $k_{37}$、$k_{41}$、$k_{43}$）

- Web 検索（3 通りの語句、arXiv を含む）で、次数 37、41、43 の第 1 layer の類数を GRH 下で決めた原稿は見つからなかった。確認できたのは検索結果の範囲で、arXiv の全件を調べたわけではない。
- Miller (2015, arXiv:1410.2921) の Table 1 では、GRH 下の $k_{17}$、$k_{19}$、$k_{23}$ は Coates–Liang–Mihăilescu (2012) に帰されている。
- **arXiv:2604.15858（Ming-Xing Luo「Module Lattice Security (Part I): Unconditional Verification of Weber's Conjecture for $k\le12$」、2026-04-17 投稿、cs.CR と quant-ph）：確認した（2026-11-26）。結論：この主張には頼れない。ノートの序文の「知られている場合」の一覧は変えない。**
  - 読み方：最初は arXiv の PDF を取得する道具を通して抜き書きを読んだ。その後、利用者が用意した PDF（v2、2026-05-07 付、24 ページ）の本文を直接読み、下の各点を該当箇所で確かめた（Table 1、Lemma 3.1、Theorem 3.1、Lemma 3.3 とその証明、Lemma 3.4、Theorem 3.2、Definition 3.2、Proposition 3.2、Example 3.1、Theorem 3.3、Algorithm 1、Proposition 3.3）。全 24 ページを通読したわけではない。
  - 主張：導手 $2^k$ の実円分体の類数 $h_k^+$ は $k=9,\dots,12$ で 1（Theorem 3.3）。出発点は Miller の $h_8^+=1$。
  - 議論の骨組み：奇素数 $\ell\mid h_k^+$ なら、位数最大の偶指標 $\chi$ の固有空間が零でない。$\psi=\chi^{-1}\omega_\ell$（$\omega_\ell$ は $\ell$ を法とする Teichmüller 指標）について $\ell\mid\mathrm{Nm}(B_{1,\psi})$（Theorem 3.2、Herbrand–Ribet として引用）。Lemma 3.3 (iii) で、これを $\ell\mid\mathrm{Nm}(B_{1,\chi^{-1}})$ と同値とする。そこで $10^9$ より大きく $\ell\equiv\pm1 \bmod 2^{k-1}$ で、ある「Bernoulli 数のノルム」を割る素数の有限集合 $S_k$ を作り、1 個ずつ調べる。
  - 問題 1（決定的）：$\psi=\chi^{-1}\omega_\ell$ は $\ell$ に依存する（導手 $2^k\ell$）。したがって「有界な 1 個の整数の素因数は有限個」という $S_k$ の有限性の議論は成り立たない。これを避けるための Lemma 3.3 (iii) は、証明がなく（「Stickelberger 元と岩澤理論を使う標準的な結果」とあるだけ）、しかも書かれているままでは意味をなさない：$\chi$ は偶指標なので $B_{1,\chi^{-1}}=0$ である。
  - 問題 2：実際に素因数分解している数は、$h_k^+$ ではなく相対類数（マイナス部分）の因子である。本文の例は $k=7$ で $|N|=692092928=2^{15}\cdot21121$。$21121=359057/17=h^-(\mathbb{Q}(\zeta_{128}))/h^-(\mathbb{Q}(\zeta_{64}))$ である（$17\cdot21121=359057$ はこちらで計算。$h^-$ の値 17、359057 は記憶による既知の値で、原典では未確認）。つまり導手 $2^k$ の奇指標の $B_1$ のノルムを見ている。「$\ell\mid h_k^+$ なら $\ell$ は $\mathbb{Q}(\zeta_{2^k})$ の相対類数を割る」という定理は、$\ell\ne2$ では知られていない（鏡映定理が結ぶのは $\mathbb{Q}(\zeta_{2^k\ell})$ の中の成分）。
  - 問題 3：文献の記述が合わない。Table 1 は $k=9,\dots,12$ を「Fukuda–Komatsu、条件つき」としているが、Fukuda–Komatsu の結果は「$10^9$ 未満の素数は $h_n$ を割らない」（無条件、すべての $n$）で、類数 1 の証明ではない。GRH 下で知られているのは Miller の導手 512（$k=9$）まで（MathSciNet の MR3304606 のレビュー）。$k=10,11,12$ は GRH 下でも知られていない。
  - 問題 1 の補足（PDF で確認）：Lemma 3.3 の証明の最後は「ノルムの $\ell$ 進付値が等しい。完全な証明は Stickelberger 元と岩澤理論を使うが、Weber 予想の文脈では標準的な結果（Fukuda–Komatsu）」とあるだけ。直後に「$\chi^{-1}$ は導手 $2^k$、位数 $2^{k-2}$ なので Lemma 3.4 の評価が使える」とあるが、Lemma 3.4 は奇指標についての補題で、$\chi^{-1}$ は偶指標である。Theorem 3.2 の証明は、Herbrand の定理が $\mathbb{Q}(\zeta_\ell)$ についてのものだと自分で述べたうえで、「高い円分体の場合を使った」と書くだけである。
  - 問題 4：Theorem 3.1 (ii)「奇素数 $\ell\mid h_k^+$ なら $\ell\equiv\pm1 \bmod 2^{k-1}$」は、証明が引用（Fukuda–Komatsu の §3）に任されている。この形の定理があるかは未確認。私の記憶にある既知の結果は、小さい法（16 や 64）についてのもので、$2^{k-1}$（$k=12$ で 2048）ではない。
  - 問題 3 の補足：序文は「Miller は $K_8^+$ の単数群全体を計算し Sinnott の指数公式で証明した」「Fukuda–Komatsu は GRH の下で Bach の上界により $k\le12$ まで確かめ、Schoof が独立に確かめた」と書く。前者は Miller の論文（原論文で確認済み。判別式の評価と単項素イデアルによる）の内容と合わない。後者に当たる結果は、確認した文献の範囲にはない。
  - そのほか：計算のコードも出力も示されていない。$k=12$ では 726 桁の数の素因数分解が要るとされ、本文が挙げる素因数分解の能力（約 250 桁）を超える。証明の文が途中で切れている箇所がある。
  - 外部の評価：同じ著者の連作（Part IV まで）を取り上げた解説記事（postquantum.com、2026-05-25）は、Part I には具体的な批判も確認も載せていない。査読は受けていない。数論の側からの確認は見つからなかった。
