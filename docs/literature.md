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
| $k_{3^3}$ | なし | Masley (1978) |
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
