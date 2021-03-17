# R, BUGS, Stanによる階層モデルのあてはめ

日本生態学会第68回大会における自由集会[「階層モデリングの実践：統計モデルを生態研究と管理・保全に活用する」](https://esj.ne.jp/meeting/abst/68/W02.html)での講演[「R, BUGS, Stanによる階層モデルのあてはめ」](https://esj.ne.jp/meeting/abst/68/W02-3.html)の発表資料および実行コード (R Markdown形式)

## 発表ファイル

- figshareで公開しています。[doi:10.6084/m9.figshare.14229572.v1](https://doi.org/10.6084/m9.figshare.14229572.v1)

## ファイル

- esj68.Rmd: 発表スライドのR Markdownファイル
    - nmix.txt: *N*混合モデルのBUGSモデルファイル
    - nmix.stan: *N*混合モデルのStanモデルファイル
    - nmix_rs.stan: *N*混合モデルのStanモデルファイル（`reduce_sum`による並列化。今回は使用せず）
    - nmix_bp.stan: *N*混合モデルのStanモデルファイル（2変量ポアソン分布によるモデル化）
- occ.Rmd: サイト占有モデルのR Markdownファイル（発表で割愛した分）
    - occ.txt:サイト占有モデルのBUGSモデルファイル
    - occ.stan: サイト占有モデルのStanモデルファイル

## 参考文献

-  Dennis E.B., Morgan B.J.T., Ridout M.S. (2015) Computational aspects of N-mixture models. Biometrics 71:237–246. [doi:10.1111/biom.12246](https://doi.org/10.1111/biom.12246)
- Fiske I., Chandler R. (2011) unmarked: An R Package for Fitting Hierarchical Models of Wildlife Occurrence and Abundance. Journal of Statistical Software 43(10): 1–23. [http://www.jstatsoft.org/v43/i10/](http://www.jstatsoft.org/v43/i10/)
- Kéry M., Royle J.A. (2015) Applied Hierarchical Modeling in Ecology: Analysis of distribution, abundance and species richness in R and BUGS: Volume 1:Prelude and Static Models. Academic Press. (日本語訳: 深谷肇一・飯島勇人・伊東宏樹(監訳), 飯島勇人・伊東宏樹・奥田武弘・長田穣・川森愛・柴田泰宙・高木俊・辰巳晋一・仁科一哉・深澤圭太・深谷肇一・正木隆(訳) (2021) 「[生態学のための階層モデリング―RとBUGSによる分布・個体数量・種の豊かさの統計解析―](https://www.kyoritsu-pub.co.jp/bookdetail/9784320058149)」 共立出版)
- Kéry M., Royle A., Meredith M. (2020) AHMbook: Functions and Data for the Book 'Applied Hierarchical Modeling in Ecology' Vols 1 and 2. [https://CRAN.R-project.org/package=AHMbook](https://CRAN.R-project.org/package=AHMbook)
- Plummer M. (2017) JAGS Version 4.3.0 user mannual. [https://mcmc-jags.sourceforge.io/](https://mcmc-jags.sourceforge.io/)
- Stan Development Team. (2020) Stan Modeling Language Users Guide and Reference Manual, 2.25.0. [https://mc-stan.org/](https://mc-stan.org/)
