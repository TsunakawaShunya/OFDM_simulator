# to do
- 変調方式の確認(IEEE802.11aとかの)
- 各パラメータの確認
- 平均二乗誤差の計算が怪しいから確認
- simulator.hの193行目からのインパルス応答h_の作り方に式(22)の時間方向の相関を導入する

# 進捗
- 完璧な推定値を用いた等化の場合，理論値を一致することを確認

# 実行方法
```
    # コンパイル
    $ g++ -std=c++17 -I/Volumes/USB1/eigen-3.4.0 -I/Users/shunya/Downloads/boost_1_87_0 -o power_ber power_ber.cpp
    # 実行
    $ ./power_ber
```

# ChatGPTへの質問
OFDMにおける時変伝送路でのビット誤り率をシミュレーションするコードをC++で書いています．下記に理論の詳細を記載するので，下記をよく読んでC++でEigenライブラリを使用してインパする応答h_{l,q}を生成するメソッドを完成させてください．よく熟考し，一歩一歩着実に実装してみてください．

ここから理論の説明に入ります．
時間方向でl番目の第qパスにおけるインパルス応答をh_{l,q}とするとh_{l,q}は
E{h_{l,q}}=0,E{h_{q_1,l_1}h_{q_2,l_2}^{*}}=\quzai_{q_1}^2\delta_{q_1,q_2}J_0(2\pi\abs{l_1-l_2}T_s f_d)を満たす複素ガウス確率の標本です．
ただし，遅延プロファイル\quzai_qについて，\sum^{Q-1}_{q}\quzai_q^2=1を満たしています．

今，\eqvec{h_q}=[h_{0,q},\cdots,h_{L-1,q}]^Tとするとき，\eqvec{h_q}=\eqmat{A_q}\eqvce{x}によって\eqvec{h_q}を生成する方法について考えています．
\eqvec{x}\in\mathbb{C}^LはE{\eqvec{x}}=\eqvec{0},E{\eqvec{x}\eqvec{x}^H}=\eqmat{I}となる複素ガウス確率変数の標本．ただし\eqmat{I}は単位行列です．

\eqmat{A_q}は導出の過程は省略するが，下記のように求められる．
1. E{h_q h_q^H}=\eqmat{R}とおく．Rはh_qにおける共分散行列で(l_1, l_2)成分はJ_0(2\pi\abs{l_1-l_2}T_s f_d)である．
2. \eqmat{R} = \eqmat{U_q}\eqmat{\Lambda_q}\eqmat{U_q}^{-1}と固有値分解する．\eqmat{U_q}はユニタリー行列で，\eqmat{\Lambda_q}は固有値ベクトルを対角行列にした行列である．
3. \eqmat{A_q} = \eqmat{U_q}\sqrt{\eqmat{\Lambda_q}}

上記のように求められた\eqmat{A_q}と複素ガウス確率変数の標本として乱数で与えられる\eqvec{x}を使って\eqvec{h_q}=\eqmat{A_q}\eqvec{x}を生成します．