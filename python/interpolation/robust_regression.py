# -*- coding: utf-8 -*-
"""!
@file robust_regression.py

@brief ロバスト推定(外れ値を含むデータへの多項式フィッティング)
       - 最小二乗法(比較用)
       - 最小メディアン法(LMedS : Least Median of Squares)
       - RANSAC(RANdom SAmple Consensus)
       - M推定(Huber, Tukey biweight．反復重み付き最小二乗法(IRLS)で解く)

@author Makoto Fujisawa
@date 2026-11 (Python版)

C++版(src/interpolation/robust_regression.cpp)と同じ手順．
乱数生成器が違う(C++版は std::mt19937 と標準の分布関数)ので，データや結果の数値はC++版と少し違う．
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math
import random

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *
from rx_funcs import *

# このスクリプトが置かれているフォルダ(データファイルの入出力先の基準にする)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# -----------------------------------------------------------------------------
# LU分解(最小自乗法での係数計算に用いる．leastsquare.pyと同じ)
# -----------------------------------------------------------------------------
def LUDecomp(A, n):
    """!
    LU分解(ピボット交換なし)
     - 行列A(n×n)を下三角行列(L)と上三角行列(U)に分解する
     - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
     - LとUを一つの行列にまとめた形で結果を返す
    @param[inout] A n×nの係数行列．LU分解した結果を格納する．
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    for i in range(n):
        # l_ijの計算(i >= j)
        for j in range(i+1):
            lu = A[i][j]
            for k in range(j):
                lu -= A[i][k]*A[k][j]  # l_ik * u_kj
            A[i][j] = lu

        # u_ijの計算(i < j)
        for j in range(i+1, n):
            lu = A[i][j]
            for k in range(i):
                lu -= A[i][k]*A[k][j]  # l_ik * u_kj
            A[i][j] = lu/A[i][i]

    return 0


def LUSolver(A, b, x, n):
    """!
    LU分解した行列A(n×n)から前進代入・後退代入によりA・x=bを解く
    @param[in] A LU分解された行列
    @param[in] b 右辺ベクトル
    @param[out] x 結果ベクトル
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    # 前進代入(forward substitution)
    #  LY=bからYを計算
    for i in range(n):
        bly = b[i]
        for j in range(i):
            bly -= A[i][j]*x[j]
        x[i] = bly/A[i][i]

    # 後退代入(back substitution)
    #  UX=YからXを計算
    for i in range(n-1, -1, -1):
        yux = x[i]
        for j in range(i+1, n):
            yux -= A[i][j]*x[j]
        x[i] = yux

    return 0


# -----------------------------------------------------------------------------
# 共通の関数
# -----------------------------------------------------------------------------
def polynominal1d(x, c, n):
    """!
    1次元n次多項式を計算
    @param[in] x 値が必要な位置
    @param[in] c 多項式の係数
    @param[in] n 多項式の次数
    @return 多項式を計算した結果
    """
    fx = 0
    xb = 1
    for i in range(n+1):
        fx += c[i]*xb
        xb *= x
    return fx


def weighted_leastsquare(fi, xi, w, m, n, c):
    """!
    重み付き最小二乗法による多項式フィッティング
     - 1次元データ(xi,fi)に対して，Σ w_k (f(x_k) - fi_k)^2 を最小にする1次元n次多項式の係数を求める
     - w_k = 1 なら通常の最小二乗法(leastsquare.py の leastsquare_interpolation と同じ)
     - w_k = 0 のデータは使われない(RANSAC，最小メディアン法でインライアだけを使うときなど)
    @param[in] fi 関数値を格納した配列
    @param[in] xi 関数値に対応する位置を格納した配列
    @param[in] w 各データの重み(Noneまたは空の場合は全て1)
    @param[in] m データ数
    @param[in] n フィッティングする多項式の次数
    @param[out] c 係数ベクトル(リストなので破壊的に更新される)
    @return 0:成功,1:失敗(重みが正のデータ数が係数の数より少ない)
    """
    dim_b = n+1  # 基底ベクトルの次元数(1次元の場合)

    # 重みが正のデータ数を確認
    cnt = sum(1 for k in range(m) if not w or w[k] > 0)
    if cnt < dim_b:
        return 1

    # Ac=y, bは基底ベクトル
    del c[:]
    c.extend(zeros(dim_b))  # 結果の係数ベクトル
    b = zeros(dim_b)        # 基底ベクトル
    y = zeros(dim_b)        # 右辺項
    A = zeros2(dim_b, dim_b)

    # 係数行列Aと右辺項yの計算
    for k in range(m):  # データ分だけ反復
        wk = 1.0 if not w else w[k]
        if wk <= 0:
            continue

        # 多項式の基底ベクトルの計算
        xb = 1
        for i in range(dim_b):
            b[i] = xb
            xb *= xi[k]

        # 係数行列Aと右辺項bの計算(重みを掛ける)
        for i in range(dim_b):
            for j in range(dim_b):
                A[i][j] += wk*b[i]*b[j]
            y[i] += wk*b[i]*fi[k]

    # Ac=yをLU分解で解く
    LUDecomp(A, dim_b)
    LUSolver(A, y, c, dim_b)

    return 0


def residuals(fi, xi, m, c, n):
    """!
    各データの残差 r_k = fi_k - f(x_k) を計算
    @param[in] fi,xi データ
    @param[in] m データ数
    @param[in] c 多項式の係数
    @param[in] n 多項式の次数
    @return 残差のリスト
    """
    return [fi[k]-polynominal1d(xi[k], c, n) for k in range(m)]


def median(v):
    """!
    中央値(メディアン)を計算
    @param[in] v 値の配列(コピーして並べ替えるので元の配列は変わらない)
    @return 中央値
    """
    s = sorted(v)
    m = len(s)
    return s[m//2] if m % 2 == 1 else 0.5*(s[m//2-1]+s[m//2])


def random_sample(m, s, rng):
    """!
    m個のデータから重複なしでs個をランダムに選ぶ
    @param[in] m データ数
    @param[in] s 選ぶ個数
    @param[in] rng 乱数生成器(random.Random)
    @return 選んだデータのインデックスのリスト
    """
    idx = []
    while len(idx) < s:
        k = rng.randint(0, m-1)
        if k not in idx:
            idx.append(k)
    return idx


# -----------------------------------------------------------------------------
# ロバスト推定法
# -----------------------------------------------------------------------------
def lmeds(fi, xi, m, n, max_trial, c, inlier, seed=1):
    """!
    最小メディアン法(LMedS : Least Median of Squares)
     - 残差の2乗の「和」ではなく「中央値」を最小にする係数を求める
     - 中央値はデータの半分までが外れ値でも影響を受けにくい(崩壊点50%)
     - 中央値の最小化は微分で解けないので，ランダムに選んだn+1点を通る多項式を候補として，
       中央値が最小になるものを探す
     - 最後に，求めた中央値から残差のばらつきσを推定して，|r| <= 2.5σ のデータ(インライア)だけで最小二乗法をやり直す
    @param[in] fi,xi データ
    @param[in] m データ数
    @param[in] n フィッティングする多項式の次数
    @param[in] max_trial 候補を作る回数
    @param[out] c 係数ベクトル(リストなので破壊的に更新される)
    @param[out] inlier 各データがインライアなら1，外れ値なら0(リストなので破壊的に更新される)
    @param[in] seed 乱数のシード
    @return 最小となった残差2乗の中央値
    """
    rng = random.Random(seed)
    s = n+1  # 多項式を決めるのに必要なデータ数
    best_med = 1e30
    best_c = None
    cs = []
    for t in range(max_trial):
        # ランダムに選んだs点だけを使って多項式を決める(s点を通る多項式)
        idx = random_sample(m, s, rng)
        w = zeros(m)
        for k in idx:
            w[k] = 1.0
        if weighted_leastsquare(fi, xi, w, m, n, cs):
            continue

        # 全データの残差の2乗の中央値
        r = residuals(fi, xi, m, cs, n)
        med = median([rk*rk for rk in r])
        if med < best_med:
            best_med = med
            best_c = list(cs)

    # ロバストな標準偏差の推定(Rousseeuw & Leroy の式)と，それを使ったインライアの判定
    sigma = 1.4826*(1.0+5.0/(m-s))*math.sqrt(best_med)
    r = residuals(fi, xi, m, best_c, n)
    del inlier[:]
    inlier.extend([1 if abs(rk) <= 2.5*sigma else 0 for rk in r])

    # インライアだけで最小二乗法
    weighted_leastsquare(fi, xi, [float(v) for v in inlier], m, n, c)

    return best_med


def ransac(fi, xi, m, n, max_trial, threshold, c, inlier, seed=1):
    """!
    RANSAC(RANdom SAmple Consensus)
     - ランダムに選んだn+1点を通る多項式を候補として，残差が閾値以下のデータ(インライア)の数が最大の候補を選ぶ
     - 最後にインライアだけで最小二乗法をやり直す
     - 閾値(どのくらいのずれまでを正しいデータとみなすか)はあらかじめ与える必要がある
    @param[in] fi,xi データ
    @param[in] m データ数
    @param[in] n フィッティングする多項式の次数
    @param[in] max_trial 候補を作る回数
    @param[in] threshold インライアとみなす残差の大きさ
    @param[out] c 係数ベクトル(リストなので破壊的に更新される)
    @param[out] inlier 各データがインライアなら1，外れ値なら0(リストなので破壊的に更新される)
    @param[in] seed 乱数のシード
    @return インライアの数
    """
    rng = random.Random(seed)
    s = n+1  # 多項式を決めるのに必要なデータ数
    best_cnt = -1
    best_err = 1e30
    best_c = None
    cs = []
    for t in range(max_trial):
        # ランダムに選んだs点だけを使って多項式を決める
        idx = random_sample(m, s, rng)
        w = zeros(m)
        for k in idx:
            w[k] = 1.0
        if weighted_leastsquare(fi, xi, w, m, n, cs):
            continue

        # インライアの数を数える(同じ数なら残差の2乗和が小さい方)
        r = residuals(fi, xi, m, cs, n)
        cnt = 0
        err = 0.0
        for rk in r:
            if abs(rk) <= threshold:
                cnt += 1
                err += rk*rk
        if cnt > best_cnt or (cnt == best_cnt and err < best_err):
            best_cnt = cnt
            best_err = err
            best_c = list(cs)

    # インライアだけで最小二乗法
    r = residuals(fi, xi, m, best_c, n)
    del inlier[:]
    inlier.extend([1 if abs(rk) <= threshold else 0 for rk in r])
    weighted_leastsquare(fi, xi, [float(v) for v in inlier], m, n, c)

    return best_cnt


def weight_huber(u, k=1.345):
    """!
    M推定で使う重み関数(Huber)
     - |u| <= k では最小二乗法と同じ(重み1)，それより大きいと重みが 1/|u| に比例して小さくなる
    @param[in] u 残差をスケールσで割った値
    @param[in] k 調整定数(正規分布のとき効率95%になる 1.345 がよく使われる)
    """
    return 1.0 if abs(u) <= k else k/abs(u)


def weight_tukey(u, k=4.685):
    """!
    M推定で使う重み関数(Tukey biweight)
     - |u| <= k では (1-(u/k)^2)^2，それより大きいと重み0(完全に無視する)
    @param[in] u 残差をスケールσで割った値
    @param[in] k 調整定数(正規分布のとき効率95%になる 4.685 がよく使われる)
    """
    if abs(u) > k:
        return 0.0
    a = 1.0-(u/k)*(u/k)
    return a*a


def m_estimation(fi, xi, m, n, wfunc, k, c, w, max_iter, eps):
    """!
    M推定(反復重み付き最小二乗法 IRLS : Iteratively Reweighted Least Squares)
     - Σ ρ(r_k/σ) を最小にする係数を求める(ρ(u) = u^2/2 なら最小二乗法)
     - 残差が大きいデータほど重みを小さくした重み付き最小二乗法を，重みを更新しながら繰り返す
     - スケールσは残差の中央絶対偏差(MAD)から推定する : σ = median(|r - median(r)|)/0.6745
     - 初期値は通常の最小二乗法(外れ値が多いと初期値が悪くなり，うまくいかないことがある)
    @param[in] fi,xi データ
    @param[in] m データ数
    @param[in] n フィッティングする多項式の次数
    @param[in] wfunc 重み関数(weight_huber, weight_tukey)
    @param[in] k 重み関数の調整定数(Huber : 1.345，Tukey : 4.685)
    @param[out] c 係数ベクトル(リストなので破壊的に更新される)
    @param[out] w 最終的な各データの重み(リストなので破壊的に更新される)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (0:成功/1:収束しなかった, 反復回数, 最後の係数の変化量)
    """
    # 初期値 : 通常の最小二乗法
    del w[:]
    w.extend([1.0]*m)
    weighted_leastsquare(fi, xi, w, m, n, c)

    ret = 1
    e = 0.0
    it = 0
    while it < max_iter:
        # 残差とスケールσ(MAD)
        r = residuals(fi, xi, m, c, n)
        med = median(r)
        sigma = median([abs(rk-med) for rk in r])/0.6745
        if sigma < 1e-12:  # 残差がほぼ0(全データにぴったり合っている)
            ret = 0
            break

        # 重みを更新して重み付き最小二乗法
        for i in range(m):
            w[i] = wfunc(r[i]/sigma, k)
        c_old = list(c)
        weighted_leastsquare(fi, xi, w, m, n, c)

        # 係数の変化量で収束判定
        e = math.sqrt(sum((c[i]-c_old[i])**2 for i in range(n+1)))
        it += 1
        if e < eps:
            ret = 0
            break

    return ret, it, e


# -----------------------------------------------------------------------------
# データの作成と結果の表示
# -----------------------------------------------------------------------------
def MakeSamplingPointsWithOutliers(x0, x1, m, func, nwidth, num_outlier, owidth, seed=12345):
    """!
    外れ値を含むサンプリング点(データ点)の生成(1次元)
     - 等間隔の点で関数値に[-0.5nwidth, 0.5nwidth]の一様ノイズを加え，
       そのうち num_outlier 個のデータに [owidth/2, owidth] の大きさの外れ値(符号はランダム)を加える
     - 結果を再現できるように乱数のシードを指定する
    @param[in] x0,x1 サンプリング範囲
    @param[in] m データ数
    @param[in] func 関数値を与える関数
    @param[in] nwidth ノイズの幅
    @param[in] num_outlier 外れ値の数
    @param[in] owidth 外れ値の大きさ
    @param[in] seed 乱数のシード
    @return (xi, yi, outlier) outlierは外れ値を加えたデータなら1
    """
    rng = random.Random(seed)
    xi = [x0+(x1-x0)*k/(m-1.0) for k in range(m)]
    yi = [func(x)+(rng.random()-0.5)*nwidth for x in xi]
    outlier = [0]*m
    for k in random_sample(m, num_outlier, rng):
        o = owidth*(0.5+0.5*rng.random())
        yi[k] += -o if rng.random() < 0.5 else o
        outlier[k] = 1
    return xi, yi, outlier


def PrintResult(name, c, n, func, x, x0, x1):
    """!
    フィッティング結果の表示
     - 係数，位置xでの値と真値との誤差，[x0,x1]での真の関数との最大誤差
    """
    fx = polynominal1d(x, c, n)
    emax = max(abs(polynominal1d(x0+(x1-x0)*i/200.0, c, n)-func(x0+(x1-x0)*i/200.0)) for i in range(201))
    print(name)
    print("  c = " + ", ".join(fmt(ci) for ci in c))
    print("  f(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-func(x)))
          + ",  max error in [" + fmt(x0) + "," + fmt(x1) + "] = " + fmt(emax))


def PrintInlier(inlier, outlier):
    """!
    インライア/外れ値の判定結果の表示(外れ値と判定されたデータの番号)
    """
    # inlier=1かつoutlier=1，またはinlier=0かつoutlier=0なら判定ミス
    wrong = sum(1 for k in range(len(inlier)) if inlier[k] == outlier[k])
    rej = " ".join(str(k) for k in range(len(inlier)) if not inlier[k])
    print("  rejected : " + rej + "  (misclassified : " + str(wrong) + ")")


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    func = FuncExp
    x0, x1 = 0, 1

    d = 3             # フィッティングする多項式の次数
    m = 21            # データ数
    num_outlier = 5   # 外れ値の数
    x = 0.5           # 値を調べる位置
    set_precision(6)

    # 外れ値を含むデータ(ノイズ幅0.05，外れ値の大きさ0.5～1.0)
    xi, yi, outlier = MakeSamplingPointsWithOutliers(x0, x1, m, func, 0.05, num_outlier, 1.0)
    OutputSamplingPoints(xi, yi)  # サンプリング点の画面出力
    print("outliers : " + " ".join(str(k) for k in range(m) if outlier[k]))
    print("ground truth f(" + fmt(x) + ") = " + fmt(func(x)))
    print()

    c, w, inlier = [], [], []

    # 通常の最小二乗法(比較用)
    weighted_leastsquare(yi, xi, None, m, d, c)
    PrintResult("[least squares]", c, d, func, x, x0, x1)

    # 最小メディアン法
    med = lmeds(yi, xi, m, d, 500, c, inlier)
    PrintResult("[LMedS] (median of r^2 = " + fmt(med) + ")", c, d, func, x, x0, x1)
    PrintInlier(inlier, outlier)

    # RANSAC(閾値はノイズの大きさから決める)
    cnt = ransac(yi, xi, m, d, 500, 0.1, c, inlier)
    PrintResult("[RANSAC] (inliers = " + str(cnt) + ")", c, d, func, x, x0, x1)
    PrintInlier(inlier, outlier)

    # M推定(Huber)
    ret, it, e = m_estimation(yi, xi, m, d, weight_huber, 1.345, c, w, 100, 1e-8)
    PrintResult("[M-estimation (Huber)] (iter = " + str(it) + ")", c, d, func, x, x0, x1)

    # M推定(Tukey biweight)
    ret, it, e = m_estimation(yi, xi, m, d, weight_tukey, 4.685, c, w, 100, 1e-8)
    PrintResult("[M-estimation (Tukey)] (iter = " + str(it) + ")", c, d, func, x, x0, x1)
    print("  weights : " + ", ".join(fmt(wk) for wk in w))

    # # グラフ描画用にデータ出力
    # dat = os.path.join(SCRIPT_DIR, "dat")
    # if not os.path.exists(dat): os.makedirs(dat)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "robust_data.txt"))  # サンプリング点のファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, func, os.path.join(dat, "robust_ground_truth.txt"))
    # weighted_leastsquare(yi, xi, None, m, d, c)
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, d), os.path.join(dat, "robust_ls.txt"))
    # lmeds(yi, xi, m, d, 500, c, inlier)
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, d), os.path.join(dat, "robust_lmeds.txt"))
    # ransac(yi, xi, m, d, 500, 0.1, c, inlier)
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, d), os.path.join(dat, "robust_ransac.txt"))
    # m_estimation(yi, xi, m, d, weight_tukey, 4.685, c, w, 100, 1e-8)
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, d), os.path.join(dat, "robust_mest.txt"))

    return 0


if __name__ == "__main__":
    main()
