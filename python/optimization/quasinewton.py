# -*- coding: utf-8 -*-
"""!
@file quasinewton.py

@brief 準ニュートン法

@author Makoto Fujisawa
@date 2019-07 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *
from rx_funcs import *

SQRT5 = 2.23606797749979


# -----------------------------------------------------------------------------
# 最適化関数(最小値探索)
# -----------------------------------------------------------------------------
def goldensection(func, x0, d, xl, xr, max_iter, eps):
    """!
    黄金分割探索法(golden-section method)
     - 多変数で探索方向が決まっている場合
     - BFGSでの探索方向に対する係数α決定のために用いる
    @param[in] func 関数値を与える関数
    @param[in] x0 探索初期地点ベクトル
    @param[in] d 探索方向ベクトル
    @param[in] xl,xr 初期探索範囲
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解ans, 実際の反復数, 実際の誤差)
    """
    eta = 2.0/(SQRT5+3.0)
    beta = xr-xl
    tau = eta*beta
    n = len(x0)

    x = zeros(4)
    x[0] = xl
    x[1] = xl+tau
    x[2] = xr-tau
    x[3] = xr

    # 探索方向上の2点x1,x2での関数値を計算
    x1 = zeros(n)
    x2 = zeros(n)
    for i in range(n):
        x1[i] = x0[i]+x[1]*d[i]
        x2[i] = x0[i]+x[2]*d[i]
    f1 = func(x1)
    f2 = func(x2)

    k = 0
    for k in range(max_iter):
        if f1 < f2:  # 区間[x2,x3]に最小値なし
            # 分割区間更新
            tx2 = x[2]
            x[2] = x[1]
            x[1] = tx2-tau
            x[3] = tx2

            # 関数値の更新
            for i in range(n):
                x1[i] = x0[i]+x[1]*d[i]
            f2 = f1
            f1 = func(x1)
        else:  # 区間[x0,x1]に最小値なし
            # 分割区間更新
            tx1 = x[1]
            x[1] = x[2]
            x[2] = tx1+tau
            x[0] = tx1

            # 関数値の更新
            for i in range(n):
                x2[i] = x0[i]+x[2]*d[i]
            f1 = f2
            f2 = func(x2)

        beta -= tau      # 新しい[x0,x3]の長さ(現在の長さから[x0,x1](or[x2,x3])の長さtauを引く)
        tau = eta*beta   # 新しい[x0,x1](or[x2,x3])の長さ

        if beta < eps:
            break

    return 0, x[1], k, beta


def quasinewton_bfgs(func, dfunc, x0, max_iter, eps):
    """!
    BFGS法(準ニュートン法の一種, Broyden-Fletcher-Goldfarb-Shanno method)
    @param[in] func 関数値を与える関数
    @param[in] dfunc 関数勾配値を与える関数
    @param[in] x0 初期探索地点
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解xout, 実際の反復数, 実際の誤差)
    """
    n = len(x0)
    H = zeros2(n, n)  # ヘッセ行列の逆行列H^-1
    for i in range(n):
        H[i][i] = 1.0

    x = list(x0)
    g = dfunc(x)
    gnorm = 1.0

    k = 0
    for k in range(max_iter):
        print("f(" + fmt(x) + ") = " + fmt(func(x)))
        # print(fmt(x) + ", " + fmt(gnorm))

        # 探索方向dの計算(d = -H∇f)
        p = zeros(n)
        for i in range(n):
            for j in range(n):
                p[i] -= H[i][j]*g[j]

        # 黄金分割法で探索方向d上の最小値までの距離を調べてalphaとして設定
        gmax_iter = 20
        geps = 1e-4
        ret, alpha, gmax_iter, geps = goldensection(func, x, p, 0.0, 3.0, gmax_iter, geps)

        # xの更新
        s = zeros(n)
        for i in range(n):
            s[i] = alpha*p[i]
            x[i] += s[i]

        # yの計算
        gprev = g
        g = dfunc(x)
        y = zeros(n)
        for i in range(n):
            y[i] = g[i]-gprev[i]

        gnorm = 0.0
        for i in range(n):
            gnorm += g[i]*g[i]
        if math.sqrt(gnorm) < eps:
            break  # 勾配ベクトル∇fの大きさで収束判定

        # H*yの計算(行列×縦ベクトル⇒縦ベクトル)
        Hy = zeros(n)
        for i in range(n):
            for j in range(n):
                Hy[i] += H[i][j]*y[j]

        # (s^T)*yと(y^T)*(H*y)の計算(どちらもベクトル同士の内積)
        sy = 0.0
        yHy = 0.0
        for i in range(n):
            sy += s[i]*y[i]
            yHy += y[i]*Hy[i]

        # H^-1の更新
        for i in range(n):
            for j in range(n):
                H[i][j] += s[i]*s[j]/sy + yHy*s[i]*s[j]/(sy*sy) - Hy[i]*s[j]/sy - s[i]*Hy[j]/sy

    return 0, x, k, gnorm


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 初期値を格納する配列x0の宣言
    x0 = zeros(2, -1.0)

    # 値を返す関数と導関数を返す関数の指定
    func = Func4
    dfunc = DFunc4

    # 準ニュートン法で数値解(関数funcが最小値を取るx)を求める
    max_iter = 100
    eps = 1e-6
    ret, x, max_iter, eps = quasinewton_bfgs(func, dfunc, x0, max_iter, eps)

    # 最終計算結果の表示
    print("(x,y) = (" + fmt(x[0]) + "," + fmt(x[1]) + "), " + " f(x,y) = " + fmt(func(x)))

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
