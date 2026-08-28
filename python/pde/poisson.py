# -*- coding: utf-8 -*-
"""!
@file poisson.py

@brief ポアソン方程式のソルバー

@author Makoto Fujisawa
@date 2019-11 (Python版)
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

# 境界値
alpha = 0.0
beta = 0.0


# -----------------------------------------------------------------------------
# ポアソン方程式のソルバー
# -----------------------------------------------------------------------------
def cg_solver(A, b, x, n, max_iter, eps):
    """!
    共役勾配法によりA・x=bを解く
    @param[in] A n×n正値対称行列
    @param[in] b 右辺ベクトル
    @param[out] x 結果ベクトル(リストなので破壊的に更新される)
    @param[in] n 行列の大きさ
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    if n <= 0:
        return 1, max_iter, eps

    r = zeros(n)
    p = zeros(n)
    y = zeros(n)
    for i in range(n):
        x[i] = 0.0

    # 第0近似解に対する残差の計算
    for i in range(n):
        ax = 0.0
        for j in range(n):
            ax += A[i][j]*x[j]
        r[i] = b[i]-ax
        p[i] = r[i]

    rr0 = dot(r, r, n)

    e = 0.0
    k = 0
    for k in range(max_iter):
        # y = AP の計算
        for i in range(n):
            y[i] = dot(A[i], p, n)

        # alpha = r*r/(P*AP)の計算
        a = rr0/dot(p, y, n)

        # 解x、残差rの更新
        for i in range(n):
            x[i] += a*p[i]
            r[i] -= a*y[i]

        # (r*r)_(k+1)の計算
        rr1 = dot(r, r, n)

        # 収束判定 (||r||<=eps)
        e = math.sqrt(rr1)
        if e < eps:
            k += 1
            break

        # βの計算とPの更新
        bt = rr1/rr0
        for i in range(n):
            p[i] = r[i]+bt*p[i]

        # (r*r)_(k+1)を次のステップのために確保しておく
        rr0 = rr1

    return 0, k+1, e


def setbc_d(f, n):
    """!
    境界条件設定(ディリクレ境界条件)
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    """
    f[0] = alpha
    f[n] = beta


def poisson1d_central(f, x0, x1, n, g):
    """!
    中心差分で1次元ポアソン方程式を解く
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] x0,x1 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @param[in] g 右辺項を与える関数
    @return 1
    """
    h = (x1-x0)/n  # 空間刻み幅

    # 線形システムの係数行列と右辺項ベクトルの計算
    b = zeros(n-1)
    x = zeros(n-1)
    A = zeros2(n-1, n-1)  # n-1×n-1の2次元配列
    for i in range(1, n-2):
        xi = x0+h*(i+1)
        A[i][i-1] = -1
        A[i][i] = 2
        A[i][i+1] = -1
        b[i] = -h*h*g(xi)
    A[0][0] = 2
    A[0][1] = -1
    b[0] = -h*h*g(x0+h)+f[0]
    A[n-2][n-3] = -1
    A[n-2][n-2] = 2
    b[n-2] = -h*h*g(x1-h)+f[n]
    for i in range(n-1):
        x[i] = f[i+1]

    # CG法で線形システムを解く
    max_iter = 100
    eps = 1e-6
    ret, max_iter, eps = cg_solver(A, b, x, n-1, max_iter, eps)
    print("cg solver : max_iter = " + str(max_iter) + ", eps = " + fmt(eps))

    # 結果を配列fに戻す
    for i in range(n-1):
        f[i+1] = x[i]
    setbc_d(f, n)  # 境界条件

    return 1


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    global alpha, beta

    func = FuncPdeX3
    a, b = 0.0, 1.0
    alpha, beta = 0.0, 0.0
    n = 10
    func_t = FuncPdeX3T

    f = zeros(n+1)
    setbc_d(f, n)

    # 中心差分+CG法でポアソン方程式を解く
    poisson1d_central(f, a, b, n, func)

    # 結果の出力
    h = (b-a)/n
    avg_error = 0.0
    for i in range(n+1):
        x = a+h*i

        e = abs(f[i]-func_t(x))
        print("f(" + fmt(x) + ") = " + fmt(f[i]) + ", error = " + fmt(e))
        # print(fmt(x) + ", " + fmt(f[i]) + ", " + fmt(func_t(x)))
        avg_error += e

    print("avg. error = " + fmt(avg_error/(n-1)))

    return 0


if __name__ == "__main__":
    main()
