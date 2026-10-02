# -*- coding: utf-8 -*-
"""!
@file gradientdecent.py

@brief 最急降下法

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


# -----------------------------------------------------------------------------
# 最適化関数(最小値探索)
# -----------------------------------------------------------------------------
def gradientdecent(dfunc, x0, alpha, max_iter, eps):
    """!
    最急降下法(gradient decent method, steepest decent method)
    @param[in] dfunc 関数勾配値を与える関数
    @param[in] x0 初期探索地点
    @param[in] alpha 移動量の係数
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解xout, 実際の反復数, 実際の誤差)
    """
    x = x0
    dx = 0.0

    k = 0
    for k in range(max_iter):
        print("x(" + str(k) + ") = " + fmt(x))
        dx = alpha*dfunc(x)  # 関数の勾配に係数αを掛ける
        x -= dx              # 勾配方向に探索点を移動させる
        if abs(dx) < eps:    # 勾配の大きさで収束判定
            break

    return 0, x, k, abs(dx)


def gradientdecent_nd(dfunc, x0, alpha, max_iter, eps):
    """!
    最急降下法(gradient decent method, steepest decent method)
     - n次元(n>1)版
     - C++版ではオーバーロードで同じ名前にしていたが，
       Pythonには関数のオーバーロードがないので名前を変えている
    @param[in] dfunc 関数勾配値を与える関数
    @param[in] x0 初期探索地点(リスト)
    @param[in] alpha 移動量の係数
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解xout, 実際の反復数, 実際の誤差)
    """
    n = len(x0)
    x = list(x0)

    norm_dx = 1.0  # 勾配ベクトルのノルム(収束判定用)
    k = 0
    for k in range(max_iter):
        # print(str(k) + " : " + fmt(x) + ", e = " + fmt(norm_dx))
        print(fmt(x) + ", " + fmt(norm_dx))

        dx = dfunc(x)
        norm_dx = 0.0
        for i in range(n):
            dx[i] *= alpha    # 関数の勾配ベクトルに係数αを掛ける
            x[i] -= dx[i]     # 勾配方向に探索点を移動させる
            norm_dx += dx[i]*dx[i]

        # 勾配ベクトルのノルムで収束判定
        norm_dx = math.sqrt(norm_dx)
        if norm_dx < eps:
            break

    return 0, x, k, norm_dx


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    alpha = 0.5
    max_iter = 100
    eps = 1e-6

    line = input("alpha = ")
    if IsNumeric(line):
        alpha = float(line)

    # 1次元の場合
    x0 = 0.0
    ret, x, max_iter, eps = gradientdecent(DFunc2, x0, alpha, max_iter, eps)
    print("x = " + fmt(x))

    # # 多次元の場合
    # x0 = zeros(2, 1)
    # ret, x, max_iter, eps = gradientdecent_nd(DFunc4, x0, alpha, max_iter, eps)
    # print("x,y = " + fmt(x[0]) + ", " + fmt(x[1]))
    # print("f(x,y) = " + fmt(Func4(x)))

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
