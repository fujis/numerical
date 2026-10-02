# -*- coding: utf-8 -*-
"""!
@file goldensection.py

@brief 黄金分割法

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
def goldensection(func, xl, xr, max_iter, eps):
    """!
    黄金分割探索法(golden-section method)
    @param[in] func 関数値を与える関数
    @param[in] xl,xr 初期探索範囲
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解ans, 実際の反復数, 実際の誤差)
    """
    e = (SQRT5-1.0)/(SQRT5+1.0)
    # e = 2.0/(SQRT5+3.0)  # 上の式の変形(√5を事前計算してないならこちらの方が速い)
    b = xr-xl
    d = e*b
    set_precision(7)

    x = zeros(4)
    x[0] = xl
    x[1] = xl+d
    x[2] = xr-d
    x[3] = xr

    # x[1],x[2]における関数値
    f1 = func(x[1])
    f2 = func(x[2])

    k = 0
    for k in range(max_iter):
        print(str(k) + " : " + fmt(x[0]) + " " + fmt(x[1]) + " " + fmt(x[2]) + " " + fmt(x[3]) + ", b = " + fmt(b))

        if f1 < f2:  # 区間[x2,x3]に最小値なし
            x2 = x[2]
            x[2] = x[1]
            x[1] = x2-d
            x[3] = x2
            f2 = f1
            f1 = func(x[1])
        else:  # 区間[x0,x1]に最小値なし
            x1 = x[1]
            x[1] = x[2]
            x[2] = x1+d
            x[0] = x1
            f1 = f2
            f2 = func(x[2])

        b -= d   # 新しい[x0,x3]の長さ(現在の長さから[x0,x1](or[x2,x3])の長さdを引く)
        d = e*b  # 新しい[x0,x1](or[x2,x3])の長さ

        if b < eps:
            break

    ans = 0.5*(x[1]+x[2])

    return 0, ans, k, b


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    x1 = 0.0
    x2 = 2.0

    max_iter = 100
    eps = 1e-6
    ret, x, max_iter, eps = goldensection(Func2, x1, x2, max_iter, eps)

    print("x = " + fmt(x))
    print("f(x) = " + fmt(Func2(x)))

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
