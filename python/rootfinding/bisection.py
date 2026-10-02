# -*- coding: utf-8 -*-
"""!
@file bisection.py

@brief 二分法による求根問題の解法

@author Makoto Fujisawa
@date 2012-06,2019-09 (Python版)
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
# 二分法による求根問題の解法
# -----------------------------------------------------------------------------
def bisection(func, xl, xr, max_iter, eps):
    """!
    2分法(bisection method)
     - C++版では解x,反復数,誤差を参照渡しで返していたが，
       Pythonでは複数の値をまとめて戻り値で返す
    @param[in] func 関数値を与える関数
    @param[in] xl,xr 探索範囲
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解x, 実際の反復数, 実際の誤差)
    """
    f = func(xl)
    fmid = func(xr)

    # 探索範囲の境界での値の符号が異なる場合のみ探索
    if f*fmid >= 0.0:
        return 0, 0.0, max_iter, eps

    dx = abs(xr-xl)
    xmid = 0.0
    k = 0
    for k in range(max_iter):
        xmid = 0.5*(xl+xr)  # 中点
        dx *= 0.5

        # 中点での関数値を求める
        fmid = func(xmid)

        # 確認用の画面出力
        print(str(k) + " : [" + fmt(xl) + ", " + fmt(xr) + "], fmid = " + fmt(fmid))

        # 収束判定
        if dx < eps or fmid == 0.0:
            break

        # 新しい区間
        if f*fmid < 0:
            xr = xmid
        else:
            xl = xmid
            f = fmid

    return 0, xmid, k, dx


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 探索範囲
    x1 = -1.0
    x2 = 1.0

    # 二分法でf(x)=0を解く
    max_iter = 100
    eps = 1e-6
    ret, x, max_iter, eps = bisection(Func1, x1, x2, max_iter, eps)

    set_precision(15)
    # 結果の画面表示
    print("x = " + fmt(x))
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    # print(fmt(math.log2(2.0/1.0e-6)))

    return 0


if __name__ == "__main__":
    main()
