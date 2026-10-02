# -*- coding: utf-8 -*-
"""!
@file gauss.py

@brief 数値積分法
       ガウス型積分公式

@author Makoto Fujisawa
@date 2019-08 (Python版)
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
# 数値積分法
# -----------------------------------------------------------------------------
def gauss2(func, a, b):
    """!
    ガウス・ルジャンドル公式による積分計算(n=2)
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @return 積分値S
    """
    x = zeros(2)
    w = zeros(2)

    # 分点と重みの計算(n=2)
    x[0] = -math.sqrt(3.0)
    x[1] = -x[0]
    w[0] = w[1] = 1

    # Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
    S = 0.0
    for i in range(2):
        S += w[i]*func((b-a)*x[i]/2+(a+b)/2)
    S *= (b-a)/2

    return S


def gauss3(func, a, b):
    """!
    ガウス・ルジャンドル公式による積分計算(n=3)
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @return 積分値S
    """
    # 分点と重みの計算(n=3)
    x = zeros(3)
    w = zeros(3)
    x[0] = -math.sqrt(3.0/5.0)
    x[1] = 0
    x[2] = -x[0]
    w[0] = w[2] = 5.0/9.0
    w[1] = 8.0/9.0

    # Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
    S = 0.0
    for i in range(3):
        S += w[i]*func((b-a)*x[i]/2+(a+b)/2)
    S *= (b-a)/2

    return S


def gauss4(func, a, b):
    """!
    ガウス・ルジャンドル公式による積分計算(n=4)
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @return 積分値S
    """
    x = zeros(4)
    w = zeros(4)

    # 分点と重みの計算(n=4)
    tmp = 2.0*math.sqrt(6.0/5.0)
    x[0] = -math.sqrt((3.0+tmp)/7.0)
    x[1] = -math.sqrt((3.0-tmp)/7.0)
    x[2] = -x[1]
    x[3] = -x[0]
    w[0] = w[3] = (18.0-math.sqrt(30))/36.0
    w[1] = w[2] = (18.0+math.sqrt(30))/36.0

    # Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
    S = 0.0
    for i in range(4):
        S += w[i]*func((b-a)*x[i]/2+(a+b)/2)
    S *= (b-a)/2

    return S


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 指数関数の[0,1]での積分 (解析解はe-1=1.718281828459045235360287471352...)
    # func = FuncExp
    # a, b = 0.0, 1.0
    # t = math.exp(1.0)-1  # 真値

    # 講義で示した例題(2017年度筑波大学前期日程入試問題)
    func = FuncT17
    a, b = 0.5, 2.0
    t = -99.0/8.0+18.0*math.log(2.0)  # 真値

    set_precision(10)

    print("")
    s = gauss2(func, a, b)
    print("gauss2 int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    s = gauss3(func, a, b)
    print("gauss3 int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    s = gauss4(func, a, b)
    print("gauss4 int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
