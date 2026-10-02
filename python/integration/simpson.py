# -*- coding: utf-8 -*-
"""!
@file simpson.py

@brief 数値積分法
       シンプソン公式

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
def simpson_integration(func, a, b, n):
    """!
    シンプソン法(シンプソンの1/3公式)
     - 台形法の誤差を評価してより精度を高めた方法
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] n 積分区間の分割数(2区間で1つの多項式なのでここでの分割数はデータ点数/2)
    @return 積分値
    """
    h = (b-a)/(2*n)  # 分割区間の横幅
    # f : 分割区間の縦幅

    S = func(a)+func(b)

    # 奇数項
    for i in range(1, n+1):
        f = func(a+(2*i-1)*h)
        S += 4*f
    # 偶数項
    for i in range(1, n):
        f = func(a+(2*i)*h)
        S += 2*f
    S *= h/3

    return S


def simpson38_integration(func, a, b, n):
    """!
    シンプソン法(シンプソンの3/8公式)
     - 台形法の誤差を評価してより精度を高めた方法
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] n 積分区間の分割数(3区間で1つの多項式なのでここでの分割数はデータ点数/3)
    @return 積分値
    """
    h = (b-a)/(3*n)  # 分割区間の横幅
    # f : 分割区間の縦幅

    S = func(a)+func(b)

    # 3の倍数-2の項
    for i in range(1, n+1):
        f = func(a+(3*i-2)*h)
        S += 3*f
    # 3の倍数-1の項
    for i in range(1, n+1):
        f = func(a+(3*i-1)*h)
        S += 3*f
    # 3の倍数の項
    for i in range(1, n):
        f = func(a+(3*i)*h)
        S += 2*f
    S *= 3*h/8

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
    n = 120  # シンプソン公式のためにnは2と3の公倍数にしておくこと

    s = simpson_integration(func, a, b, n//2)
    print("simpson int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))

    s = simpson38_integration(func, a, b, n//3)
    print("simpson(3/8) int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))

    print("ground truth = " + fmt(t))

    # 円の面積の計算(上半分の積分-下半分の積分)
    n = 120
    r = 1.0
    a, b = -r, r
    t = RX_PI*r*r
    s = simpson_integration(FuncCircleTop, a, b, n)-simpson_integration(FuncCircleBottom, a, b, n)
    print("area of circle = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
