# -*- coding: utf-8 -*-
"""!
@file integration.py

@brief 数値積分法
       区分求積法と台形公式

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
def segment_integration(func, a, b, n):
    """!
    区分求積法
     - 積分区間を矩形で分割する方法
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] n 積分区間の分割数
    @return 積分値
    """
    h = (b-a)/n  # 分割区間の横幅

    S = 0
    for i in range(n):
        f = func(a+(i+1)*h)  # xの大きい方の辺の長さを縦幅とする
        S += f*h

    return S


def trapezoidal_integration(func, a, b, n):
    """!
    台形法(台形公式)
     - 積分区間を台形で分割する方法
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] n 積分区間の分割数
    @return 積分値
    """
    h = (b-a)/n  # 分割区間の横幅
    # f1, f2 : 分割区間の縦幅(台形の長辺と短辺)

    f2 = func(a)
    S = 0
    for i in range(n):
        f1 = f2
        f2 = func(a+(i+1)*h)
        S += (f1+f2)*h/2

    return S


def trapezoidal_integration2(func, a, b, y1f, y2f, n, m):
    """!
    2重積分(台形公式)
    @param[in] func 関数値を与える関数
    @param[in] a,b x方向積分範囲
    @param[in] y1f,y2f y方向積分範囲を与える関数
    @param[in] n,m x,y方向の積分区間の分割数
    @return 積分値
    """
    h1 = (b-a)/n  # x方向刻み幅
    S = 0.0
    for i in range(n+1):
        xi = a+i*h1
        h2 = (y2f(xi)-y1f(xi))/m  # y方向刻み幅
        y1 = y1f(xi)
        if abs(h2) < 1e-10:
            continue

        # 台形公式によるy方向積分(Fiの計算)
        f2 = func(xi, y1)
        Fi = 0.0
        for j in range(n):
            f1 = f2
            f2 = func(xi, y1+(j+1)*h2)
            Fi += (f1+f2)*h2/2

        if i == 0 or i == n:
            S += Fi
        else:
            S += 2*Fi
    S *= h1/2

    return S


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 指数関数の[0,1]での積分 (解析解はe-1=1.718281828459045235360287471352...)
    func = FuncExp
    a, b = 0.0, 1.0
    t = math.exp(1.0)-1  # 真値

    # 講義で示した例題(2017年度筑波大学前期日程入試問題)
    # func = FuncT17
    # a, b = 0.5, 2.0
    # t = -99.0/8.0+18.0*math.log(2.0)  # 真値

    set_precision(10)
    n = 10

    s = segment_integration(func, a, b, n)
    print("segment int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))

    s = trapezoidal_integration(func, a, b, n)
    print("trapezoidal int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))

    print("ground truth = " + fmt(t))

    print("")

    # 円の面積の計算(上半分の積分-下半分の積分)
    # n = 100
    # r = 1.0
    # a, b = -r, r
    # t = RX_PI*r*r
    # s = trapezoidal_integration(FuncCircleTop, a, b, n)-trapezoidal_integration(FuncCircleBottom, a, b, n)
    # print("area of circle = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    # print("ground truth = " + fmt(t))

    # 重積分
    # a, b = 1, 2
    # t = 54  # 真値
    # s = trapezoidal_integration2(FuncP2, a, b, FuncY1, FuncY2, n, n)
    # print("trapezoidal int2 f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    # print("ground truth = " + fmt(t))

    # # 球の体積(こちらは上半分の積分x2で計算)
    # a, b = -sr, sr
    # t = (4*RX_PI*sr*sr*sr)/3.0  # 真値
    # s = 2*trapezoidal_integration2(FuncSphere, a, b, FuncSphereY1, FuncSphereY2, n, n)
    # print("trapezoidal int2 f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    # print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
