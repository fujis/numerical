# -*- coding: utf-8 -*-
"""!
@file romberg.py

@brief 数値積分法
       ロンバーグ法

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


def setw(x, w=10):
    """!
    数値を幅wで右詰めした文字列にする(C++の setw(w) に相当)
    @param[in] x 数値
    @param[in] w 幅
    @return 文字列
    """
    return fmt(x).rjust(w)


# -----------------------------------------------------------------------------
# 数値積分法
# -----------------------------------------------------------------------------
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


def romberg(func, a, b, k_max, eps):
    """!
    ロンバーグ法
     - 台形法の誤差を評価してより精度を高めた方法
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] k_max 最大分割回数(n=2^k_maxまで分割)
    @param[in] eps 許容誤差
    @return (積分値, 実際の分割回数, 実際の誤差)
    """
    h = b-a    # 分割幅(初期分割幅はn=1のもの)
    e = 0.0    # 誤差
    I = zeros(k_max+1)

    I[0] = (h/2)*(func(a)+func(b))  # I_0,0の計算

    # 計算の途中経過確認用出力
    print(setw(I[0]))

    n = 1
    l = 0
    k = 0
    for k in range(1, k_max+1):
        h = h/2
        n *= 2  # 分割幅を1/2にしていく
        I[k] = trapezoidal_integration(func, a, b, n)  # 台形公式でI_k,0を計算

        # 計算の途中経過確認用出力
        print(setw(I[k]), end="")

        # 収束判定
        e = abs(I[k]-I[k-1])
        if e < eps:
            l = k  # 結果が格納されている位置
            print("")
            break

        # 漸化式の計算
        m4 = 4  # 4^m
        m_broke = False
        for m in range(1, k+1):
            i = k-m  # I_k,mを格納する配列上の位置
            I[i] = (m4*I[i+1]-I[i])/(m4-1)  # I_k,mを計算
            m4 *= 4
            l = i  # 結果が格納されている位置

            # 計算の途中経過確認用出力
            print("   " + setw(I[i]), end="")

            # 収束判定2(ロンバーグ表で1つ上の値と比較)
            if i >= 1:
                e = abs(I[i]-I[i-1])
                if e < eps:
                    m_broke = True
                    break
        print("")
        if m_broke:
            break

    return I[l], k, e


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

    kmax = 5  # n=2^kmaxまで分割
    eps = 1.0e-6
    s, kmax, eps = romberg(func, a, b, kmax, eps)
    print("romberg int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    print("n = " + str(int(pow(2.0, float(kmax)))) + ", eps = " + fmt(eps))

    print("ground truth = " + fmt(t))

    # 円の面積の計算(上半分の積分-下半分の積分)
    kmax1, kmax2 = 7, 7  # n=2^kmaxまで分割
    eps1, eps2 = 1.0e-6, 1.0e-6
    r = sr
    a, b = -r, r
    t = RX_PI*r*r
    s1, kmax1, eps1 = romberg(FuncCircleTop, a, b, kmax1, eps1)
    s2, kmax2, eps2 = romberg(FuncCircleBottom, a, b, kmax2, eps2)
    s = s1-s2
    print("area of circle = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
