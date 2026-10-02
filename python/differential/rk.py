# -*- coding: utf-8 -*-
"""!
@file rk.py

@brief ルンゲ・クッタ法

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
import rx_funcs  # 方程式の係数(λなど)を書き換えるためにモジュールとしても読み込む
from rx_funcs import *


# -----------------------------------------------------------------------------
# デバッグ用変数
# -----------------------------------------------------------------------------
TF = None  # 真値を与える関数


# -----------------------------------------------------------------------------
# 常微分方程式の近似解
# -----------------------------------------------------------------------------
def rk3(func, y0, a, b, n):
    """!
    3段3次のルンゲ・クッタ法(クッタの3次公式)(3次精度)
     - y(n+1)=y(n)+(h/6)(k1+4*k2+k3)
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @return x=bでの解
    """
    h = (b-a)/n  # 刻み幅

    x = a   # xの初期値
    y = y0  # yの初期値
    print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力
    for i in range(n):
        k1 = func(x, y)                 # k1の算出
        k2 = func(x+h/2, y+(h/2)*k1)    # k2の算出
        k3 = func(x+h, y-h*k1+2*h*k2)   # k3の算出

        y = y+(h/6)*(k1+4*k2+k3)  # yの更新
        x = x+h                   # xの更新

        # 出力用
        # print("y(" + fmt(x) + ") = " + fmt(y))
        print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def rk4(func, y0, a, b, n):
    """!
    4段4次のルンゲ・クッタ法(4次精度)
     - y(n+1)=y(n)+(h/6)(k1+2*k2+2*k3+k4)
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @return x=bでの解
    """
    h = (b-a)/n  # 刻み幅

    x = a   # xの初期値
    y = y0  # yの初期値
    # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力
    for i in range(n):
        k1 = func(x, y)                 # k1の算出
        k2 = func(x+h/2, y+(h/2)*k1)    # k2の算出
        k3 = func(x+h/2, y+(h/2)*k2)    # k3の算出
        k4 = func(x+h, y+h*k3)          # k4の算出

        y = y+(h/6)*(k1+2*k2+2*k3+k4)  # yの更新
        x = x+h                        # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    global TF

    # 常微分方程式dy/dx=-λy (解析解はy = C e^-λx = y0 e^-λx)
    func = FuncOdeY
    a, b = 0.0, 1.0  # 範囲[a,b]
    y0 = 1.0         # 初期値
    rx_funcs.lmbda = 25
    TF = lambda x: FuncOdeY_true(x, y0)  # 真値

    # # 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
    # func = FuncOdeXY
    # a, b = 0.0, 1.0  # 範囲[a,b]
    # y0 = 1.0         # 初期値
    # TF = lambda x: FuncOdeXY_true(x, y0)  # 真値

    set_precision(10)
    n = 20
    t = TF(b)  # 真値

    # 3段3次のルンゲ・クッタ法(3次精度)
    print("[RK3]")
    y = rk3(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # 4段4次のルンゲ・クッタ法(4次精度)
    print("[RK4]")
    y = rk4(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
