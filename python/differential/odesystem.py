# -*- coding: utf-8 -*-
"""!
@file odesystem.py

@brief 連立常微分方程式(system of ordinary differential equations)

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
import rx_funcs  # 方程式の係数を書き換えるためにモジュールとしても読み込む
from rx_funcs import *


# -----------------------------------------------------------------------------
# デバッグ用変数
# -----------------------------------------------------------------------------
TF = None  # 真値を与える関数


# -----------------------------------------------------------------------------
# 常微分方程式の近似解
# -----------------------------------------------------------------------------
def eular(func, y0, a, b, n):
    """!
    前進オイラー法(多変数版)
     - y(n+1)=y(n)+hf(x(n),y(n))
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @return x=bでの解
    """
    h = (b-a)/n  # 刻み幅

    x = a  # xの初期値
    y = list(y0)
    # print(fmt(x) + ", " + fmt(y))  # グラフ描画用
    for k in range(n):
        f = func(x, y)
        for i in range(len(y0)):
            y[i] = y[i]+h*f[i]  # yの更新
        x = x+h  # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y))  # グラフ描画用

    return y


def heun(func, y0, a, b, n):
    """!
    ホイン法(多変数版)
     - y(n+1)=y(n)+h/2(f(x(n),y(n))+f(x(n+1),y(n+1)))
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @return x=bでの解
    """
    h = (b-a)/n  # 刻み幅

    x = a  # xの初期値
    y = list(y0)
    Yi = list(y0)
    # print(fmt(x) + ", " + fmt(y))  # グラフ描画用
    for k in range(n):
        f = func(x, y)
        for i in range(len(y0)):
            Yi[i] = y[i]+h*f[i]  # Yiの計算

        # f(x,y)とY(i+1)の平均でyを更新
        f1 = func(x+h, Yi)
        for i in range(len(y0)):
            y[i] = y[i]+h*(f[i]+f1[i])/2.0

        x = x+h  # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y))  # グラフ描画用

    return y


def rk4(func, y0, a, b, n):
    """!
    4段4次のルンゲ・クッタ法(多変数版)
     - y(n+1)=y(n)+(h/6)(k1+2*k2+2*k3+k4)
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @return x=bでの解
    """
    h = (b-a)/n  # 刻み幅

    x = a  # xの初期値
    y = list(y0)
    yk = list(y0)
    # print(fmt(x) + ", " + fmt(y))  # グラフ描画用
    for k in range(n):
        k1 = func(x, y)

        for i in range(len(y0)):
            yk[i] = y[i]+(h/2)*k1[i]
        k2 = func(x+h/2, yk)

        for i in range(len(y0)):
            yk[i] = y[i]+(h/2)*k2[i]
        k3 = func(x+h/2, yk)

        for i in range(len(y0)):
            yk[i] = y[i]+h*k3[i]
        k4 = func(x+h, yk)

        # yの更新
        for i in range(len(y0)):
            y[i] = y[i]+(h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])

        x = x+h  # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y))  # グラフ描画用

    return y


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # # Lotka-Volterraモデル
    # func = FuncOdeLV
    # a, b = 0.0, 1000.0  # 範囲[a,b] → 時刻(日数)
    # n = 1000
    # y0 = zeros(2)  # 初期値
    # y0[0] = 300; y0[1] = 300

    # 単振り子モデル
    func = FuncOdePendulum
    a, b = 0.0, 10.0  # 範囲[a,b] → 時刻(秒数)
    n = 100
    y0 = zeros(2)  # 初期値
    # y0[0] = RX_PI/4.0; y0[1] = 0
    y0[0] = 5.0
    y0[1] = 0

    set_precision(10)

    # 多変数版オイラー法(1次精度)
    # print("[Eular method]")
    # y = eular(func, y0, a, b, n)
    # print("")

    # 多変数版ホイン法(2次精度)
    # print("[Heun method]")
    # y = heun(func, y0, a, b, n)
    # print("")

    # 多変数版RK4(4次精度)
    print("[RK4]")
    y = rk4(func, y0, a, b, n)
    print("")

    return 0


if __name__ == "__main__":
    main()
