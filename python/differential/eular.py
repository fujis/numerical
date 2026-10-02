# -*- coding: utf-8 -*-
"""!
@file eular.py

@brief オイラー法，ホイン法(改良オイラー法)

@author Makoto Fujisawa
@date 2019-10 (Python版)
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
def eular(func, y0, a, b, n):
    """!
    オイラー法(1次精度)
     - y(n+1)=y(n)+hf(x(n),y(n))
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
        fi = func(x, y)
        y = y+h*fi  # yの更新
        x = x+h     # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def heun(func, y0, a, b, n):
    """!
    ホイン法(2次精度)
     - y(n+1)=y(n)+h/2(f(x(n),y(n))+f(x(n+1),y(n+1)))
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
        # Y(i+1) = f(x(i+1)+y(i+1))をオイラー法で求める
        fi = func(x, y)
        Yi = y+h*fi  # Yiの計算

        # f(x,y)とY(i+1)の平均でyを更新
        y = y+h*(fi+func(x+h, Yi))/2.0
        x = x+h      # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def backward_eular(func, dfunc, y0, a, b, n, max_iter, eps):
    """!
    後退オイラー法(1次精度)
     - y(n+1)=y(n)+hf(x(n+1),y(n+1))
     - ニュートン法を使って計算
    @param[in] func 関数f(x,y)の値を与える関数
    @param[in] dfunc 関数f(x,y)のyについての偏微分を与える関数
    @param[in] y0 初期値(y0=g(a))
    @param[in] a,b 計算範囲
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @param[in] max_iter ニュートン法の最大反復数
    @param[in] eps ニュートン法の許容誤差
    @return x=bでの解
    """
    k_avg = 0
    h = (b-a)/n  # 刻み幅
    x = a   # xの初期値
    y = y0  # yの初期値
    # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力
    for i in range(n):
        yk = y  # y_(i+1)計算用の一時変数

        # ニュートン法でy_(i+1)を求める
        k = 0
        for k in range(max_iter):
            g = yk-h*func(x+h, yk)-y
            dg = 1-h*dfunc(x+h, yk)
            yk = yk-g/dg
            if abs(g/dg) < eps or abs(g) < eps:
                break
        k_avg += k

        y = yk  # yの更新
        # y = y/(1+25*h)
        x = x+h     # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力
    print("average number of iterations : " + fmt(k_avg/float(n)))

    return y


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    global TF

    # 常微分方程式dy/dx=-λy (解析解はy = C e^-λx = y0 e^-λx)
    func = FuncOdeY
    dfunc = DyFuncOdeY
    a, b = 0.0, 1.0  # 範囲[a,b]
    y0 = 1.0         # 初期値
    # λの値を変更(Pythonではモジュール側の変数を直接書き換える必要がある)
    rx_funcs.lmbda = 25
    TF = lambda x: FuncOdeY_true(x, y0)  # 真値

    # # 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
    # func = FuncOdeXY
    # dfunc = DyFuncOdeXY
    # a, b = 0.0, 1.0  # 範囲[a,b]
    # y0 = 1.0         # 初期値
    # TF = lambda x: FuncOdeXY_true(x, y0)  # 真値

    set_precision(10)
    n = 10
    t = TF(b)  # 真値

    # オイラー法(1次精度)
    print("[Eular method]")
    y = eular(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # ホイン法(2次精度)
    print("[Heun method]")
    y = heun(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # 後退オイラー法(1次精度)
    print("[backward Eular method]")
    y = backward_eular(func, dfunc, y0, a, b, n, 30, 1e-6)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
