# -*- coding: utf-8 -*-
"""!
@file predictor-corrector.py

@brief アダムス・バッシュホース法とアダムス・ムルトン法を追加した予測子修正子法

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
def predictorcorrector_abm3(func, y0, a, b, n):
    """!
    アダムス・バッシュホース法(3点)+アダムスムルトン法(3点)による予測子修正子法
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
    f = [0, 0, 0, 0]  # f_(i-2), f_(i-1), f_(i), f_(i+1)
    for i in range(n):
        f[2] = func(x, y)
        if i < 2:  # 最初の方はオイラー法を使う
            # 前進オイラー法でy_(i+1)の予測値を計算
            f[3] = func(x+h, y+h*f[2])

            # 改良オイラー法で解を修正
            y = y+h*(f[2]+f[3])/2.0
        else:
            # アダムス・バッシュフォース法でy_(i+1)の予測値を計算
            y1 = y+h*(5*f[0]-16*f[1]+23*f[2])/12.0
            f[3] = func(x+h, y1)

            # アダムス・ムルトン法で解を修正
            y = y+h*(-f[1]+8*f[2]+5*f[3])/12.0

        f[0] = f[1]
        f[1] = f[2]
        x = x+h  # xの更新

        # 出力用
        # print("y(" + fmt(x) + ") = " + fmt(y))
        print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def predictorcorrector_abm4(func, y0, a, b, n):
    """!
    アダムス・バッシュホース法(4点)+アダムスムルトン法(4点)による予測子修正子法
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
    f = [0, 0, 0, 0, 0]  # f_(i-3), f_(i-2), f_(i-1), f_(i), f_(i+1)
    for i in range(n):
        f[3] = func(x, y)
        if i < 3:  # 最初の方はオイラー法を使う
            # 前進オイラー法でy_(i+1)の予測値を計算
            f[4] = func(x+h, y+h*f[3])

            # 改良オイラー法で解を修正
            y = y+h*(f[3]+f[4])/2.0
        else:
            # アダムス・バッシュフォース法でy_(i+1)の予測値を計算
            y1 = y+h*(-9*f[0]+37*f[1]-59*f[2]+55*f[3])/24.0
            f[4] = func(x+h, y1)

            # アダムス・ムルトン法で解を修正
            y = y+h*(f[1]-5*f[2]+19*f[3]+9*f[4])/24.0

        f[0] = f[1]
        f[1] = f[2]
        f[2] = f[3]
        x = x+h  # xの更新

        # 出力用
        # print("y(" + fmt(x) + ") = " + fmt(y))
        print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def predictorcorrector_fbe(func, y0, a, b, n):
    """!
    前進オイラー+改良オイラーによる予測子修正子法
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
    f = [0, 0]  # f_(i), f_(i+1)
    for i in range(n):
        f[0] = func(x, y)

        # 前進オイラー法でy_(i+1)の予測値を計算
        y1 = y+h*f[0]
        f[1] = func(x+h, y1)

        # 改良オイラー法で解を修正
        y = y+h*(f[0]+f[1])/2.0

        x = x+h  # xの更新

        # 出力用
        # print("y(" + fmt(x) + ") = " + fmt(y))
        print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

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

    # 予測子修正子法(4次精度)
    print("[Predictor-Corrector method with Adams-Bashforth(4) + Adams-Moulton(4)]")
    y = predictorcorrector_abm4(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # 予測子修正子法(3次精度)
    print("[Predictor-Corrector method with Adams-Bashforth(3) + Adams-Moulton(3)]")
    y = predictorcorrector_abm3(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # 予測子修正子法(2次精度)
    print("[Predictor-Corrector method with Forward Eular + Modified Eular]")
    y = predictorcorrector_fbe(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
