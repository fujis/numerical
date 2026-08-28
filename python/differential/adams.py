# -*- coding: utf-8 -*-
"""!
@file adams.py

@brief アダムス・バッシュホース法

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
def adamsbashforth3(func, y0, a, b, n):
    """!
    アダムス・バッシュホース法(3点)
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
    fi1 = 0
    fi2 = 0  # f_i, f_(i-1), f_(i-2)
    for i in range(n):
        fi = func(x, y)
        if i < 2:  # 最初の方はオイラー法を使う
            y = y+h*fi  # yの更新
        else:
            y = y+h*(5*fi2-16*fi1+23*fi)/12.0  # yの更新

        fi2 = fi1
        fi1 = fi
        x = x+h  # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


def adamsbashforth4(func, y0, a, b, n):
    """!
    アダムス・バッシュホース法(4点)
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
    fi1 = 0
    fi2 = 0
    fi3 = 0  # f_i, f_(i-1), f_(i-2), f_(i-3)
    for i in range(n):
        fi = func(x, y)
        if i < 3:  # 最初の方はオイラー法を使う
            y = y+h*fi  # yの更新
        else:
            y = y+h*(-9*fi3+37*fi2-59*fi1+55*fi)/24.0  # yの更新

        fi3 = fi2
        fi2 = fi1
        fi1 = fi
        x = x+h  # xの更新

        # 出力用
        print("y(" + fmt(x) + ") = " + fmt(y))
        # print(fmt(x) + ", " + fmt(y) + ", " + fmt(TF(x)) + ", " + fmt(abs(y-TF(x))))  # グラフ描画用に真値も出力

    return y


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    global TF

    # # 常微分方程式dy/dx=-λy (解析解はy = C e^-λx = y0 e^-λx)
    # func = FuncOdeY
    # a, b = 0.0, 1.0  # 範囲[a,b]
    # y0 = 1.0         # 初期値
    # rx_funcs.lmbda = 25
    # TF = lambda x: FuncOdeY_true(x, y0)  # 真値

    # 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
    func = FuncOdeXY
    a, b = 0.0, 1.0  # 範囲[a,b]
    y0 = 1.0         # 初期値
    TF = lambda x: FuncOdeXY_true(x, y0)  # 真値

    set_precision(10)
    n = 10
    t = TF(b)  # 真値

    # アダムス・バッシュホース法(3点,3次精度)
    print("[Adams-Bashforth(3)]")
    y = adamsbashforth3(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    # アダムス・バッシュホース法(4点,4次精度)
    print("[Adams-Bashforth(4)]")
    y = adamsbashforth4(func, y0, a, b, n)
    print("y(" + fmt(b) + ") = " + fmt(y) + ",  error = " + fmt(abs(y-t)))
    print("")

    print("ground truth = " + fmt(t))

    return 0


if __name__ == "__main__":
    main()
