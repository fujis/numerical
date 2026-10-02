# -*- coding: utf-8 -*-
"""!
@file monte-carlo.py

@brief 数値積分法
       モンテカルロ積分

@author Makoto Fujisawa
@date 2019-08 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math
import random
import time

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *
from rx_funcs import *


# -----------------------------------------------------------------------------
# 数値積分法
# -----------------------------------------------------------------------------
def monte_carlo(func, a, b, n):
    """!
    モンテカルロ法(1次元)
    @param[in] func 関数値を与える関数
    @param[in] a,b 積分範囲
    @param[in] n サンプル点個数
    @return 積分値
    """
    s = 0
    for i in range(n):
        rnd = random.random()  # [0,1]の乱数
        x = a+(b-a)*rnd
        f = func(x)
        s += f
    return (b-a)*s/n


# Xorshift用の内部状態(C++版では関数内のstatic変数)
_xorshift_y = int(time.time()) & 0xFFFFFFFF


def rand_xorshift():
    """!
    Xorshift(32bit)による乱数生成
     - https://ja.wikipedia.org/wiki/Xorshift のコード例そのまま
     - オリジナルは G. Marsaglia, "Xorshift RNGs". Journal of Statistical Software, 8(14), 2013.
     - Pythonの整数は桁あふれしないので，32bitに収まるよう毎回 0xFFFFFFFF でマスクする
    """
    global _xorshift_y
    y = _xorshift_y
    y = (y ^ (y << 13)) & 0xFFFFFFFF
    y = y ^ (y >> 17)
    y = (y ^ (y << 5)) & 0xFFFFFFFF
    _xorshift_y = y
    return y


def monte_carlo_2d(iofunc, x1, x2, y1, y2, n):
    """!
    モンテカルロ法(面積)
     - 関数f(x)=1として2重積分で面積を計算する方法
    @param[in] iofunc サンプル点が範囲内かどうかを判定する関数(範囲内で1,外で0を返す)
    @param[in] x1,x2,y1,y2 積分範囲
    @param[in] n 総サンプリング数
    @return 積分値
    """
    x = zeros(2)
    n_in = 0
    for i in range(n):
        rnd1 = random.random()  # [0,1]の乱数2
        rnd2 = random.random()  # [0,1]の乱数2
        # rnd1 = float(rand_xorshift() % 0xFFFFFF)/0xFFFFFF  # [0,1]の乱数1
        # rnd2 = float(rand_xorshift() % 0xFFFFFF)/0xFFFFFF  # [0,1]の乱数1
        x[0] = x1+(x2-x1)*rnd1
        x[1] = y1+(y2-y1)*rnd2
        if iofunc(x):
            n_in += 1  # 範囲内のサンプル点個数を数える
    V = (x2-x1)*(y2-y1)
    return V*(float(n_in)/float(n))


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
    n = 100

    # モンテカルロ積分
    s = monte_carlo(func, a, b, 100)
    print("n=" + str(n).rjust(4) + " : ", end="")
    print("monte carlo int f(x) = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
    print("ground truth = " + fmt(t))
    print("")

    # 円の面積の計算(nを変えて実行)
    r = sr
    t = RX_PI*r*r
    random.seed()  # 乱数のシード値を時間によって変える
    n = 10
    for i in range(7):
        s = monte_carlo_2d(FuncCircle, -r, r, -r, r, n)
        print("n=" + str(n).rjust(10) + " : ", end="")
        print("area of circle = " + fmt(abs(s)) + ",  error = " + fmt(abs(abs(s)-t)))
        n *= 10
    print("ground truth = " + fmt(t))

    # 誤差の平均値を求めるためのコード
    # m = 6
    # avg_error = zeros(7)
    # for j in range(10):
    #     random.seed()
    #     n = 10
    #     for i in range(m):
    #         s = monte_carlo_2d(FuncCircle, -r, r, -r, r, n)
    #         avg_error[i] += abs(abs(s)-t)
    #         n *= 10
    # n = 10
    # for i in range(m):
    #     avg_error[i] /= m
    #     print(str(n) + ", " + fmt(avg_error[i]))
    #     n *= 10

    return 0


if __name__ == "__main__":
    main()
