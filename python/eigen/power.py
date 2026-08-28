# -*- coding: utf-8 -*-
"""!
@file power.py

@brief 固有値・固有ベクトル
       べき乗法(power method)

@author Makoto Fujisawa
@date 2019-12 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *


def eigen_power(A, x, n, max_iter, eps):
    """!
    べき乗法による固有値の算出
    @param[in] A 正方行列(n×n)
    @param[in] x 絶対値最大固有値に対応する固有ベクトルの初期値(大きさnのリスト)
    @param[in] n 行列のサイズ(n×n)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 絶対値最大固有値λ, 固有ベクトルx, 実際の反復数, 実際の誤差)
    """
    x1 = zeros(n)  # x^(k+1)格納用
    l2 = dot(x, x, n)
    e = 0.0
    lambda_ = 0.0
    k = 0
    for k in range(max_iter):
        # |x^(k)|=1となるように正規化
        l = math.sqrt(l2)         # |x^(k)|の計算
        x = mul_sv(1.0/l, x, n)   # |x^(k)|で割る(1/|x|を掛ける)

        # x^(k+1) = A x^(k)の計算
        x1 = mul_mv(A, x, n)

        # 固有値λの計算(|x^(k)|=1で正規化されていることが前提)
        lambda_ = dot(x, x1, n)

        print(str(k) + " : lambda = " + fmt(lambda_) + ",  v = " + fmt(x))

        # 収束判定
        l2 = dot(x1, x1, n)
        e = abs(l2-lambda_*lambda_)
        if e < eps*eps:
            break

        x = x1

    x = x1
    x = mul_sv(1.0/math.sqrt(l2), x, n)  # 最後に固有ベクトルを正規化しておく

    return 1, lambda_, x, k, math.sqrt(e)


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込んで解く
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    A = ReadMatrix(fn, ",")

    n = len(A)  # n×n行列
    print("A(" + str(n) + " x " + str(n) + ") = ")
    OutputMatrix(A, n, n)
    print("")

    v = zeros(n, 1.0)  # vを1で初期化
    set_precision(8)

    max_iter = 100
    eps = 1e-6
    ret, lambda_, v, max_iter, eps = eigen_power(A, v, n, max_iter, eps)

    # 絶対値最大固有値と固有ベクトルの表示
    print("lambda_max = " + fmt(lambda_))
    print("v = " + fmt(v))

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
