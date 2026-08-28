# -*- coding: utf-8 -*-
"""!
@file jacobi.py

@brief ヤコビ(Jacobi)反復法

@author Makoto Fujisawa
@date 2012-06,2019-09 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *


# -----------------------------------------------------------------------------
# 反復法による線形システムソルバ
# -----------------------------------------------------------------------------
def JacobiIteration(A, n, max_iter, eps):
    """!
    ヤコビ反復法(Jacobi iterative method)
     - 解が収束するのは
         ・対角有利(diagonal dominant, 対角要素の絶対値>その行の他の要素の絶対値の和)
         ・係数行列が対称(symmetric)かつ正定(positive definite)
       のどちらかの場合
    @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の拡大行列．n+1列目に解が入る．
    @param[in] n n元連立一次方程式
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    x = zeros(n)  # 初期値はすべて0とする
    y = zeros(n)

    e = 0.0  # 誤差平均
    k = 0
    for k in range(max_iter):
        e = 0.0
        # 現在の値を代入して，次の解候補を計算
        for i in range(n):
            y[i] = A[i][n]
            for j in range(n):
                y[i] -= (A[i][j]*x[j] if j != i else 0.0)
            y[i] /= A[i][i]

            e += abs(y[i]-x[i])

        # 確認のため現在の解を画面出力
        set_precision(8)
        s = str(k) + " : "
        for i in range(n):
            s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
        # for i in range(n): s += ", " + fmt(abs(y[i]-x[i]))
        # s += ", " + fmt(e/n)
        print(s)

        # 収束判定
        converged = True
        for l in range(n):
            # if abs((y[l]-x[l])/y[l]) > eps:  # 相対誤差の場合
            if abs(y[l]-x[l]) > eps:  # 絶対誤差の場合
                converged = False
                break
        if converged:  # すべての解が許容誤差以下なら反復終了
            break

        x, y = y, x

    # 解をA_(i,n)に格納
    for i in range(n):
        A[i][n] = y[i]

    return 0, k, e/n


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込んで解く
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    A = ReadMatrix(fn, ",")

    n = len(A)     # n元連立一次方程式
    m = len(A[0])

    # 読み込んだ行列を確認用に画面表示
    print("A(" + str(n) + " x " + str(m) + ") = ")
    OutputMatrix(A, n, n+1)
    print("")

    # ヤコビ反復法で線形システムを解く
    max_iter = 100
    eps = 1e-6
    ret, max_iter, eps = JacobiIteration(A, n, max_iter, eps)

    # 結果の画面表示
    s = ""
    for i in range(n):
        s += "x" + str(i) + " = " + fmt(A[i][n]) + ("" if i == n-1 else ", ")
    print(s)

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
