# -*- coding: utf-8 -*-
"""!
@file lu.py

@brief LU分解

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
# 三角分解
# -----------------------------------------------------------------------------
def LUDecomp(A, n):
    """!
    LU分解(ピボット交換なし)
     - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
     - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
     - LとUを一つの行列にまとめた形で結果を返す
    @param[inout] A n×nの係数行列．LU分解した結果を格納する．
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    for i in range(n):
        # l_ijの計算(i >= j)
        for j in range(i+1):
            lu = A[i][j]
            for k in range(j):
                lu -= A[i][k]*A[k][j]  # l_ik * u_kj
            A[i][j] = lu

        # u_ijの計算(i < j)
        for j in range(i+1, n):
            lu = A[i][j]
            for k in range(i):
                lu -= A[i][k]*A[k][j]  # l_ik * u_kj
            A[i][j] = lu/A[i][i]

    return 0


def LUSolver(A, b, x, n):
    """!
    LU分解した行列A(n×n)から前進代入・後退代入によりA・x=bを解く
    @param[in] A LU分解された行列
    @param[in] b 右辺ベクトル
    @param[out] x 結果ベクトル
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    # 前進代入(forward substitution)
    #  LY=bからYを計算
    for i in range(n):
        bly = b[i]
        for j in range(i):
            bly -= A[i][j]*x[j]
        x[i] = bly/A[i][i]

    # 後退代入(back substitution)
    #  UX=YからXを計算
    for i in range(n-1, -1, -1):
        yux = x[i]
        for j in range(i+1, n):
            yux -= A[i][j]*x[j]
        x[i] = yux

    return 0


def LUInverse(mat, inv):
    """!
    LU分解で逆行列を計算
    @param[in]  mat 計算したい行列
    @param[out] inv 逆行列
    @return 0:成功,1:失敗(n <= 0)
    """
    n = len(mat)
    if n <= 0:
        return 1

    # 行列を1回だけLU分解
    LUDecomp(mat, n)

    b = zeros(n)
    x = zeros(n)

    # 行ごとに逆行列を算出
    for j in range(n):
        for i in range(n):
            b[i] = 1.0 if i == j else 0.0

        LUSolver(mat, b, x, n)
        for i in range(n):
            inv[i][j] = x[i]

    return 0


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込んで解く
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    Ab = ReadMatrix(fn, ",")
    n = len(Ab)  # n元連立一次方程式

    # 読み込んだ行列を確認用に画面表示
    print("A(" + str(n) + " x " + str(n+1) + ") = ")
    OutputMatrix(Ab, n, n+1)
    print("")

    # 左辺の行列の抽出
    A = zeros2(n, n)
    for i in range(n):
        A[i] = Ab[i][0:n]

    # 右辺項の抽出
    b = zeros(n)
    for i in range(n):
        b[i] = Ab[i][n]

    # LU分解
    LUDecomp(A, n)

    # LU分解結果を画面表示(L行列とU行列をひとつにまとめた行列として表示)
    print("LU = ")
    OutputMatrix(A, n, n)
    print("")

    # LU分解を用いて連立1次方程式を解く
    x = zeros(n)
    LUSolver(A, b, x, n)

    # 結果の画面表示
    s = ""
    for i in range(n):
        s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
    print(s)

    # LU分解を用いて逆行列を算出する
    invA = zeros2(n, n)
    for i in range(n):
        A[i] = Ab[i][0:n]
    LUInverse(A, invA)

    # 逆行列を画面表示
    print("inv(A) = ")
    OutputMatrix(invA, n, n)
    print("")

    # 逆行列チェック
    C = zeros2(n, n)
    MulMatrix(Ab, invA, C, n)

    # 元の行列と計算した逆行列を掛けた結果の行列を画面表示(これがほぼ単位行列ならOK)
    print("I = ")
    OutputMatrix(C, n, n)
    print("")

    return 0


if __name__ == "__main__":
    main()
