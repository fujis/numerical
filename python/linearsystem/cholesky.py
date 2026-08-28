# -*- coding: utf-8 -*-
"""!
@file cholesky.py

@brief 修正/不完全コレスキー分解

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
# コレスキー分解
# -----------------------------------------------------------------------------
def CholeskyDecomp(A, L, n):
    """!
    コレスキー分解
     - 正定値対称行列A(n×n)を下三角行列(L:Lower triangular matrix)とその転置(L^T)に分解する
     - L: i >= j (i>jとi==jに分けて処理)
    @param[in] A n×nの対称行列
    @param[out] L 下三角行列
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    for j in range(n):
        # i == jについて解く
        ll = A[j][j]
        for k in range(j):
            ll -= L[j][k]*L[j][k]
        L[j][j] = math.sqrt(ll)

        for i in range(j+1, n):
            ll = A[i][j]
            for k in range(j):
                ll -= L[i][k]*L[j][k]
            L[i][j] = ll/L[j][j]

    return 0


def ModifiedCholeskyDecomp(A, L, d, n):
    """!
    修正コレスキー分解(modified Cholesky decomposition)
     - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
     - l_ii = 1とした場合
     - L: i > jの要素が非ゼロで対角成分は1
    @param[in] A n×nの対称行列
    @param[out] L 対角成分が1の下三角行列
    @param[out] d 対角行列(対角成分のみ)
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    d[0] = A[0][0]
    L[0][0] = 1.0

    for i in range(1, n):
        # i < k の場合
        for j in range(i):
            lld = A[i][j]
            for k in range(j):
                lld -= L[i][k]*L[j][k]*d[k]
            L[i][j] = (1.0/d[j])*lld

        # i == k の場合
        ld = A[i][i]
        for k in range(i):
            ld -= L[i][k]*L[i][k]*d[k]
        d[i] = ld
        L[i][i] = 1.0

    return 0


def IncompleteCholeskyDecomp(A, L, d, n):
    """!
    不完全コレスキー分解(incomplete Cholesky decomposition)
     - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
     - l_ii = 1とした場合
     - L: i > jの要素が非ゼロで対角成分は1
     - 行列Aの値が0である要素に対応する部分を飛ばす
    @param[in] A n×nの対称行列
    @param[out] L 対角成分が1の下三角行列
    @param[out] d 対角行列(対角成分のみ)
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    d[0] = A[0][0]
    L[0][0] = 1.0

    for i in range(1, n):
        # i < k の場合
        for j in range(i):
            if abs(A[i][j]) < 1.0e-10:
                continue

            lld = A[i][j]
            for k in range(j):
                lld -= L[i][k]*L[j][k]*d[k]
            L[i][j] = (1.0/d[j])*lld

        # i == k の場合
        ld = A[i][i]
        for k in range(i):
            ld -= L[i][k]*L[i][k]*d[k]
        d[i] = ld
        L[i][i] = 1.0

    return 0


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込んで解く
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    Ab = ReadMatrix(fn, ",")
    n = len(Ab)  # n元連立一次方程式

    # 左辺の行列の抽出
    A = zeros2(n, n)
    for i in range(n):
        A[i] = Ab[i][0:n]

    # 読み込んだ行列を確認用に画面表示
    print("A(" + str(n) + " x " + str(n+1) + ") = ")
    OutputMatrix(A, n, n)
    print("")

    # コレスキー分解
    L = zeros2(n, n)
    d = []
    CholeskyDecomp(A, L, n)

    # 修正コレスキー分解，不完全コレスキー分解
    # d = zeros(n)
    # ModifiedCholeskyDecomp(A, L, d, n)
    # IncompleteCholeskyDecomp(A, L, d, n)

    print("L = ")
    OutputMatrix(L, n, n)

    if not d:  # コレスキー分解のチェック
        # 分解結果のチェック用
        LL = zeros2(n, n)
        MulMatrix(L, transpose(L, n), LL, n)

        # 分解した結果を掛け合わせた行列を画面表示(これが元の行列Aと同じならOK)
        print("LL^T = ")
        OutputMatrix(LL, n, n)
        print("")
    else:  # 修正コレスキー分解，不完全コレスキー分解のチェック
        D = zeros2(n, n)
        for i in range(n):
            D[i][i] = d[i]
        print("D = ")
        OutputMatrix(D, n, n)
        print("")

        # 分解結果のチェック用
        LD = zeros2(n, n)
        LDL = zeros2(n, n)
        MulMatrix(L, D, LD, n)
        MulMatrix(LD, transpose(L, n), LDL, n)

        # 分解した結果を掛け合わせた行列を画面表示(これが元の行列Aと同じならOK)
        print("LDL^T = ")
        OutputMatrix(LDL, n, n)
        print("")

    return 0


if __name__ == "__main__":
    main()
