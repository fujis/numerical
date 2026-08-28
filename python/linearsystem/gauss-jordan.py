# -*- coding: utf-8 -*-
"""!
@file gauss-jordan.py

@brief ガウス-ジョルダン(Gauss-Jordan)の消去法

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
# ガウス・ジョルダン法による逆行列計算
#  - ガウスの消去法を拡張して逆行列計算に用いる
# -----------------------------------------------------------------------------
def GaussJordan(A, n):
    """!
    ガウス・ジョルダン法(ピボット選択なし)
    @param[inout] A n×2nの拡張行列
    @param[in] n n元連立一次方程式
    """
    # 拡張行列の右半分を単位行列にする
    for i in range(n):
        for j in range(n):
            A[i][j+n] = 1.0 if i == j else 0.0

    # ガウス・ジョルダン法(Gauss-Jordan method)で逆行列計算
    for k in range(n):
        akk = A[k][k]
        # 対角要素を1にするために，k行目のすべての要素をa_kkで割る
        for j in range(2*n):
            A[k][j] /= akk

        # k列目の非対角要素を0にする
        for i in range(n):
            if i == k:
                continue
            aik = A[i][k]
            for j in range(2*n):
                A[i][j] -= A[k][j]*aik

        print("k = " + str(k))
        OutputMatrix(A, n, 2*n)
        print("")

    return 0


def Pivoting(A, n, k):
    """!
    ピボット選択(Pivoting)
     - 行入れ替えだけの部分的ピボッティング
    @param[inout] A n×2nの拡張行列
    @param[in] n n元連立一次方程式
    @param[in] k 対象行
    @return 交換した行
    """
    # k行目以降でk列目の絶対値が最も大きい要素を持つ行を検索
    p = k              # 絶対値が最大の行
    am = abs(A[k][k])  # 最大値
    for i in range(k+1, n):
        if abs(A[i][k]) > am:
            p = i
            am = abs(A[i][k])
    # k != pならば行を交換
    if k != p:
        A[k], A[p] = A[p], A[k]

    return p


def GaussJordanWithPivoting(A, n):
    """!
    ガウス・ジョルダン法(ピボット選択あり)
    @param[inout] A n×2nの拡張行列
    @param[in] n n元連立一次方程式
    """
    # 拡張行列の右半分を単位行列にする
    for i in range(n):
        for j in range(n):
            A[i][j+n] = 1.0 if i == j else 0.0

    # ガウス・ジョルダン法(Gauss-Jordan method)で逆行列計算
    ptable = [k for k in range(n)]
    for k in range(n):
        # ピボット選択
        p = Pivoting(A, n, k)
        if p != k:
            ptable[k], ptable[p] = ptable[p], ptable[k]

        akk = A[k][k]
        # 対角要素を1にするために，k行目のすべての要素をa_kkで割る
        for j in range(2*n):
            A[k][j] /= akk

        # k列目の非対角要素を0にする
        for i in range(n):
            if i == k:
                continue
            aik = A[i][k]
            for j in range(2*n):
                A[i][j] -= A[k][j]*aik

    return 0


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込む
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    A = ReadMatrix(fn, ",")
    A0 = [list(row) for row in A]  # チェック用に元の行列を確保(Pythonでは中身までコピーする必要あり)

    n = len(A)

    # 読み込んだ行列を確認用に画面表示
    print("A(" + str(n) + " x " + str(n) + ") = ")
    OutputMatrix(A, n, n)
    print("")

    # 拡張行列を作成
    for i in range(n):
        A[i] = A[i][0:n] + [0.0]*n  # 各行をn要素からなる左半分と2n要素にするための右半分に分ける
        for j in range(n):
            A[i][j+n] = 1.0 if i == j else 0.0

    # ガウスジョルダンで逆行列を求める
    # GaussJordan(A, n)              # ピボッティングなし
    GaussJordanWithPivoting(A, n)    # ピボッティングあり

    # 拡張行列A(nx2n)全体の画面表示
    OutputMatrix(A, n, 2*n)
    print("")

    # 拡張行列から逆行列部分だけを抽出
    invA = zeros2(n, n)
    for i in range(n):
        invA[i] = A[i][n:]

    # 逆行列部分のみの画面表示
    print("A^-1 = ")
    OutputMatrix(invA, n, n)
    print("")

    # 逆行列チェック
    C = zeros2(n, n)
    MulMatrix(A0, invA, C, n)

    # 元の行列と計算した逆行列を掛けた結果の行列を画面表示(これがほぼ単位行列ならOK)
    print("A A^-1 = ")
    OutputMatrix(C, n, n)
    print("")

    return 0


if __name__ == "__main__":
    main()
