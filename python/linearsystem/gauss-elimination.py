# -*- coding: utf-8 -*-
"""!
@file gauss-elimination.py

@brief ガウスの消去法(Gaussian elimination)

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
# 関数
# -----------------------------------------------------------------------------
def GaussElimination(A, n):
    """!
    ガウスの消去法(ピボット交換なし)
    @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
    @param[in] n n元連立一次方程式
    """
    # 前進消去(forward elimination)
    #  - 対角要素をのぞいた左下要素をすべて0にする(上三角行列にする)
    for k in range(n-1):
        akk = A[k][k]
        for i in range(k+1, n):
            aik = A[i][k]
            for j in range(k, n+1):  # 確認のため左下要素が0になるようにj=kとしたが，実際にはj=k+1でよい
                A[i][j] = A[i][j]-aik*(A[k][j]/akk)
        print("k = " + str(k))
        OutputMatrix(A, n, n+1)
        print("")

    # 後退代入(back substitution)
    #  - x_nの解はb_n/a_nn，x_nをさらにn-1行の式に代入することでx_(n-1)を求める．
    #  - この作業を1行目まで続けることですべての解を得る．
    A[n-1][n] = A[n-1][n]/A[n-1][n-1]
    for i in range(n-2, -1, -1):
        ax = 0.0
        for j in range(i+1, n):
            ax += A[i][j]*A[j][n]
        A[i][n] = (A[i][n]-ax)/A[i][i]

        print("i = " + str(i))
        OutputMatrix(A, n, n+1)
        print("")

    return 0


def Pivoting(A, n, k):
    """!
    ピボット選択(Pivoting)
     - 行入れ替えだけの部分的ピボッティング
    @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列
    @param[in] n n元連立一次方程式
    @param[in] k 対象行
    """
    # k行目以降でk列目の絶対値が最も大きい要素を持つ行を検索
    p = k              # 絶対値が最大の行
    am = abs(A[k][k])  # 最大値
    for i in range(k+1, n):
        if abs(A[i][k]) > am:
            p = i
            am = abs(A[i][k])
    # k != pならば行を交換(ピボット選択)
    if k != p:
        A[k], A[p] = A[p], A[k]


def GaussEliminationWithPivoting(A, n):
    """!
    ガウスの消去法(行に関するピボット選択(部分ピボッティング)あり)
    @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
    @param[in] n n元連立一次方程式
    """
    # 前進消去(forward elimination)
    #  - 対角要素をのぞいた左下要素をすべて0にする(上三角行列にする)
    for k in range(n-1):
        # ピボット選択
        Pivoting(A, n, k)

        akk = A[k][k]
        for i in range(k+1, n):
            aik = A[i][k]
            for j in range(k, n+1):  # 確認のため左下要素が0になるようにj=kとしたが，実際にはj=k+1でよい
                A[i][j] = A[i][j]-aik*(A[k][j]/akk)
        print("k = " + str(k))
        OutputMatrix(A, n, n+1)
        print("")

    # 後退代入(back substitution)
    #  - x_nの解はb_n/a_nn，x_nをさらにn-1行の式に代入することでx_(n-1)を求める．
    #  - この作業を1行目まで続けることですべての解を得る．
    A[n-1][n] = A[n-1][n]/A[n-1][n-1]
    for i in range(n-2, -1, -1):
        ax = 0.0
        for j in range(i+1, n):
            ax += A[i][j]*A[j][n]
        A[i][n] = (A[i][n]-ax)/A[i][i]

    return 0


def ScalingForGauss(A, n):
    """!
    スケーリング処理
     - ガウスの消去法の前処理として使用
    @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列
    @param[in] n n元連立一次方程式
    @return 1:成功
    """
    # ガウスの消去法において桁落ちを防ぐために，
    # 各係数を各行の最大値で割ることで正規化する．
    for i in range(n):
        # 行要素で絶対値最大のものを検索
        am = abs(A[i][0])
        for j in range(1, n):  # 右辺項b_i=a_i,nは含まない
            if abs(A[i][j]) > am:
                am = abs(A[i][j])

        # すべての行要素を最大値で割る
        for j in range(n+1):
            A[i][j] /= am

    return 1


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込む場合
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    A = ReadMatrix(fn, ",")

    # 2次元リストに初期値として直接値を設定する場合
    # A = [[2, 1, 3, 9],
    #      [1, 3, 2, 1],
    #      [3, 4, 3, 4]]

    # 行列のサイズの取得と確認のための画面表示
    n = len(A)     # n元連立一次方程式
    m = len(A[0])
    print("A(" + str(n) + " x " + str(m) + ") = ")
    OutputMatrix(A, n, n+1)
    print("")

    # 最後に結果をチェックするために行列を退避しておく
    #  - Pythonでは A0 = A としても中身がコピーされない(同じリストを指すだけ)ので，
    #    各行をコピーして新しいリストを作る必要がある
    A0 = [list(row) for row in A]

    # ガウスの消去法で線形システムを解く
    # ScalingForGauss(A, n)              # スケーリング処理
    # GaussElimination(A, n)             # ピボッティングなし
    GaussEliminationWithPivoting(A, n)   # ピボッティングあり

    # 結果の表示
    s = ""
    for i in range(n):
        s += "x" + str(i) + " = " + fmt(A[i][n]) + ("" if i == n-1 else ", ")
    print(s)

    # 結果のチェック
    x = zeros(n)
    y = zeros(n)
    for i in range(n):
        x[i] = A[i][n]
    MulMatrixVector(A0, x, y, n)
    ok = True
    for i in range(n):
        # 桁落ち/丸め誤差があるので実数型で
        # if A0[i][n] != y[i] のような比較をしないように！
        if abs(A0[i][n]-y[i]) > 1e-6:
            ok = False
            print("x_" + str(i) + " is wrong result. ")
    if ok:
        print("solution check : all clear")

    return 0


if __name__ == "__main__":
    main()
