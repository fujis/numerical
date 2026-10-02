# -*- coding: utf-8 -*-
"""!
@file jacobi.py

@brief 固有値・固有ベクトル
       ヤコビ(Jacobi)法

@author Makoto Fujisawa
@date 2012-06 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *


def eigen_jacobi(a, v, n, max_iter, eps):
    """!
    Jacobi法による固有値の算出
    @param[inout] a 実対称行列．計算後，対角要素に固有値が入る
    @param[out] v 固有ベクトル(aと同じサイズ)
    @param[in] n 行列のサイズ(n×n)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    bim = zeros(n)
    bjm = zeros(n)

    for i in range(n):
        for j in range(n):
            v[i][j] = 1.0 if i == j else 0.0

    e = 0.0
    k = 0
    for k in range(max_iter):
        # 非対角要素で絶対値が最大の要素をa_ijとして0にする
        i = 0
        j = 0
        x = 0.0
        for ia in range(n):
            for ja in range(n):
                if ia != ja and abs(a[ia][ja]) > x:
                    i = ia
                    j = ja
                    x = abs(a[ia][ja])

        # sinθ,cosθを算出
        aii = a[i][i]
        ajj = a[j][j]
        aij = a[i][j]

        alpha = (aii-ajj)/2.0
        beta = math.sqrt(alpha*alpha+aij*aij)

        ct = math.sqrt((1.0+abs(alpha)/beta)/2.0)                  # cosθ
        st = (1.0 if (aii-ajj) >= 0.0 else -1.0)*aij/(2.0*beta*ct)  # sinθ

        # A = PAPの計算
        for m in range(n):
            if m == i or m == j:
                continue

            aim = a[i][m]
            ajm = a[j][m]

            bim[m] = aim*ct+ajm*st
            bjm[m] = -aim*st+ajm*ct

        bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st
        bij = 0.0

        bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct
        bji = 0.0

        for m in range(n):
            a[i][m] = a[m][i] = bim[m]
            a[j][m] = a[m][j] = bjm[m]
        a[i][i] = bii
        a[i][j] = bij
        a[j][j] = bjj
        a[j][i] = bji

        # V = PVの計算
        for m in range(n):
            vmi = v[m][i]
            vmj = v[m][j]

            bim[m] = vmi*ct+vmj*st
            bjm[m] = -vmi*st+vmj*ct
        for m in range(n):
            v[m][i] = bim[m]
            v[m][j] = bjm[m]

        # 非対角要素の絶対値の和で収束を判定
        e = 0.0
        for ja in range(n):
            for ia in range(n):
                if ia != ja:
                    e += abs(a[ja][ia])
        if e < eps:
            k += 1
            break

    return 1, k, e


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

    V = zeros2(n, n)

    max_iter = 100
    eps = 1e-6
    ret, max_iter, eps = eigen_jacobi(A, V, n, max_iter, eps)

    # 固有値の表示
    s = "e = ("
    for i in range(n):
        s += fmt(A[i][i]) + (")" if i == n-1 else ", ")
    print(s)

    # 固有ベクトルの表示
    for j in range(n):
        s = "v" + str(j) + " = ("
        for i in range(n):
            s += fmt(V[i][j]) + (")" if i == n-1 else ", ")
        print(s)
    print("")

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
