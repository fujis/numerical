# -*- coding: utf-8 -*-
"""!
@file pca.py

@brief 主成分分析(PCA:Principal Component Analysis)
       共分散行列の固有値・固有ベクトルをQR法+逆反復法で求める

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

# QR分解・QR法・逆反復法の関数は qr.py にあるものと全く同じなのでそちらから読み込む
#  - C++版では同じコードが qr.cpp と pca.cpp の両方に書かれている
from qr import lu_decomp, lu_solver, qr_decomposition, inverse_iteration, qr


def covariance_matrix2(data):
    """!
    2次元座標データから共分散行列を作成
    @param[in] data 2次元座標データ(rxPoint2のリスト)
    @return (A, c) 2x2の共分散行列とデータの平均値(重心)
    """
    m = len(data)
    n = 2

    # データの平均値(重心)の計算
    c = rxPoint2(0.0, 0.0)
    for p in data:
        c.x += p.x
        c.y += p.y
    c.x /= m
    c.y /= m

    # 共分散行列の計算
    A = zeros2(n, n)
    for p in data:
        A[0][0] += (p.x-c.x)*(p.x-c.x)
        A[1][1] += (p.y-c.y)*(p.y-c.y)
        A[0][1] += (p.x-c.x)*(p.y-c.y)
        A[1][0] += (p.y-c.y)*(p.x-c.x)
    A[0][0] /= m
    A[1][1] /= m
    A[0][1] /= m
    A[1][0] /= m

    return A, c


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから点データを読み込む
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "points.txt")
    data = Read2d(fn)
    if data is None:
        print("file read error!")
        return 1
    m = len(data)
    n = 2

    # 共分散行列の計算
    A, c = covariance_matrix2(data)
    A0 = [list(row) for row in A]

    print("number of points = " + str(m))
    print("center = (" + fmt(c.x) + ", " + fmt(c.y) + ")")
    print("covariance matrix = ")
    OutputMatrix(A, n, n)
    print("")

    # QR法で固有値を求める
    max_iter = 100
    eps = 1e-6
    A, max_iter, eps = qr(A, n, max_iter, eps, verbose=False)

    # 固有値の表示
    s = "e = ("
    lambda_ = zeros(n)
    for i in range(n):
        s += fmt(A[i][i]) + (")" if i == n-1 else ", ")
        lambda_[i] = A[i][i]
    print(s)
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    # 逆反復法で固有ベクトルを計算
    max_iter = 100
    eps = 1e-6
    v, max_iter, eps = inverse_iteration(A0, lambda_, n, max_iter, eps)

    # 固有ベクトルの表示
    #  - 第1主成分(最大固有値に対応する固有ベクトル)がデータの分散が最大となる方向
    for i in range(n):
        s = "v" + str(i+1) + " = ("
        for j in range(n):
            s += fmt(v[i][j]) + (")" if j == n-1 else ", ")
        print(s)
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))

    return 0


if __name__ == "__main__":
    main()
