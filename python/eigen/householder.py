# -*- coding: utf-8 -*-
"""!
@file householder.py

@brief 固有値・固有ベクトル
       ハウスホルダー(Householder)法

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


def HouseholderTransformation(a, n):
    """!
    ハウスホルダー変換で実対称行列を三重対角行列に変換
    @param[inout] a 元の行列(n×n)．変換された三重対角行列を格納
    @param[in] n 行列のサイズ
    """
    u = zeros(n)
    v = zeros(n)
    q = zeros(n)

    for k in range(n-2):
        # sの計算
        s = 0.0
        for i in range(k+1, n):
            s += a[i][k]*a[i][k]
        s = -(1 if a[k+1][k] >= 0 else -1)*math.sqrt(s)

        # |x-y|の計算
        alpha = math.sqrt(2.0*s*(s-a[k+1][k]))
        if abs(alpha) < 1e-8:
            continue

        # uの計算
        u[k+1] = (a[k+1][k]-s)/alpha
        for i in range(k+2, n):
            u[i] = a[i][k]/alpha

        # Auの計算
        q[k] = alpha/2.0
        for i in range(k+1, n):
            q[i] = 0.0
            for j in range(k+1, n):
                q[i] += a[i][j]*u[j]

        # v=2(Au-uu^T(Au))の計算
        alpha = 0.0
        for i in range(k+1, n):
            alpha += u[i]*q[i]
        v[k] = 2.0*q[k]
        for i in range(k+1, n):
            v[i] = 2.0*(q[i]-alpha*u[i])

        # A = PAP = A-uv^T-vu^Tの計算
        a[k][k+1] = a[k+1][k] = s
        for i in range(k+2, n):
            a[k][i] = 0.0
            a[i][k] = 0.0
        for i in range(k+1, n):
            a[i][i] = a[i][i]-2*u[i]*v[i]
            for j in range(i+1, n):
                a[i][j] = a[i][j]-u[i]*v[j]-v[i]*u[j]
                a[j][i] = a[i][j]


def NCS(lambda_, A, n):
    """!
    スツルム関数列の符号変化回数をカウント
     - λが大きいときのf_k(λ)の計算でオーバーフローが発生するのを防ぐために，
       f_k でなく g_k = f_k/f_(k-1) を使う．
     - g_k < 0 となる回数を数える
    @param[in] lambda_ 固有値候補
    @param[in] A 三重対角行列(n×n)
    @param[in] n 行列のサイズ
    @return 符号変化回数
    """
    g0 = lambda_-A[0][0]
    ncs = 1 if g0 < 0 else 0
    for k in range(1, n):
        if g0 == 0:
            g0 = RX_FEQ_EPS

        a = A[k][k]
        b = A[k-1][k]
        g = lambda_-a-(b*b)/g0
        if g < 0:
            ncs += 1

        g0 = g

    return ncs


def CFunc(lambda_, A, n):
    """!
    三重対角行列の特性方程式の漸化式(スツルム関数列)を計算
     - ニュートン法で使うため，勾配値も同時に計算
    @param[in] lambda_ 固有値候補(変数)
    @param[in] A 三重対角行列(n×n)
    @param[in] n 行列のサイズ
    @return (f, df) スツルム関数列f_nの値とその勾配値
    """
    if n == 0:
        return 1.0, 0.0
    if n == 1:
        return lambda_-A[0][0], 1.0

    f0 = 1.0                # f_0
    f1 = lambda_-A[0][0]    # f_1
    df0 = 0.0
    df1 = 1.0

    f = f1
    df = df1
    for k in range(1, n):
        a = A[k][k]
        b = A[k-1][k]
        f = (lambda_-a)*f1-b*b*f0
        df = (lambda_-a)*df1-b*b*df0+f1

        f0 = f1
        f1 = f
        df0 = df1
        df1 = df

    return f, df


def Newton(x, eps, A, n):
    """!
    ニュートン法でf_n = 0となるλを探索
    @param[in] x 探索開始位置
    @param[in] eps 収束判定用許容誤差
    @param[in] A 三重対角行列(n×n)
    @param[in] n 行列のサイズ
    @return 根
    """
    iter = 0
    while True:
        x0 = x

        f, df = CFunc(x0, A, n)

        x = x0-f/df
        iter += 1
        if not (abs(x-x0) > eps and iter < 100):
            break

    return x


def EigenHouseholder(H, n, eps):
    """!
    ハウスホルダー法による固有値の算出
     - 実対称行列を三重対角行列に変換し，スツルムの定理と2分法,ニュートン法で固有値を求める
    @param[inout] H 実対称行列．計算後，対角要素に固有値が入る
    @param[in] n 行列のサイズ(n×n)
    @param[in] eps 許容誤差
    @return 1
    """
    #
    # ハウスホルダー変換で三重対角行列に変換
    #
    HouseholderTransformation(H, n)

    print("H(" + str(n) + " x " + str(n) + ") = ")
    OutputMatrix(H, n, n)
    print("")

    #
    # 三重対角行列の固有値をスツルム(Sturm)の定理と2分法,ニュートン法で算出
    #

    # 三重対角行列の対角成分と非対角成分
    a = zeros(n)
    b = zeros(n-1)
    for i in range(n-1):
        a[i] = H[i][i]
        b[i] = H[i+1][i]
    a[n-1] = H[n-1][n-1]

    # 探索範囲[alpha, beta]をゲルシュゴリン(Gerschgorin)の定理より算出
    alpha = a[0]-abs(b[0])
    beta = a[0]+abs(b[0])
    for i in range(1, n-1):
        alpha = RX_MIN(alpha, a[i]-abs(b[i])-abs(b[i-1]))
        beta = RX_MAX(beta, a[i]+abs(b[i])+abs(b[i-1]))
    alpha = RX_MIN(alpha, a[n-1]-abs(b[n-2]))
    beta = RX_MAX(beta, a[n-1]+abs(b[n-2]))

    # # ファイル出力
    # fo = open("_debug.txt", "w")
    # x1 = alpha
    # while x1 < beta:
    #     f, df = CFunc(x1, H, n)
    #     fo.write(fmt(x1) + ", " + fmt(f) + ", " + str(NCS(x1, H, n)) + "\n")
    #     x1 += 0.01
    # fo.close()

    ermax = abs(beta-alpha)/n

    m = 0
    p = 0
    xsl = zeros(n)
    xsr = zeros(n)
    nsl = [0]*n
    nsr = [0]*n

    # 探索範囲境界でのスツルム関数列の符号変化回数を計算
    xsl[p] = alpha
    nsl[p] = NCS(alpha, H, n)
    xsr[p] = beta
    nsr[p] = NCS(beta, H, n)

    # 二分法で符号変化回数の差が1となる領域を探索
    eigen = zeros(n)
    while p >= 0:
        if (nsl[p]-nsr[p] == 1) and (xsr[p]-xsl[p] < ermax):
            # 符号変化回数の差が1となったらニュートン法で固有値を探索
            x = 0.5*(xsl[p]+xsr[p])
            eigen[m] = Newton(x, eps, H, n)
            p -= 1
            m += 1

            if p < 0:
                break

        xl = xsl[p]
        nl = nsl[p]
        xr = xsr[p]
        nr = nsr[p]
        xm = 0.5*(xl+xr)
        nm = NCS(xm, H, n)
        p -= 1

        # 符号変化回数に変化があれば，中点を新たな探索境界にする
        if nl > nm:
            p += 1
            xsl[p] = xl
            nsl[p] = nl
            xsr[p] = xm
            nsr[p] = nm
        if nm > nr:
            p += 1
            xsl[p] = xm
            nsl[p] = nm
            xsr[p] = xr
            nsr[p] = nr

    for i in range(n):
        for j in range(n):
            H[i][j] = eigen[i] if i == j else 0.0

    return 1


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

    # ハウスホルダー法で固有値を算出
    eps = 1e-6
    EigenHouseholder(A, n, eps)

    # 固有値の表示
    s = "e = ("
    for i in range(n):
        s += fmt(A[i][i]) + (")" if i == n-1 else ", ")
    print(s)

    return 0


if __name__ == "__main__":
    main()
