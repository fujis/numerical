# -*- coding: utf-8 -*-
"""!
@file conjugate-gradient.py

@brief 共役勾配法(conjugate gradient method)

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
# (前処理付き)共役勾配法(CG法)による線形システムソルバ
# -----------------------------------------------------------------------------
def ModifiedCholeskyDecomp2(A, L, d, n):
    """!
    修正コレスキー分解(modified Cholesky decomposition)
     - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
     - l_ii * d_i = 1とした場合 <- この部分がcholesky.pyと異なるので注意
     - L: i > jの要素が非ゼロで対角成分は1
    @param[in] A n×nの対称行列
    @param[out] L 対角成分が1の下三角行列
    @param[out] d 対角行列(対角成分のみ)
    @param[in] n 行列の大きさ
    @return 0:成功,1:失敗
    """
    if n <= 0:
        return 1

    L[0][0] = A[0][0]
    d[0] = 1.0/L[0][0]

    for i in range(1, n):
        for j in range(i+1):
            lld = A[i][j]
            for k in range(j):
                lld -= L[i][k]*L[j][k]*d[k]
            L[i][j] = lld
        d[i] = 1.0/L[i][i]

    return 0


def IncompleteCholeskyDecomp2(A, L, d, n):
    """!
    不完全コレスキー分解(incomplete Cholesky decomposition)
     - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
     - l_ii * d_i = 1とした場合 <- この部分がcholesky.pyと異なるので注意
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

    L[0][0] = A[0][0]
    d[0] = 1.0/L[0][0]

    for i in range(1, n):
        for j in range(i+1):
            if abs(A[i][j]) < 1.0e-10:
                continue

            lld = A[i][j]
            for k in range(j):
                lld -= L[i][k]*L[j][k]*d[k]
            L[i][j] = lld

        d[i] = 1.0/L[i][i]

    return 0


def CGSolver(A, b, x, n, max_iter, eps):
    """!
    共役勾配法によりA・x=bを解く
    @param[in] A n×n正値対称行列
    @param[in] b 右辺ベクトル
    @param[out] x 結果ベクトル(リストなので破壊的に更新される)
    @param[in] n 行列の大きさ
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    if n <= 0:
        return 1, max_iter, eps
    set_precision(8)

    r = zeros(n)
    p = zeros(n)
    y = zeros(n)
    for i in range(n):
        x[i] = 0.0

    # 第0近似解に対する残差の計算
    for i in range(n):
        ax = 0.0
        for j in range(n):
            ax += A[i][j]*x[j]
        r[i] = b[i]-ax
        p[i] = r[i]

    rr0 = dot(r, r, n)

    e = 0.0
    k = 0
    for k in range(max_iter):
        # y = AP の計算
        for i in range(n):
            y[i] = dot(A[i], p, n)

        # alpha = r*r/(P*AP)の計算
        pap = dot(p, y, n)
        if abs(pap) < 1e-6:
            break  # 右辺項bが全て0かつxも全て0だと0割になるのでその対策(b=0のときのx=0は自明な解)
        alpha = rr0/pap

        # 解x、残差rの更新
        for i in range(n):
            x[i] += alpha*p[i]
            r[i] -= alpha*y[i]

        # 確認のため現在の解を画面出力
        s = str(k) + " : "
        for i in range(n):
            s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
        print(s)

        # (r*r)_(k+1)の計算
        rr1 = dot(r, r, n)

        # 収束判定 (||r||<=eps)
        e = math.sqrt(rr1)
        if e < eps:
            k += 1
            break

        # βの計算とPの更新
        beta = rr1/rr0
        for i in range(n):
            p[i] = r[i]+beta*p[i]

        # (r*r)_(k+1)を次のステップのために確保しておく
        rr0 = rr1

    return 0, k+1, e


def ICRes(L, d, r, u, n):
    """!
    (LDL^T)^-1 r の計算
    @param[in] L,d IC分解で得られた下三角行列と対角行列(対角成分のみのベクトル)
    @param[in] r 残差ベクトル
    @param[out] u (LDL^T)^-1 rを計算した結果
    @param[in] n 行列の大きさ
    """
    y = zeros(n)
    for i in range(n):
        rly = r[i]
        for j in range(i):
            rly -= L[i][j]*y[j]
        y[i] = rly/L[i][i]

    for i in range(n-1, -1, -1):
        lu = 0.0
        for j in range(i+1, n):
            lu += L[j][i]*u[j]
        u[i] = y[i]-d[i]*lu


def ICCGSolver(A, b, x, n, max_iter, eps):
    """!
    不完全コレスキー分解による前処理付共役勾配法によりA・x=bを解く
    @param[in] A n×n正値対称行列
    @param[in] b 右辺ベクトル
    @param[out] x 結果ベクトル(リストなので破壊的に更新される)
    @param[in] n 行列の大きさ
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    if n <= 0:
        return 1, max_iter, eps
    set_precision(8)

    r = zeros(n)
    p = zeros(n)
    y = zeros(n)
    r2 = zeros(n)

    # 初期値を設定
    for i in range(n):
        x[i] = 0.0

    # コレスキー分解
    d = zeros(n)
    L = zeros2(n, n)
    # ModifiedCholeskyDecomp2(A, L, d, n)
    IncompleteCholeskyDecomp2(A, L, d, n)

    # 第0近似解に対する残差の計算
    for i in range(n):
        ax = 0.0
        for j in range(n):
            ax += A[i][j]*x[j]
        r[i] = b[i]-ax

    # p_0 = (LDL^T)^-1 r_0 の計算
    ICRes(L, d, r, p, n)

    # C++版ではここだけfloat(単精度)を使っているが，Pythonのfloatは倍精度になる
    rr0 = dot(r, p, n)

    e = 0.0
    k = 0
    for k in range(max_iter):
        s = str(k)

        # y = AP の計算
        for i in range(n):
            y[i] = dot(A[i], p, n)

        # alpha = r*r/(P*AP)の計算
        pap = dot(p, y, n)
        if abs(pap) < 1e-6:
            break  # 右辺項bが全て0かつxも全て0だと0割になるのでその対策(b=0のときのx=0は自明な解)
        alpha = rr0/pap

        # 解x、残差rの更新
        for i in range(n):
            x[i] += alpha*p[i]
            r[i] -= alpha*y[i]

        # 確認のため現在の解を画面出力
        # s = str(k) + " : "
        # for i in range(n):
        #     s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
        # print(s)

        # (r*r)_(k+1)の計算
        ICRes(L, d, r, r2, n)
        rr1 = dot(r, r2, n)

        # 収束判定 (||r||<=eps)
        e = math.sqrt(rr1)
        s += ", " + fmt(e/n)
        print(s)
        if e < eps:
            k += 1
            break

        # βの計算とPの更新
        beta = rr1/rr0
        for i in range(n):
            p[i] = r2[i]+beta*p[i]

        # (r*r)_(k+1)を次のステップのために確保しておく
        rr0 = rr1

    return 0, k, e


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

    # 共役勾配法を用いて連立1次方程式を解く
    max_iter = 100
    eps = 1e-6
    x = zeros(n)
    # ret, max_iter, eps = CGSolver(A, b, x, n, max_iter, eps)
    # max_iter = 100; eps = 1e-6
    ret, max_iter, eps = ICCGSolver(A, b, x, n, max_iter, eps)

    # 結果の画面表示
    s = ""
    for i in range(n):
        s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
    print(s)

    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
