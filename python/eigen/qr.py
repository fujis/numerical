# -*- coding: utf-8 -*-
"""!
@file qr.py

@brief 固有値・固有ベクトル
       QR法

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


def HouseholderTransformationForQR(a, n):
    """!
    ハウスホルダー変換で実対称行列をヘッセンベルグ行列に変換
    @param[inout] a 元の行列(n×n)．変換された三重対角行列を格納
    @param[in] n 行列のサイズ
    """
    b = zeros2(n, n)
    p = zeros2(n, n)
    q = zeros2(n, n)
    u = zeros(n)

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

        # Pの計算
        for i in range(k+1, n):
            for j in range(i, n):
                if j == i:
                    p[i][j] = 1.0-2.0*u[i]*u[i]
                else:
                    p[i][j] = -2.0*u[i]*u[j]
                    p[j][i] = p[i][j]

        # PAの計算
        for i in range(k+1, n):
            for j in range(k+1, n):
                q[i][j] = 0.0
                for m in range(k+1, n):
                    q[i][j] += p[i][m]*a[m][j]

        # A = PAP^Tの計算
        for i in range(k+1):
            b[i][k] = a[i][k]
        b[k+1][k] = s
        for i in range(k+2, n):
            b[i][k] = 0.0
        for j in range(k+1, n):
            for i in range(k+1):
                b[i][j] = 0.0
                for m in range(k+1, n):
                    b[i][j] += a[i][m]*p[j][m]
            for i in range(k+1, n):
                b[i][j] = 0.0
                for m in range(k+1, n):
                    b[i][j] += q[i][m]*p[j][m]

        for i in range(n):
            for j in range(n):
                a[i][j] = b[i][j]


def EigenQR(h, n, max_iter, eps):
    """!
    QR法による固有値の算出
    @param[inout] h 元の行列(ヘッセンベルグ行列)．計算後，対角要素に固有値が入る
    @param[in] n 行列のサイズ(n×n)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    r = zeros2(n, n)
    q = zeros2(n, n)
    t = zeros2(n, n)
    u = zeros(n)
    v = zeros(n)

    # R = H : ヘッセンベルグ行列
    for i in range(n):
        for j in range(n):
            r[i][j] = h[i][j]

    e = 0.0
    l = 0
    for l in range(max_iter):
        # Q=I (単位行列)
        for i in range(n):
            for j in range(n):
                q[i][j] = 1.0 if i == j else 0.0

        for k in range(n-1):
            # sinθ,cosθの計算
            alpha = math.sqrt(r[k][k]*r[k][k]+r[k+1][k]*r[k+1][k])
            if abs(alpha) < 1e-8:
                continue

            c = r[k][k]/alpha
            s = -r[k+1][k]/alpha

            # Rの計算
            for j in range(k+1, n):
                u[j] = c*r[k][j]-s*r[k+1][j]
                v[j] = s*r[k][j]+c*r[k+1][j]
            r[k][k] = alpha
            r[k+1][k] = 0.0
            for j in range(k+1, n):
                r[k][j] = u[j]
                r[k+1][j] = v[j]

            # Qの計算
            for j in range(k+1):
                u[j] = c*q[k][j]
                v[j] = s*q[k][j]
            q[k][k+1] = -s
            q[k+1][k+1] = c
            for j in range(k+1):
                q[k][j] = u[j]
                q[k+1][j] = v[j]

        # RQの計算
        for i in range(n):
            for j in range(n):
                rq = 0.0
                for m in range(n):
                    rq += r[i][m]*q[j][m]

                t[i][j] = rq

        # 収束判定
        e = 0.0
        for i in range(n):
            e += abs(t[i][i]-h[i][i])
        if e < eps:
            l += 1
            break

        for i in range(n):
            for j in range(n):
                r[i][j] = t[i][j]
                h[i][j] = t[i][j]

    return 1, l, e


def lu_decomp(A, n):
    """!
    LU分解(ピボット交換なし)
     - 行列A(n×n)を下三角行列(L)と上三角行列(U)に分解する
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
            A[i][j] = div(lu, A[i][i])  # 特異行列だと0割になるのでdiv()を使う(C++と同じくinf/nanになる)

    return 0


def lu_solver(A, b, x, n):
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
        x[i] = div(bly, A[i][i])  # 特異行列だと0割になるのでdiv()を使う(C++と同じくinf/nanになる)

    # 後退代入(back substitution)
    #  UX=YからXを計算
    for i in range(n-1, -1, -1):
        yux = x[i]
        for j in range(i+1, n):
            yux -= A[i][j]*x[j]
        x[i] = yux

    return 0


def qr_decomposition(A, n):
    """!
    ハウスホルダー変換でQR分解 (A=QR)
    @param[in] A 元の行列(n×n)
    @param[in] n 行列のサイズ
    @return (Q, R) 分解された直行行列Qと上三角行列R
    """
    u = zeros(n)
    H = zeros2(n, n)             # 行列H
    R = [list(row) for row in A]  # A^(0)
    Q = zeros2(n, n)
    unit(Q, n)                   # Qを単位行列で初期化

    for k in range(n-1):
        # a_kk^(k+1)の計算 (ベクトルxの大きさ)
        akk = 0.0
        for i in range(k, n):
            akk += R[i][k]*R[i][k]
        akk = math.sqrt(akk)

        # u=(x1-y1)/|x1-y1|の計算
        l = 0.0
        for i in range(k, n):
            u[i] = R[i][k]-(akk if i == k else 0.0)
            l += u[i]*u[i]
        if abs(l) < 1e-10:
            break
        l = math.sqrt(l)
        for i in range(k, n):
            u[i] /= l

        # ハウスホルダー行列H^(k)の計算(H=I-2uu^T)
        unit(H, n)  # Hを単位行列で初期化
        for i in range(k, n):
            for j in range(k, n):
                H[i][j] -= 2*u[i]*u[j]

        # A^(k+1) = H^(k) A^(k)
        R = mul_mm(H, R, n)

        # Q^(k+1) = Q^(k) (H^(k))^T
        Q = mul_mm(Q, transpose(H, n), n)

    return Q, R


def inverse_iteration(A, lambda_, n, max_iter, eps):
    """!
    逆反復法による固有ベクトルの算出
    @param[in] A 元の行列
    @param[in] lambda_ 固有値のリスト
    @param[in] n 行列のサイズ(n×n)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (固有ベクトルのリストv, 平均反復数, 平均誤差)
    """
    if n <= 0:
        return None, max_iter, eps

    LU = [list(row) for row in A]  # LU分解した行列を格納するための配列
    y = zeros(n, 1.0)
    y1 = zeros(n)  # y^(k), y^(k+1)格納用
    sum_e = 0.0
    sum_k = 0
    v = [None]*n

    for l in range(n):
        # (A-λI)を計算
        for i in range(n):
            for j in range(n):
                LU[i][j] = A[i][j]-(lambda_[l] if i == j else 0.0)

        # (A-λI)をLU分解
        lu_decomp(LU, n)

        for i in range(n):
            y[i] = 1.0
        l2 = dot(y, y, n)
        mu0 = 0.0
        e = 0.0
        k = 0
        for k in range(max_iter):
            # |y^(k)|=1となるように正規化
            ly = math.sqrt(l2)        # |y^(k)|の計算
            y = mul_sv(1.0/ly, y, n)  # |y^(k)|で割る(1/|y|を掛ける)

            # y^(k+1) = B y^(k)の計算(By^(k+1)=y^(k)をLU分解で解く)
            lu_solver(LU, y, y1, n)

            # 固有値1/|λ-λi|の計算(|y^(k)|=1で正規化されていることが前提)
            mu = dot(y, y1, n)

            # 収束判定
            e = abs(div(mu-mu0, mu))
            if e < eps:
                break

            y = list(y1)
            l2 = dot(y, y, n)
            mu0 = mu

        v[l] = mul_sv(div(1.0, math.sqrt(dot(y, y, n))), y, n)  # 最後に正規化しておく
        sum_e += e
        sum_k += k

    return v, sum_k//n, sum_e/n


def qr(A, n, max_iter, eps, verbose=True):
    """!
    QR法による固有値の算出
    @param[in] A 元の行列
    @param[in] n 行列のサイズ(n×n)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @param[in] verbose 反復ごとの確認用表示のON/OFF
    @return (計算後の行列A(対角要素に固有値が入る), 実際の反復数, 実際の誤差)
    """
    if n <= 0:
        return A, max_iter, eps
    lambda_ = zeros(n)  # 収束判定用
    e = eps*n
    k = 0
    for k in range(max_iter):
        # A^(k) ⇐ Q^(k) R^(k)
        Q, R = qr_decomposition(A, n)

        # A^(k+1) ⇐ R^(k) Q^(k)
        A = mul_mm(R, Q, n)

        # 確認用表示
        if verbose:
            s = str(k) + " : lambda = "
            for i in range(n):
                s += fmt(A[i][i]) + ("" if i == n-1 else ", ")
            print(s)

        # 収束判定
        e = 0.0
        for i in range(n):
            e += abs(lambda_[i]-A[i][i])
        if e/n < eps:
            break

        for i in range(n):
            lambda_[i] = A[i][i]

    return A, k, e/n


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
    A0 = [list(row) for row in A]
    print("")

    # QR分解
    Q, R = qr_decomposition(A, n)
    print("Q = ")
    OutputMatrix(Q, n, n)
    print("R = ")
    OutputMatrix(R, n, n)

    QR = mul_mm(Q, R, n)
    print("QR = ")
    OutputMatrix(QR, n, n)

    QQ = mul_mm(Q, transpose(Q, n), n)
    print("QQ = ")
    OutputMatrix(QQ, n, n)
    print("")

    # QR法で固有値を求める
    max_iter = 100
    eps = 1e-6
    A, max_iter, eps = qr(A, n, max_iter, eps)

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
    for i in range(n):
        s = "v" + str(i+1) + " = ("
        for j in range(n):
            s += fmt(v[i][j]) + (")" if j == n-1 else ", ")
        print(s)
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))

    print("")

    # HouseholderTransformationForQR(A, n)
    #
    # max_iter = 100
    # eps = 1e-6
    # ret, max_iter, eps = EigenQR(A, n, max_iter, eps)
    #
    # # 固有値の表示
    # s = "e = ("
    # for i in range(n):
    #     s += fmt(A[i][i]) + (")" if i == n-1 else ", ")
    # print(s)
    #
    # print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    # print("")

    return 0


if __name__ == "__main__":
    main()
