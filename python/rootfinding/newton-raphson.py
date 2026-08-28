# -*- coding: utf-8 -*-
"""!
@file newton-raphson.py

@brief ニュートン・ラフソン法(Newton-Raphson method)
       もしくは単にニュートン法

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
from rx_funcs import *


# -----------------------------------------------------------------------------
# ニュートン法による求根問題の解法
# -----------------------------------------------------------------------------
def newton(func, dfunc, x, max_iter, eps):
    """!
    ニュートン・ラフソン法(Newton-Raphson method)
    @param[in] func 関数値を与える関数
    @param[in] dfunc 導関数値を与える関数
    @param[in] x 探索開始位置
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解x, 実際の反復数, 実際の誤差)
    """
    dx = 0.0
    k = 0
    for k in range(max_iter):
        # 現在の位置xにおける関数値と導関数の計算
        f = func(x)
        df = dfunc(x)

        # 確認用の画面出力
        print(str(k) + " : f(" + fmt(x) + ") = " + fmt(f))

        # 導関数の結果から次の位置を計算
        x = x-f/df

        # 収束判定
        dx = abs(f/df)
        if dx < eps or abs(f) < eps:
            break

    return 0, x, k, dx


def secant(func, x0, x1, max_iter, eps):
    """!
    セカント法(secant method, 割線法)
    @param[in] func 関数値を与える関数
    @param[in] x0,x1 探索開始位置
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解x0, 解x1, 実際の反復数, 実際の誤差)
    """
    f0 = func(x0)
    f1 = func(x1)

    dx = 0.0
    k = 0
    for k in range(max_iter):
        # 確認用の画面出力
        print(str(k) + " : f(" + fmt(x1) + ") = " + fmt(f1))

        # 次の位置を計算
        if abs(f1-f0) < 1e-10:  # ゼロ除算防止
            break
        x2 = x1 - f1*(x1-x0)/(f1-f0)

        # 関数値の計算(次のループのため＆収束判定用)
        f0 = f1
        f1 = func(x2)

        # 収束判定
        dx = abs(x2-x1)
        if dx < eps or abs(f1) < eps:
            break

        x0 = x1
        x1 = x2

    return 0, x0, x1, k, dx


def Pivoting(A, n, k):
    """!
    ピボット選択(Pivoting)
     - 行入れ替えだけの部分的ピボッティング
    @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列
    @param[in] n n元連立一次方程式
    @param[in] k 対象行
    """
    # k行目以降でk列目の絶対値が最も大きい要素を持つ行を検索
    p = k                # 絶対値が最大の行
    am = abs(A[k][k])    # 最大値
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


class FUNCTION:
    """!
    関数とその導関数(勾配)をまとめたもの
     - C++版の struct FUNCTION に相当
    """
    def __init__(self, func, dfunc):
        self.func = func
        self.dfunc = dfunc


def newton_nd(funcs, x, n, max_iter, eps):
    """!
    ニュートン・ラフソン法(n次元版)
     - C++版ではオーバーロードで同じ newton という名前にしていたが，
       Pythonには関数のオーバーロードがないので名前を変えている
    @param[in] funcs 関数群(FUNCTIONのリスト)
    @param[inout] x 解(探索開始位置を渡す，リストなので破壊的に更新される)
    @param[in] n 次元数
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    f = zeros(n)
    J = zeros2(n, n+1)  # 最後の列に右辺項ベクトルを入れるためにn+1にしている

    d = 0.0
    k = 0
    for k in range(max_iter):
        # 現在の位置xにおける関数値の計算
        for i in range(n):
            f[i] = -funcs[i].func(x)

        # ヤコビ行列の計算
        for i in range(n):
            df = funcs[i].dfunc(x)
            for j in range(n):
                J[i][j] = df[j]

        # 確認用の画面出力
        s = str(k) + " : "
        for j in range(n):
            s += "f" + str(j) + "("
            for i in range(n):
                s += fmt(x[i]) + (") = " if i == n-1 else ",")
            s += fmt(f[j]) + ("" if j == n-1 else ", ")
        print(s)

        # 線形システムを解いてδを計算
        for i in range(n):
            J[i][n] = f[i]  # 上の関数値の計算のところでマイナスをすでに付けてあることに注意
        GaussEliminationWithPivoting(J, n)

        # xを更新
        for i in range(n):
            x[i] += J[i][n]

        # 収束判定
        d = 0.0
        for i in range(n):
            d += abs(f[i])
        if d < eps:
            break

    return 0, k, d


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # # 探索開始位置
    # x = -1
    #
    # # ニュートン法でf(x)=0を解く
    # max_iter = 100
    # eps = 1e-6
    # ret, x, max_iter, eps = newton(Func1, DFunc1, x, max_iter, eps)
    # # x0 = x - 0.1
    # # ret, x0, x, max_iter, eps = secant(Func1, x0, x, max_iter, eps)
    #
    # # 結果の画面表示
    # set_precision(12)
    # print("x = " + fmt(x))
    # print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    # print("")

    # 多次元のニュートン法
    xv = [1.0, 0.0]
    funcs = []
    funcs.append(FUNCTION(Func4, DFunc4))
    funcs.append(FUNCTION(Func4a, DFunc4a))

    max_iter = 100
    eps = 1e-6
    ret, max_iter, eps = newton_nd(funcs, xv, 2, max_iter, eps)

    # 結果の画面表示
    print("(x1,x2) = (" + fmt(xv[0]) + ", " + fmt(xv[1]) + ")")
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    return 0


if __name__ == "__main__":
    main()
