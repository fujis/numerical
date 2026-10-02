# -*- coding: utf-8 -*-
"""!
@file dka.py

@brief DKA法(ワイヤストラス法+アバースの初期値)
       代数方程式(多項式からなる方程式)の求根問題

@author Makoto Fujisawa
@date 2012-06 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math
import cmath  # 複素数(Pythonでは複素数はcomplex型として標準で使える)

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *
from rx_funcs import *


def cstr(x):
    """!
    複素数の表示用文字列を作る
     - C++版でcomplex型に対して<<オペレータをオーバーロードしていたものに相当
    @param[in] x 複素数
    @return 表示用文字列
    """
    s = fmt(x.real)
    if abs(x.imag) > 1e-6:
        s += " + " + fmt(x.imag) + "i"
    return s


# -----------------------------------------------------------------------------
# 代数方程式
# -----------------------------------------------------------------------------
def func(x, b, n):
    """!
    係数から代数方程式の値を計算
     - C++版ではテンプレート関数．実数でも複素数でもそのまま使える
    @param[in] x 変数
    @param[in] b 係数
    @param[in] n 方程式の次数
    @return 代数方程式の値
    """
    f = b[n]
    for i in range(n):
        f += b[i]*x**(n-i)
    return f


def func_h(x, b, n):
    """!
    ホーナー法で代数方程式の値を計算
    @param[in] x 変数
    @param[in] b 係数
    @param[in] n 方程式の次数
    @return 代数方程式の値
    """
    f = b[0]
    for i in range(1, n+1):
        f = b[i]+f*x
    return f


def dfunc_h(x, b, n):
    """!
    ホーナー法で代数方程式の導関数値を計算
    @param[in] x 変数
    @param[in] b 係数
    @param[in] n 方程式の次数
    @return 代数方程式の導関数値
    """
    df = n*b[0]
    for i in range(1, n):
        df = (n-i)*b[i]+df*x
    return df


# -----------------------------------------------------------------------------
# DKA法による求根問題の解法
# -----------------------------------------------------------------------------
def newton(b, n, x, max_iter, eps):
    """!
    ニュートン・ラフソン法(Newton-Raphson method)
     - Aberthの初期値を求めるために用いる
    @param[in] b 多項式の係数
    @param[in] n 方程式の次数
    @param[in] x 探索開始位置
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解x, 実際の反復数, 実際の誤差)
    """
    dx = 0.0
    k = 0
    for k in range(max_iter):
        # 現在の位置xにおける関数値と導関数の計算
        f = func_h(x, b, n)
        df = dfunc_h(x, b, n)

        if abs(df) < 1e-16:  # ゼロ除算防止(C++版では暗黙にinfになる)
            break

        # 導関数の結果から次の位置を計算
        x = x-f/df

        # 収束判定
        dx = abs(f/df)
        if dx < eps or abs(f) < eps:
            break

    return 0, x, k, dx


def horner(a, b, n):
    """!
    ホーナー法(組立除法)
     - P(x) = a0 x^n + a1 x^(n-1) + ... + a_(n-1) x + a_n を (x-b) で割ったときの商と余りを返す
     - 商はn-1次の多項式の係数として返す
    @param[in] a 代数方程式の係数
    @param[in] b 割る1次式の係数(x-b)
    @param[in] n 方程式の次数(配列aの大きさはn+1,配列cはn)
    @return (c, rm) 商であるn-1次の多項式の係数とあまり
    """
    if n <= 1:
        return None, 0.0
    rm = a[0]  # 最終的に余りになる
    c = zeros(n)
    for i in range(1, n+1):
        c[i-1] = rm
        rm *= b
        rm += a[i]
    return c, rm


def aberth(c, n, max_iter, eps):
    """!
    Aberthの方法で初期値を算出
    @param[in] c 多項式の係数(複素数のリスト)
    @param[in] n 方程式の次数
    @param[in] max_iter 半径計算のための最大反復数
    @param[in] eps 半径計算のための許容誤差
    @return 解探索のための初期値z(複素数のリスト)
    """
    # 半径算出のための方程式の係数
    cd = [1.0+0j]*(n+1)  # 係数c'
    c1n = -c[1].real/n
    cd[0] = c[0]
    # 係数格納用の一時的な変数
    a = [c[i].real for i in range(n+1)]
    # zの多項式をz+c1/nで割っていくことでwの多項式の係数を求める
    for i in range(n, 1, -1):
        tmp, rm = horner(a, c1n, i)
        cd[i] = rm
        a = tmp
    cd[1] = a[1]+c1n

    # 多項式S(w)の係数
    b = zeros(len(cd))
    b[0] = abs(cd[0])
    for i in range(1, n+1):
        b[i] = -abs(cd[i])

    # Aberthの初期値の半径をニュートン法で算出
    m = 0  # 係数cの中で0でないものの個数
    for i in range(n+1):
        m += 1 if abs(c[i].real) > 1e-6 else 0
    rmax = 0.0  # 半径の最大値
    for i in range(1, n+1):
        ri = pow(m*abs(c[i].real)/abs(c[0].real), 1.0/(i+1.0))
        if ri > rmax:
            rmax = ri
    print("r_max = " + fmt(rmax))
    r = rmax
    ret, r, max_iter, eps = newton(b, n, r, max_iter, eps)
    print("r = " + fmt(r))
    r = 10

    # Aberthの初期値
    zc = -c[1]/(c[0]*n)
    z = [0j]*n
    s = ""
    for j in range(n):
        theta = (2*RX_PI/n)*j+RX_PI/(2.0*n)
        z[j] = zc+r*complex(math.cos(theta), math.sin(theta))
        s += "x" + str(j) + "(0) = " + cstr(z[j]) + ("" if j == n-1 else ", ")
    print(s)

    return z


def weierstrass(z, c, n, max_iter, eps):
    """!
    ワイヤストラス法(DK公式)
    @param[inout] z 初期値位置を受け取り，解を返す(リストなので破壊的に更新される)
    @param[in] c 多項式の係数
    @param[in] n 方程式の次数
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    e = 0.0

    s = "0, "
    for i in range(n):
        s += fmt(z[i].real) + ", " + fmt(z[i].imag) + ", "
    print(s)

    for k in range(max_iter):
        zp = list(z)

        # DK式の計算
        for j in range(n):
            f = func(z[j], c, n)
            df = c[0]
            for i in range(n):
                if i != j:
                    df *= zp[j]-zp[i]

            z[j] = zp[j]-f/df

        s = str(k+1) + ", "
        for i in range(n):
            s += fmt(z[i].real) + ", " + fmt(z[i].imag) + ("" if i == n-1 else ", ")
        print(s)

        # 誤差の算出
        e = 0.0
        for j in range(n):
            ej = abs(func(z[j], c, n))
            if ej > e:
                e = ej

        # 収束判定
        if e < eps:
            return 1, k, e

    return 0, max_iter, e


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 代数方程式の係数を読み込む
    #  - データファイルはこのスクリプトと同じフォルダに置いてある
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "algebra.txt")
    c0 = ReadAlgebra(fn, ",")
    if c0:
        n = len(c0)-1  # 代数方程式の次数
        c = [complex(c0[i]) for i in range(n+1)]
    else:
        print("file read error!")
        return 1

    # (確認用)元の方程式の表示
    s = "f(x) = "
    for i in range(n):
        s += fmt(c0[i]) + "x^" + str(n-i) + " + "
    s += fmt(c0[n]) + " = 0"
    print(s)
    print("")

    max_iter = 100
    eps = 1e-6

    # Aberthの初期値
    z = aberth(c, n, max_iter, eps)

    # ワイヤストラス法(DK公式)で解を求める
    ret, max_iter, eps = weierstrass(z, c, n, max_iter, eps)

    # 結果の画面表示
    print("solutions : ")
    for i in range(n):
        s = "x" + str(i) + " = " + cstr(z[i])

        # 解の精度のチェック(f(x)を計算してみる)
        f = complex(c0[n])
        zi = z[i]
        for j in range(n-1, -1, -1):
            f += c0[j]*zi
            zi *= z[i]
        s += "  --> f = " + cstr(f)
        if math.sqrt(f.real*f.real+f.imag*f.imag) < 1e-6:
            s += " (OK)"
        print(s)
    print("")
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))

    return 0


if __name__ == "__main__":
    main()
