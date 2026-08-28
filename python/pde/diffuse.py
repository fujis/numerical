# -*- coding: utf-8 -*-
"""!
@file diffuse.py

@brief 拡散方程式のソルバー

@author Makoto Fujisawa
@date 2019-11 (Python版)
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

# このスクリプトが置かれているフォルダ(データファイルの出力先の基準にする)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 境界値
alpha = 0.0
beta = 0.0
cn_lambda = 0.5


# -----------------------------------------------------------------------------
# 拡散方程式のソルバー
# -----------------------------------------------------------------------------
def cg_solver(A, b, x, n, max_iter, eps):
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
        a = rr0/dot(p, y, n)

        # 解x、残差rの更新
        for i in range(n):
            x[i] += a*p[i]
            r[i] -= a*y[i]

        # (r*r)_(k+1)の計算
        rr1 = dot(r, r, n)

        # 収束判定 (||r||<=eps)
        e = math.sqrt(rr1)
        if e < eps:
            k += 1
            break

        # βの計算とPの更新
        bt = rr1/rr0
        for i in range(n):
            p[i] = r[i]+bt*p[i]

        # (r*r)_(k+1)を次のステップのために確保しておく
        rr0 = rr1

    return 0, k+1, e


def setbc_d(f, n):
    """!
    境界条件設定(ディリクレ境界条件, 1D)
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    """
    f[0] = alpha
    f[n] = beta


def setbc2_d(f, n):
    """!
    境界条件設定(ディリクレ境界条件, 2D)
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    """
    for i in range(n):
        f[i][0] = alpha
        f[i][n] = beta
    for j in range(n):
        f[0][j] = alpha
        f[n][j] = beta


def diffuse1d_ftcs(f1, f0, a, dt, n, x0, xn):
    """!
    FTCS法(前進オイラー法+中心差分)による拡散方程式のソルバー
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n < 0:
        return 1
    dx = (xn-x0)/n  # 空間刻み幅
    eta = a*dt/(dx*dx)
    for i in range(1, n):
        f1[i] = f0[i]+eta*(f0[i+1]-2*f0[i]+f0[i-1])

    return 0


def diffuse1d_cn(f1, f0, a, dt, n, x0, xn):
    """!
    クランク・ニコルソン法による拡散方程式のソルバー
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n < 0:
        return 1
    dx = (xn-x0)/n  # 空間刻み幅
    eta = a*dt/(dx*dx)

    # 線形システムの係数行列と右辺項ベクトルの計算
    b = zeros(n-1)
    x = zeros(n-1)
    A = zeros2(n-1, n-1)  # n-1×n-1の2次元配列
    for i in range(n-1):
        j = i+1  # 行列としてはi行目でf[i+1]での式となる
        if i != 0:
            A[i][i-1] = -cn_lambda*eta
        A[i][i] = 1+2*cn_lambda*eta
        if i != n-2:
            A[i][i+1] = -cn_lambda*eta
        b[i] = f0[j]+(1-cn_lambda)*eta*(f0[j+1]-2*f0[j]+f0[j-1])
    b[0] += cn_lambda*eta*f0[0]
    b[n-2] += cn_lambda*eta*f0[n]
    for i in range(n-1):
        x[i] = f0[i+1]

    # CG法で線形システムを解く
    max_iter = 100
    eps = 1e-6
    cg_solver(A, b, x, n-1, max_iter, eps)
    # print("cg solver : max_iter = " + str(max_iter) + ", eps = " + fmt(eps))

    # 結果を配列fに戻す
    for i in range(n-1):
        f1[i+1] = x[i]

    return 0


def diffuse2d_ftcs(f1, f0, a, dt, n, x0, xn, y0, yn):
    """!
    FTCS法(前進オイラー法+中心差分)による拡散方程式のソルバー(2D)
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
    @param[in] x0,xn,y0,yn 計算範囲
    @return 問題なければ0を返す
    """
    if n < 0:
        return 1
    dx = (xn-x0)/n  # 空間刻み幅
    dy = (yn-y0)/n  # 空間刻み幅

    eta1 = a*dt/(dx*dx)
    eta2 = a*dt/(dy*dy)
    for i in range(1, n):
        for j in range(1, n):
            f1[i][j] = f0[i][j]+eta1*(f0[i+1][j]-2*f0[i][j]+f0[i-1][j])+eta2*(f0[i][j+1]-2*f0[i][j]+f0[i][j-1])

    return 0


def diffuse2d_cn(f1, f0, a, dt, n, x0, xn, y0, yn):
    """!
    クランク・ニコルソン法による拡散方程式のソルバー(2D)
     - 係数行列が(n-1)^2×(n-1)^2の密行列になるので，Pythonでは計算にかなり時間がかかる
       (nを10程度に下げて試すことを推奨)
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
    @param[in] x0,xn,y0,yn 計算範囲
    @return 問題なければ0を返す
    """
    if n < 0:
        return 1
    dx = (xn-x0)/n  # 空間刻み幅
    dy = (yn-y0)/n  # 空間刻み幅

    eta1 = a*dt/(dx*dx)
    eta2 = a*dt/(dy*dy)

    # 線形システムの係数行列と右辺項ベクトルの計算
    m = (n-1)*(n-1)
    b = zeros(m)
    x = zeros(m)
    A = zeros2(m, m)  # (n-1)^2×(n-1)^2の2次元配列
    for i in range(n-1):
        i1 = i+1
        for j in range(n-1):
            j1 = j+1
            k = i+j*(n-1)

            A[k][k] = 1+2*cn_lambda*eta1+2*cn_lambda*eta2

            if i != 0:
                A[k][(i-1)+(j)*(n-1)] = -cn_lambda*eta1
            if i != n-2:
                A[k][(i+1)+(j)*(n-1)] = -cn_lambda*eta1
            if j != 0:
                A[k][(i)+(j-1)*(n-1)] = -cn_lambda*eta2
            if j != n-2:
                A[k][(i)+(j+1)*(n-1)] = -cn_lambda*eta2

            b[k] = (f0[i1][j1]
                    + (1-cn_lambda)*eta1*(f0[i1+1][j1]-2*f0[i1][j1]+f0[i1-1][j1])
                    + (1-cn_lambda)*eta2*(f0[i1][j1+1]-2*f0[i1][j1]+f0[i1][j1-1]))
    for j in range(n-1):
        b[(0)+j*(n-1)] += cn_lambda*eta1*f0[0][j+1]
        b[(n-2)+j*(n-1)] += cn_lambda*eta1*f0[n][j+1]
    for i in range(n-1):
        b[i+(0)*(n-1)] += cn_lambda*eta2*f0[i+1][0]
        b[i+(n-2)*(n-1)] += cn_lambda*eta2*f0[i+1][n]
    for i in range(n-1):
        for j in range(n-1):
            x[i+j*(n-1)] = f0[i+1][j+1]

    # CG法で線形システムを解く
    max_iter = 100
    eps = 1e-6
    cg_solver(A, b, x, m, max_iter, eps)
    # print("cg solver : max_iter = " + str(max_iter) + ", eps = " + fmt(eps))

    # 結果を配列fに戻す
    for i in range(n-1):
        for j in range(n-1):
            f1[i+1][j+1] = x[i+j*(n-1)]

    return 0


def diffuse1d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn):
    """!
    1ステップ進める
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f 未知関数fの各グリッドでの値
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    """
    f0 = list(f)
    bcfunc(f0, n)
    sdfunc(f, f0, a, dt, n, x0, xn)
    bcfunc(f, n)


def diffuse1d(sdfunc, bcfunc, f, a, dt, n, x0, xn, filename, nstep=1000):
    """!
    時間ステップを進めていったときの結果をファイル出力するための関数
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピーして使う)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @param[in] filename 出力ファイル名
    @param[in] nstep 計算する時間ステップ数
    """
    f = list(f)  # C++版では値渡しなので，Pythonでもここでコピーを作る
    fo = open(filename, "w")
    fo.write("#1d,0.0,1.0\n")

    # 初期値の出力
    OutputValueToFile(f, n, x0, xn, 0, fo)

    # タイムステップを進める
    for k in range(1, nstep):
        diffuse1d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn)

        # 結果の出力
        OutputValueToFile(f, n, x0, xn, k*dt, fo)
    fo.close()

    print("diffuse1d : " + filename)


def diffuse2d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn, y0, yn):
    """!
    1ステップ進める(2D)
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f 未知関数fの各グリッドでの値
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
    @param[in] x0,xn,y0,yn x方向,y方向の計算範囲
    """
    f0 = [list(row) for row in f]
    bcfunc(f0, n)
    sdfunc(f, f0, a, dt, n, x0, xn, y0, yn)
    bcfunc(f, n)


def diffuse2d(sdfunc, bcfunc, f, a, dt, n, x0, xn, y0, yn, filename, nstep=200):
    """!
    時間ステップを進めていったときの結果をファイル出力するための関数(2D)
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピーして使う)
    @param[in] a 拡散係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
    @param[in] x0,xn,y0,yn x方向,y方向の計算範囲
    @param[in] filename 出力ファイル名
    @param[in] nstep 計算する時間ステップ数
    """
    f = [list(row) for row in f]  # C++版では値渡しなので，Pythonでもここでコピーを作る
    fo = open(filename, "w")
    fo.write("#2d,0.0,1.0\n")

    # 初期値の出力
    OutputValueToFile2D(f, n, x0, xn, y0, yn, 0, fo)

    # タイムステップを進める
    for k in range(1, nstep):
        diffuse2d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn, y0, yn)

        # 結果の出力
        OutputValueToFile2D(f, n, x0, xn, y0, yn, k*dt, fo)
    fo.close()

    print("diffuse2d : " + filename)


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    global cn_lambda

    dt = 0.01
    x0, xn = 0.0, 1.0
    a = 0.02
    n = 50
    bcfunc = setbc_d

    f = zeros(n+1)

    # 初期値設定
    for i in range(n+1):
        f[i] = 1.0
    bcfunc(f, n)

    # データ保存パス(ファイル名の最後の_ignはGitで除外ファイルにするためのもの)
    path = os.path.join(SCRIPT_DIR, "data")
    if not os.path.exists(path):
        os.makedirs(path)

    def fn(name):
        return os.path.join(path, name)

    # 前進オイラー+中心差分
    diffuse1d(diffuse1d_ftcs, bcfunc, f, a, dt, n, x0, xn, fn("diffuse1d_ftcs_ign.txt"))
    # クランク・ニコルソン法
    cn_lambda = 0.5
    diffuse1d(diffuse1d_cn, bcfunc, f, a, dt, n, x0, xn, fn("diffuse1d_cn_ign.txt"))

    dt = 0.011
    # 前進オイラー+中心差分
    diffuse1d(diffuse1d_ftcs, bcfunc, f, a, dt, n, x0, xn, fn("diffuse1d_ftcs_dt11_ign.txt"))
    # クランク・ニコルソン法
    cn_lambda = 0.5
    diffuse1d(diffuse1d_cn, bcfunc, f, a, dt, n, x0, xn, fn("diffuse1d_cn_dt11_ign.txt"))
    # クランク・ニコルソン法(dt=0.05)
    dt = 0.05
    diffuse1d(diffuse1d_cn, bcfunc, f, a, dt, n, x0, xn, fn("diffuse1d_cn_dt50_ign.txt"))

    # 2次元
    #  - クランク・ニコルソン法(diffuse2d_cn)は係数行列が(n-1)^2×(n-1)^2の密行列になるため，
    #    Pythonで実行するとn=50では非常に時間がかかる．そのため2次元計算ではnを小さくしている．
    n2 = 20
    y0, yn = 0.0, 1.0
    dt = 0.005
    bcfunc2 = setbc2_d
    f2 = zeros2(n2+1, n2+1, 1.0)
    bcfunc2(f2, n2)

    # 前進オイラー+中心差分
    diffuse2d(diffuse2d_ftcs, bcfunc2, f2, a, dt, n2, x0, xn, y0, yn, fn("diffuse2d_ftcs_ign.txt"))
    # クランク・ニコルソン法
    cn_lambda = 0.5
    diffuse2d(diffuse2d_cn, bcfunc2, f2, a, dt, n2, x0, xn, y0, yn, fn("diffuse2d_cn_ign.txt"))

    return 0


if __name__ == "__main__":
    main()
