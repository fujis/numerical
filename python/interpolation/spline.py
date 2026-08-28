# -*- coding: utf-8 -*-
"""!
@file spline.py

@brief スプライン補間

@author Makoto Fujisawa
@date 2019-09 (Python版)
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

# このスクリプトが置かれているフォルダ(データファイルの入出力先の基準にする)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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
        return 0, max_iter, eps

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
        alpha = rr0/dot(p, y, n)

        # 解x、残差rの更新
        for i in range(n):
            x[i] += alpha*p[i]
            r[i] -= alpha*y[i]

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

    return 1, k+1, e


# -----------------------------------------------------------------------------
# 補間法
# -----------------------------------------------------------------------------
def spline_interpolation(f, xi, m, x):
    """!
    3次スプライン補間
    @param[in] f 関数値を格納した配列
    @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
    @param[in] m データ数
    @param[in] x 補間した値が必要な位置x
    @return (ans, a, b, c, d) 補間値と各区間における補間係数
    """
    n = m-1  # 補間区間の数
    # 補間係数配列のメモリ確保
    a = zeros(n)
    b = zeros(n)
    c = zeros(n)
    d = zeros(n)
    h = zeros(n)  # 各補間区間の幅
    for i in range(n):
        h[i] = xi[i+1]-xi[i]

    # 位置x_iでの多項式の2階微分u_iの計算
    # 線形システムの係数行列Hと右辺項ベクトルbの計算
    H = zeros2(n-1, n-1)  # 係数行列を0で初期化
    tmp = zeros(n-1)
    v = zeros(n-1)
    for i in range(n-1):
        # 係数行列と右辺項ベクトルの要素の計算
        if i != 0:
            H[i][i-1] = h[i]
        H[i][i] = 2*(h[i]+h[i+1])
        if i != n-2:
            H[i][i+1] = h[i+1]
        v[i] = 6*((f[i+2]-f[i+1])/h[i+1]-(f[i+1]-f[i])/h[i])
    # CG法で線形システムを解く(Hは対称疎行列)
    max_iter = 100
    eps = 1e-6
    cg_solver(H, v, tmp, n-1, max_iter, eps)

    # u_0とu_nを追加
    u = zeros(n+1)
    u[0] = 0.0
    u[n] = 0.0
    for i in range(n-1):
        u[i+1] = tmp[i]

    # 係数a,b,c,dの計算
    for i in range(n):
        a[i] = (u[i+1]-u[i])/(6.0*h[i])
        b[i] = u[i]/2.0
        c[i] = (f[i+1]-f[i])/h[i]-h[i]*(2*u[i]+u[i+1])/6.0
        d[i] = f[i]

    # xが含まれる区間の探索(区間幅がすべて同じならint(x/h)で求められる)
    k = 0
    for i in range(n):
        if x >= xi[i] and x < xi[i+1]:
            k = i
            break

    # 計算済みの係数を使って3次多項式で補間値を計算
    dx = x-xi[k]
    ans = a[k]*dx*dx*dx + b[k]*dx*dx + c[k]*dx + d[k]

    return ans, a, b, c, d


def cubic_spline(x, xi, n, a, b, c, d):
    """!
    3次スプライン多項式で関数値を計算
    @param[in] x 補間した値が必要な位置
    @param[in] xi 各区間の境界位置
    @param[in] n 区間数
    @param[in] a,b,c,d 各区間における補間係数
    @return スプライン多項式を計算した結果
    """
    # xが含まれる区間の探索(区間幅がすべて同じならint(x/h)で求められる)
    k = 0
    for i in range(n):
        if x >= xi[i] and x < xi[i+1]:
            k = i
            break

    # 計算済みの係数を使って3次多項式で補間値を計算
    dx = x-xi[k]
    return a[k]*dx*dx*dx + b[k]*dx*dx + c[k]*dx + d[k]


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    func = FuncExp
    x0, x1 = 0, 1

    # func = FuncRunge
    # x0, x1 = -1, 1

    x = 0.5
    gt = func(x)  # 真値
    set_precision(6)

    # 補間用の点(サンプリング点数6)
    xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/5, func)
    OutputSamplingPoints(xi, yi)  # サンプリング点の画面出力

    # スプライン補間
    fx, a, b, c, d = spline_interpolation(yi, xi, len(xi), x)
    print("f_spline(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-gt)))
    print("ground truth = " + fmt(gt))

    # # グラフ描画用にデータ出力
    # dat = os.path.join(SCRIPT_DIR, "dat")
    # if not os.path.exists(dat): os.makedirs(dat)
    # m = 11  # データ点数
    # xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func)
    # fx, a, b, c, d = spline_interpolation(yi, xi, len(xi), x)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "spline"+TOSTR(m)+"_data.txt"))  # サンプリング点のファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: cubic_spline(t, xi, len(xi)-1, a, b, c, d),
    #                os.path.join(dat, "spline"+TOSTR(m)+".txt"))
    #
    # # 真値のグラフ作成用ファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, func, os.path.join(dat, "spline_ground_truth.txt"))

    return 0


if __name__ == "__main__":
    main()
