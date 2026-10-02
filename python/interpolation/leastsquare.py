# -*- coding: utf-8 -*-
"""!
@file leastsquare.py

@brief 最小二乗近似

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


# -----------------------------------------------------------------------------
# LU分解(最小自乗法での係数計算に用いる)
# -----------------------------------------------------------------------------
def LUDecomp(A, n):
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
            A[i][j] = lu/A[i][i]

    return 0


def LUSolver(A, b, x, n):
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
        x[i] = bly/A[i][i]

    # 後退代入(back substitution)
    #  UX=YからXを計算
    for i in range(n-1, -1, -1):
        yux = x[i]
        for j in range(i+1, n):
            yux -= A[i][j]*x[j]
        x[i] = yux

    return 0


# -----------------------------------------------------------------------------
# 補間法
# -----------------------------------------------------------------------------
def leastsquare_interpolation(fi, xi, m, n, x, c):
    """!
    最小二乗法による補間
     - 1次元データ(xi,fi)に対して，1次元n次多項式をフィッティング
    @param[in] fi 関数値を格納した配列
    @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
    @param[in] m データ数
    @param[in] n フィッティングする多項式の次数
    @param[in] x 補間した値が必要な位置x
    @param[out] c 係数ベクトル(リストなので破壊的に更新される)
    @return 補間値
    """
    # 多項式の次数から行列のサイズ(=基底ベクトルの次元数)を計算
    dim_b = n+1  # 1次元の場合

    # Ac=y, bは基底ベクトル
    del c[:]
    c.extend(zeros(dim_b))  # 結果の係数ベクトル
    b = zeros(dim_b)        # 基底ベクトル
    y = zeros(dim_b)        # 右辺項
    A = zeros2(dim_b, dim_b)

    # 係数行列Aと右辺項yの計算
    for k in range(m):  # データ分だけ反復
        # 多項式の基底ベクトルの計算(多項式以外でフィッティングする場合はこの部分を変えれば良い)
        xb = 1
        for i in range(dim_b):
            b[i] = xb
            xb *= xi[k]

        # 係数行列Aと右辺項bの計算
        for i in range(dim_b):
            for j in range(dim_b):
                A[i][j] += b[i]*b[j]
            y[i] += b[i]*fi[k]

    # Ac=yをLU分解で解く(疎行列とは限らないのでCGソルバは使わない)
    LUDecomp(A, dim_b)
    LUSolver(A, y, c, dim_b)

    # 基底ベクトルに求めた係数を掛けて行くことで位置xにおける値yを計算
    fx = 0
    xb = 1
    for i in range(dim_b):
        fx += c[i]*xb
        xb *= x
    return fx


def polynominal1d(x, c, n):
    """!
    1次元n次多項式を計算
    @param[in] x 値が必要な位置
    @param[in] c 多項式の係数
    @param[in] n 多項式の次数
    @return 多項式を計算した結果
    """
    fx = 0
    xb = 1
    for i in range(n+1):
        fx += c[i]*xb
        xb *= x
    return fx


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

    # 補間用の点(サンプリング点の数=6)
    xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/5, func)
    OutputSamplingPoints(xi, yi)  # サンプリング点の画面出力

    # 3次多項式による補間
    c = []  # 係数ベクトル
    fx = leastsquare_interpolation(yi, xi, len(xi), 3, x, c)
    print("f_ls3(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-gt)))
    print("ground truth = " + fmt(gt))

    # # グラフ描画用にデータ出力
    # dat = os.path.join(SCRIPT_DIR, "dat")
    # if not os.path.exists(dat): os.makedirs(dat)
    # d = 3   # 多項式の次数
    # m = 11  # サンプリング点数
    # xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func)
    # leastsquare_interpolation(yi, xi, len(xi), d, x, c)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "leastsquare"+TOSTR(m)+"_data.txt"))  # サンプリング点のファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, len(c)-1),
    #                os.path.join(dat, "leastsquare"+TOSTR(m)+"_d"+TOSTR(d)+".txt"))
    # # チェビシェフ節点を用いる場合
    # xi, yi = MakeChebyshevNodes(x0, x1, (x1-x0)/(m-1.0), func)
    # leastsquare_interpolation(yi, xi, len(xi), d, x, c)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "leastsquare"+TOSTR(m)+"c_data.txt"))
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, len(c)-1),
    #                os.path.join(dat, "leastsquare"+TOSTR(m)+"_d"+TOSTR(d)+"c.txt"))
    #
    # # 真値のグラフ作成用ファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, func, os.path.join(dat, "leastsquare_ground_truth.txt"))
    #
    # # ノイズ付きデータに対する最小２乗法
    # xi, yi = MakeSamplingPointsWithWhiteNoise(x0, x1, (x1-x0)/(m-1.0), func, 0.15)
    # leastsquare_interpolation(yi, xi, len(xi), d, x, c)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "leastsquare"+TOSTR(m)+"n_data.txt"))  # サンプリング点のファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, len(c)-1),
    #                os.path.join(dat, "leastsquare"+TOSTR(m)+"n_d"+TOSTR(d)+".txt"))
    #
    # # 外れ値付きデータに対する最小２乗法
    # xi, yi = MakeSamplingPointsWithWhiteNoise(x0, x1, (x1-x0)/(m-1.0), func, 0.01)
    # yi[2] = func(xi[2])+1.2  # 外れ値を意図的に追加
    # leastsquare_interpolation(yi, xi, len(xi), d, x, c)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "leastsquare"+TOSTR(m)+"o_data.txt"))  # サンプリング点のファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, lambda t: polynominal1d(t, c, len(c)-1),
    #                os.path.join(dat, "leastsquare"+TOSTR(m)+"o_d"+TOSTR(d)+".txt"))

    return 0


if __name__ == "__main__":
    main()
