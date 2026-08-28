# -*- coding: utf-8 -*-
"""!
@file interpolation.py

@brief 補間法

@author Makoto Fujisawa
@date 2019-08 (Python版)
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
# 補間法
# -----------------------------------------------------------------------------
def linear_interpolation(f, a, b, x):
    """!
    線形補間
    @param[in] f 関数値(2つ)
    @param[in] a,b 関数値fに対応する位置x
    @param[in] x 補間した値が必要な位置x
    @return 補間値
    """
    l = b-a
    return ((x-a)*f[1]+(b-x)*f[0])/l


def lagrangian_interpolation(fi, xi, m, x):
    """!
    ラグランジュ補間
     - 線形補間は2点の時のラグランジュ補間
    @param[in] fi 関数値を格納した配列
    @param[in] xi 関数値に対応する位置を格納した配列
    @param[in] m データ数
    @param[in] x 補間した値が必要な位置x
    @return 位置xでの補間値
    """
    n = m-1
    Ln = 0.0
    for i in range(n+1):
        # 補間係数(Π(x-xj)/(xi-xj) (j!=i) の計算)
        l = 1
        for j in range(n+1):
            if j == i:
                continue
            l *= (x-xi[j])/(xi[i]-xi[j])
        # 補間値
        Ln += fi[i]*l
    return Ln


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

    # 補間用の点
    xi, yi = MakeSamplingPoints(x0, x1, x1-x0, func)
    OutputSamplingPoints(xi, yi)  # サンプリング点の画面出力

    # 線形補間
    fx = linear_interpolation(yi, xi[0], xi[1], x)
    print("f_linear(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-gt)))

    # 2点ラグランジュ補間(=線形補間)
    fx = lagrangian_interpolation(yi, xi, len(xi), x)
    print("f_lagrangian(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-gt)))
    print("")

    # 補間用の点を増やす(6点⇒5区間なので5次多項式で補間)
    xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/5.0, func)
    OutputSamplingPoints(xi, yi)  # サンプリング点の画面出力

    # 6点ラグランジュ補間(5次多項式補間)
    fx = lagrangian_interpolation(yi, xi, len(xi), x)
    print("f_lagrangian(" + fmt(x) + ") = " + fmt(fx) + ",  error = " + fmt(abs(fx-gt)))
    print("ground truth = " + fmt(gt))

    # グラフ描画用にデータ出力
    dat = os.path.join(SCRIPT_DIR, "dat")
    if not os.path.exists(dat):
        os.makedirs(dat)
    m = 20  # データ点数(次数はデータ点数-1)
    xi, yi = MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func)
    OutputSamplingPoints(xi, yi, os.path.join(dat, "lagrangian"+TOSTR(m)+"_data.txt"))
    # C++版の std::bind に相当するものとしてPythonではラムダ式を使う
    OutputFunction(x0, x1, (x1-x0)/200,
                   lambda t: lagrangian_interpolation(yi, xi, len(xi), t),
                   os.path.join(dat, "lagrangian"+TOSTR(m)+".txt"))

    # # チェビシェフ節点を用いる場合
    # xi, yi = MakeChebyshevNodes(x0, x1, (x1-x0)/(m-1.0), func)
    # OutputSamplingPoints(xi, yi, os.path.join(dat, "lagrangian"+TOSTR(m)+"c_data.txt"))
    # OutputFunction(x0, x1, (x1-x0)/200,
    #                lambda t: lagrangian_interpolation(yi, xi, len(xi), t),
    #                os.path.join(dat, "lagrangian"+TOSTR(m)+"c.txt"))

    # # 真値のグラフ作成用ファイル出力
    # OutputFunction(x0, x1, (x1-x0)/200, func, os.path.join(dat, "lagrangian_ground_truth.txt"))

    return 0


if __name__ == "__main__":
    main()
