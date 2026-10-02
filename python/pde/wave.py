# -*- coding: utf-8 -*-
"""!
@file wave.py

@brief 波動方程式のソルバー

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


# -----------------------------------------------------------------------------
# 波動方程式のソルバー
# -----------------------------------------------------------------------------
def setbc_d(f, n):
    """!
    境界条件設定(ノイマン境界条件)
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    """
    f[0] = f[1]
    f[n] = f[n-1]


def wave1d_ftcs(f2, f1, f0, c, dt, n, x0, xn):
    """!
    FTCS法(前進オイラー法+中心差分)による波動方程式のソルバー
    @param[inout] f2 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f1 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップi-1での値)
    @param[in] c 波の速度を表す係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n < 0:
        return 1
    dx = (xn-x0)/n  # 空間刻み幅
    eta = c*c*dt*dt/(dx*dx)
    for i in range(1, n):
        f2[i] = 2*f1[i]-f0[i]+eta*(f1[i+1]-2*f1[i]+f1[i-1])

    return 0


def wave1d_step(sdfunc, bcfunc, f2, f1, c, dt, n, x0, xn):
    """!
    1ステップ進める
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f2,f1 未知関数fの各グリッドでの値
    @param[in] c 波の速度を表す係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    """
    f0 = zeros(len(f1))
    for i in range(len(f2)):
        f0[i] = f1[i]
        f1[i] = f2[i]
    bcfunc(f0, n)
    bcfunc(f1, n)
    sdfunc(f2, f1, f0, c, dt, n, x0, xn)
    bcfunc(f2, n)


def wave1d(sdfunc, bcfunc, f, c, dt, n, x0, xn, filename):
    """!
    時間ステップを進めていったときの結果をファイル出力するための関数
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピーして使う)
    @param[in] c 波の速度を表す係数
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @param[in] filename 出力ファイル名
    """
    f = list(f)  # C++版では値渡しなので，Pythonでもここでコピーを作る
    fo = open(filename, "w")
    fo.write("#1d,-0.4,0.8\n")

    # 初期値の出力
    OutputValueToFile(f, n, x0, xn, 0, fo)

    f0 = list(f)

    # タイムステップを進める
    for k in range(1, 4000):
        wave1d_step(sdfunc, bcfunc, f, f0, c, dt, n, x0, xn)

        # 結果の出力
        OutputValueToFile(f, n, x0, xn, k*dt, fo)
    fo.close()

    print("wave1d : " + filename)


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    dt = 0.01
    x0, xn = 0.0, 1.0
    c = 0.02
    n = 100
    bcfunc = setbc_d

    f = zeros(n+1)

    # 初期値設定
    for i in range(n+1):
        f[i] = 0.0
    f[n//2] = 0.5

    bcfunc(f, n)

    # データ保存パス(ファイル名の最後の_ignはGitで除外ファイルにするためのもの)
    path = os.path.join(SCRIPT_DIR, "data")
    if not os.path.exists(path):
        os.makedirs(path)

    # 前進オイラー+中心差分
    wave1d(wave1d_ftcs, bcfunc, f, c, dt, n, x0, xn, os.path.join(path, "wave1d_ftcs_ign.txt"))

    return 0


if __name__ == "__main__":
    main()
