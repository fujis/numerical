# -*- coding: utf-8 -*-
"""!
@file advection.py

@brief 移流方程式のソルバー

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
# 移流法のソルバー
# -----------------------------------------------------------------------------
def advect1d_central(f1, f0, u, dt, n, x0, xn):
    """!
    中心差分でdt分だけ移流を進める
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (xn-x0)/n  # 空間刻み幅
    nu = u*dt/h    # クーラン数
    for i in range(1, n):
        # ∂f/∂xの差分による近似
        gx = (f0[i+1]-f0[i-1])/2.0
        # fの値を更新
        f1[i] += -nu*gx
    return 0


def advect1d_upwind(f1, f0, u, dt, n, x0, xn):
    """!
    風上差分でdt分だけ移流を進める
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (xn-x0)/n  # 空間刻み幅
    nu = u*dt/h    # クーラン数
    for i in range(1, n):
        # ∂f/∂xの差分による近似
        if u > 0:
            gx = f0[i]-f0[i-1]
        else:
            gx = f0[i+1]-f0[i]

        # fの値を更新
        f1[i] += -nu*gx
    return 0


def advect1d_lw(f1, f0, u, dt, n, x0, xn):
    """!
    Lax-Wendroff法でdt分だけ移流を進める
     - P. Lax and B.Wendroff, "Systems of conservation laws", Commun. Pure Appl Math. 13, pp.217-237, 1960.
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (xn-x0)/n  # 空間刻み幅
    nu = u*dt/h    # クーラン数
    for i in range(1, n):
        # ∂f/∂xの差分による近似
        gx = (f0[i+1]-f0[i-1])/2.0-nu*(f0[i+1]-2*f0[i]+f0[i-1])/2.0
        # fの値を更新
        f1[i] += -nu*gx
    return 0


def advect1d_sl(f1, f0, u, dt, n, x0, xn):
    """!
    セミラグランジュ法でdt分だけ移流を進める
    @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
    @param[in] f0 未知関数fの各グリッドでの値(ステップiでの値)
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (xn-x0)/n  # 空間刻み幅
    for i in range(1, n):
        # 現在の位置からバックトレースした位置を求める
        x = x0+i*h
        x -= u*dt  # バックトレース

        # バックトレースした位置のグリッド情報計算
        ib = int((x-x0)/h)  # xが含まれるグリッド番号
        if ib < 0:
            ib = 0
        if ib >= n:
            ib = n-1
        dx = (x-ib*h)/h  # グリッド位置からの距離を[0,1]で正規化

        # fの値を線形補間で求める
        f1[i] = (1-dx)*f0[ib]+dx*f0[ib+1]
    return 0


def setbc_d(f, n):
    """!
    境界条件設定(ディリクレ境界条件)
    @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    """
    f[0] = alpha
    f[n] = beta


def advect1d_step_forwardeular(sdfunc, bcfunc, f, u, dt, n, x0, xn):
    """!
    前進オイラー法で1ステップ進める
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f 未知関数fの各グリッドでの値
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    """
    f0 = list(f)
    sdfunc(f, f0, u, dt, n, x0, xn)
    bcfunc(f, n)


def advect1d_step_heun(sdfunc, bcfunc, f, u, dt, n, x0, xn):
    """!
    ホイン法で1ステップ進める
     - 空間差分部分の計算切り替えのために(主にセミラグランジュ法のために)，
       f(k+1)=f(k)+(dt/2)(F(k)+F(k+1)) を f(k+1)=(f(k)+dt*F(k))/2+(f(k)+dt*F(k+1))/2として計算している
     - 手法比較のために無駄な処理やメモリを使っているので実際に計算する場合は元の式で計算した方が良い
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f 未知関数fの各グリッドでの値
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    """
    f1 = list(f)
    f2 = list(f)

    # 前進オイラー法で仮のf(k+1)を求める
    sdfunc(f1, f, u, dt, n, x0, xn)
    bcfunc(f, n)

    # f(k)+dt*F(t(k+1),f(k+1))を計算
    sdfunc(f2, f1, u, dt, n, x0, xn)
    bcfunc(f2, n)

    # f(k+1)=(f(k)+dt*F(k))/2+(f(k)+dt*F(k+1))/2でfを更新
    for i in range(1, n):
        f[i] = f1[i]/2+f2[i]/2
    bcfunc(f, n)


def advect1d_step_rk4(sdfunc, bcfunc, f, u, dt, n, x0, xn):
    """!
    RK4で1ステップ進める
     - 空間差分部分の計算切り替えのために(主にセミラグランジュ法のために)，
       f(i+1)=f(i)+(dt/6)(k1+k2+k3+k4) を f(i+1)=(fi+dt*k1)/6+(fi+dt*k2)/3+(fi+dt*k3)/3+(fi+dt*k4)/6 として計算している
     - 手法比較のために無駄な処理やメモリを使っているので実際にRK4で計算する場合は元の式で計算した方が良い
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[inout] f 未知関数fの各グリッドでの値
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    """
    k1 = list(f)
    k2 = list(f)
    k3 = list(f)
    k4 = list(f)
    fk = list(f)
    fk2 = list(f)
    fk3 = list(f)

    # fi+dt*k1の計算
    sdfunc(k1, f, u, dt, n, x0, xn)
    bcfunc(k1, n)

    # fi+(dt/2)*k1=fi+(dt/2)*Fの計算
    sdfunc(fk, f, u, dt/2, n, x0, xn)
    bcfunc(fk, n)
    # fi+dt*k2の計算
    sdfunc(k2, fk, u, dt, n, x0, xn)
    bcfunc(k2, n)

    # fi+(dt/2)*k2=fi+(dt/2)*F(t+dt/2,fi+(dt/2)*F)の計算
    sdfunc(fk2, fk, u, dt/2, n, x0, xn)
    bcfunc(fk2, n)
    # fi+dt*k3の計算
    sdfunc(k3, fk2, u, dt, n, x0, xn)
    bcfunc(k3, n)

    # fi+dt*k3=fi+dt*F(t+dt,fi+(dt/2)*F(t+dt/2,fi+(dt/2)*F))の計算
    sdfunc(fk3, fk2, u, dt, n, x0, xn)
    bcfunc(fk3, n)
    # fi+dt*k4の計算
    sdfunc(k4, fk3, u, dt, n, x0, xn)
    bcfunc(k4, n)

    # f(i+1)=(fi+dt*k1)/6+(fi+dt*k2)/3+(fi+dt*k3)/3+(fi+dt*k4)/6
    for i in range(1, n):
        f[i] = k1[i]/6+k2[i]/3+k3[i]/3+k4[i]/6
    bcfunc(f, n)


def advect1d(tdfunc, sdfunc, bcfunc, f, u, dt, n, x0, xn, filename):
    """!
    時間ステップを進めていったときの結果をファイル出力するための関数
    @param[in] tdfunc 時間差分を行う関数
    @param[in] sdfunc 空間差分を行う関数
    @param[in] bcfunc 境界条件設定用関数
    @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピーして使う)
    @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
    @param[in] dt 時間ステップ幅
    @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
    @param[in] x0,xn 計算範囲
    @param[in] filename 出力ファイル名
    """
    f = list(f)  # C++版では値渡しなので，Pythonでもここでコピーを作る
    fo = open(filename, "w")
    fo.write("#1d,0.0,1.0\n")

    # 初期値の出力
    OutputValueToFile(f, n, x0, xn, 0, fo)

    # タイムステップを進める
    for k in range(1, 750):
        tdfunc(sdfunc, bcfunc, f, u, dt, n, x0, xn)

        # 結果の出力
        OutputValueToFile(f, n, x0, xn, k*dt, fo)
    fo.close()

    print("advect1d : " + filename)


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    x0, xn = 0.0, 5.0  # 計算範囲(空間)
    n = 128            # 分割数(空間)
    u = 0.75           # 移流速度
    dx = (xn-x0)/n
    dt = 0.1*dx/u      # 時間ステップ幅(時間方向の刻み幅)

    f = zeros(n+1)

    # 初期値設定
    SetValueRectangle(f, n, x0, xn)
    # SetValueSin(f, n, x0, xn)

    # データ保存パス
    path = os.path.join(SCRIPT_DIR, "data")
    if not os.path.exists(path):
        os.makedirs(path)

    def fn(name):
        return os.path.join(path, name)

    # 前進オイラー+中心差分
    advect1d(advect1d_step_forwardeular, advect1d_central, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_fe+cen_ign.txt"))
    # # 前進オイラー+風上差分
    # advect1d(advect1d_step_forwardeular, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_fe+upwind_ign.txt"))
    # # 前進オイラー+Lax-Wendroff
    # advect1d(advect1d_step_forwardeular, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_fe+lw_ign.txt"))
    # # 前進オイラー+セミラグランジュ法
    # advect1d(advect1d_step_forwardeular, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_fe+sl_ign.txt"))
    # # ホイン+中心差分
    # advect1d(advect1d_step_heun, advect1d_central, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_heun+cen_ign.txt"))
    # # ホイン+風上差分
    # advect1d(advect1d_step_heun, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_heun+upwind_ign.txt"))
    # # ホイン+Lax-Wendroff
    # advect1d(advect1d_step_heun, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_heun+lw_ign.txt"))
    # # ホイン+セミラグランジュ法
    # advect1d(advect1d_step_heun, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_heun+sl_ign.txt"))
    # # RK4+中心差分
    # advect1d(advect1d_step_rk4, advect1d_central, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_rk4+cen_ign.txt"))
    # # RK4+風上差分
    # advect1d(advect1d_step_rk4, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_rk4+upwind_ign.txt"))
    # # RK4+Lax-Wendroff
    # advect1d(advect1d_step_rk4, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_rk4+lw_ign.txt"))
    # # RK4+セミラグランジュ法
    # advect1d(advect1d_step_rk4, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, fn("advect1d_rk4+sl_ign.txt"))

    return 0


if __name__ == "__main__":
    main()
