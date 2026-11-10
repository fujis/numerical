# -*- coding: utf-8 -*-
"""!
@file penalty.py

@brief ペナルティ関数法による制約付き最小化問題
        - 内点ペナルティ関数法(バリア関数法) -> コード中ではLogBarrierFunc, InverseBarrierFunc
        - 外点ペナルティ関数法 -> コード中ではPenaltyFunc

@author Makoto Fujisawa
@date 2026-11 (Python版)
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

SQRT5 = 2.23606797749979


# -----------------------------------------------------------------------------
# 制約付き最小化問題における制約を与える関数の定義
#  - 等式制約の場合は h(x) = 0
#  - 不等式制約の場合は h(x) <= 0
#  - h(x) >= 0 の不等式制約を掛けたい場合は返値にマイナスを付けて，-h(x) <= 0とする
# -----------------------------------------------------------------------------
def Constraint1(x):
    """!
    制約条件h(x)とその勾配∇h(x) : h(x) = x+y-1 (等式制約 h(x) = 0 として使用)
    @param[in] x 変数
    @return 値or勾配
    """
    return x[0]+x[1]-1.0


def DConstraint1(x):
    dg = zeros(len(x))
    dg[0] = 1.0
    dg[1] = 1.0
    return dg


def Constraint2(x):
    """!
    制約条件h(x)とその勾配∇h(x) : h(x) = x+y-0.5 (不等式制約 h(x) <= 0 として使用)
    @param[in] x 変数
    @return 値or勾配
    """
    return x[0]+x[1]-0.5


def DConstraint2(x):
    dg = zeros(len(x))
    dg[0] = 1.0
    dg[1] = 1.0
    return dg


def Constraint3(x):
    """!
    制約条件h(x)とその勾配∇h(x) : h(x) = -x (不等式制約 h(x) <= 0 すなわち x >= 0 として使用)
    @param[in] x 変数
    @return 値or勾配
    """
    return -x[0]


def DConstraint3(x):
    dg = zeros(len(x))
    dg[0] = -1.0
    dg[1] = 0.0
    return dg


class ObjectiveAndConstraints:
    """!
    目的関数と制約をまとめたクラス(C++版の構造体ObjectiveAndConstraintsに対応)
    """
    def __init__(self):
        self.objectivefunc = None           # 目的関数の値を返す関数
        self.objectivegrad = None           # 目的関数の勾配値を返す関数
        self.constraintfuncs = []           # 制約関数の値を返す関数
        self.constraintgrads = []           # 制約関数の勾配値を返す関数
        self.constraintequalities = []      # 制約が等式(1)か不等式(0)を格納する配列
        self.num_constraint = 0             # 制約の数

    def AddConstraint(self, func, grad, equality):
        """!
        制約の追加
        @param[in] func,grad 制約関数とその勾配を返す関数
        @param[in] equality 1:等式制約 h(x) = 0, 0:不等式制約 h(x) <= 0
        """
        self.constraintfuncs.append(func)
        self.constraintgrads.append(grad)
        self.constraintequalities.append(equality)
        self.num_constraint = len(self.constraintfuncs)


def Violation(oac, x):
    """!
    制約の違反量 Σ|h_i(x)| (等式制約) + Σmax(0,h_j(x)) (不等式制約)
    @param[in] oac 目的関数と制約関数
    @param[in] x 変数
    """
    v = 0.0
    for j in range(oac.num_constraint):
        h = oac.constraintfuncs[j](x)
        if oac.constraintequalities[j] == 1:
            v += math.fabs(h)   # 等式制約
        elif h > 0.0:
            v += h              # 不等式制約
    return v


# -----------------------------------------------------------------------------
# 内点ペナルティ関数法 : バリア関数を加えた目的関数 F(x) = f(x)+b(x)
#  - LogBarrierFuncやDLogBarrierFuncを既存の最適化問題の関数に引数xだけの関数として渡すために，
#    バリアパラメータや制約関数などはグローバル変数で与える形にしている
#    (外点ペナルティ関数法のPenaltyFunc,DPenaltyFuncも同様)
#  - g_oacはbarrierfunction()やpenaltyfunction()の中で引数のoacから設定する
# -----------------------------------------------------------------------------
g_mu = 1.0                          # バリアパラメータμ
g_oac = ObjectiveAndConstraints()   # 目的関数と制約関数


def LogBarrierFunc(x):
    """!
    ログバリア関数を加えた目的関数
     - ログバリア関数 : b(x) = -μ Σ log(-h_j(x))
     - 制約範囲外(h_j(x) >= 0)では大きな値を返す(探索点が範囲外に出ないように)
    @param[in] x 変数
    """
    # 目的関数の値
    f = g_oac.objectivefunc(x)

    # バリア関数の値を目的関数の値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        if h >= 0.0:
            return 1e30  # ゼロ割やlog(0)を防ぐためのエラー処理
        f -= g_mu*math.log(-h)
    return f


def DLogBarrierFunc(x):
    """!
    ログバリア関数を加えた目的関数の勾配値
     - ログバリア関数の勾配:　∇b(x) = -μ Σ ∇h_j(x)/h_j(x)
    @param[in] x 変数
    """
    # 目的関数の勾配値
    df = list(g_oac.objectivegrad(x))

    # バリア関数の勾配値を目的関数の勾配値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        dh = g_oac.constraintgrads[j](x)
        c = -g_mu/h
        for k in range(len(df)):
            df[k] += c*dh[k]
    return df


def InverseBarrierFunc(x):
    """!
    逆数バリア関数を加えた目的関数
     - 逆数バリア関数 : b(x) = μ Σ 1/h_j(x)^2
     - 制約範囲外(h_j(x) >= 0)では大きな値を返す(探索点が範囲外に出ないように)
    @param[in] x 変数
    """
    # 目的関数の値
    f = g_oac.objectivefunc(x)

    # バリア関数の値を目的関数の値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        if h >= 0.0:
            return 1e30  # ゼロ割を防ぐためのエラー処理
        f += g_mu/(h*h)
    return f


def DInverseBarrierFunc(x):
    """!
    逆数バリア関数を加えた目的関数の勾配値
     - 逆数バリア関数の勾配:　∇b(x) = -2μ Σ ∇h_j(x)/h_j(x)^3
    @param[in] x 変数
    """
    # 目的関数の勾配値
    df = list(g_oac.objectivegrad(x))

    # バリア関数の勾配値を目的関数の勾配値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        dh = g_oac.constraintgrads[j](x)
        c = -2.0*g_mu/(h*h*h)
        for k in range(len(df)):
            df[k] += c*dh[k]
    return df


# -----------------------------------------------------------------------------
# 外点ペナルティ関数法 : ペナルティ関数を加えた目的関数
#  F(x) = f(x)+ρ(Σ|h_i(x)|^α + Σ(max(0,h_j(x)))^β)  (h_i:等式制約, h_j:不等式制約)
#  - α = β = 2 とする(微分可能にするため．α = β = 1 では制約の境界で微分できない)
# -----------------------------------------------------------------------------
g_rho = 1.0     # ペナルティパラメータρ


def PenaltyFunc(x):
    """!
    ペナルティ関数を加えた目的関数
     - 制約範囲内では元の目的関数f(x)と同じ
    @param[in] x 変数
    """
    # 目的関数の値
    f = g_oac.objectivefunc(x)

    # ペナルティ関数の値を目的関数の値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        if g_oac.constraintequalities[j] == 1:     # 等式制約
            f += g_rho*h*h
        elif g_oac.constraintequalities[j] == 0:   # 不等式制約
            if h > 0.0:
                f += g_rho*h*h
    return f


def DPenaltyFunc(x):
    """!
    ペナルティ関数を加えた目的関数の勾配
     - ∇F(x) = ∇f(x) + 2ρ Σ h_i(x)∇h_i(x) + 2ρ Σ max(0,h_j(x))∇h_j(x)
    @param[in] x 変数
    """
    # 目的関数の勾配値
    df = list(g_oac.objectivegrad(x))

    # ペナルティ関数の勾配値を目的関数の勾配値に加える
    for j in range(g_oac.num_constraint):
        h = g_oac.constraintfuncs[j](x)
        if g_oac.constraintequalities[j] == 1:     # 等式制約
            dh = g_oac.constraintgrads[j](x)
            for k in range(len(df)):
                df[k] += 2.0*g_rho*h*dh[k]
        elif g_oac.constraintequalities[j] == 0:   # 不等式制約
            if h <= 0.0:
                continue
            dh = g_oac.constraintgrads[j](x)
            for k in range(len(df)):
                df[k] += 2.0*g_rho*h*dh[k]
    return df


# -----------------------------------------------------------------------------
# 制約なしの最適化関数(最小値探索) : quasinewton.py と同じ
# -----------------------------------------------------------------------------
def goldensection(func, x0, d, xl, xr, max_iter, eps):
    """!
    黄金分割探索法(golden-section method)
     - 多変数で探索方向が決まっている場合
     - BFGSでの探索方向に対する係数α決定のために用いる
    @param[in] func 関数値を与える関数
    @param[in] x0 探索初期地点ベクトル
    @param[in] d 探索方向ベクトル
    @param[in] xl,xr 初期探索範囲
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解ans, 実際の反復数, 実際の誤差)
    """
    eta = 2.0/(SQRT5+3.0)
    beta = xr-xl
    tau = eta*beta
    n = len(x0)

    x = zeros(4)
    x[0] = xl
    x[1] = xl+tau
    x[2] = xr-tau
    x[3] = xr

    # 探索方向上の2点x1,x2での関数値を計算
    x1 = zeros(n)
    x2 = zeros(n)
    for i in range(n):
        x1[i] = x0[i]+x[1]*d[i]
        x2[i] = x0[i]+x[2]*d[i]
    f1 = func(x1)
    f2 = func(x2)

    k = 0
    for k in range(max_iter):
        if f1 < f2:  # 区間[x2,x3]に最小値なし
            # 分割区間更新
            tx2 = x[2]
            x[2] = x[1]
            x[1] = tx2-tau
            x[3] = tx2

            # 関数値の更新
            for i in range(n):
                x1[i] = x0[i]+x[1]*d[i]
            f2 = f1
            f1 = func(x1)
        else:  # 区間[x0,x1]に最小値なし
            # 分割区間更新
            tx1 = x[1]
            x[1] = x[2]
            x[2] = tx1+tau
            x[0] = tx1

            # 関数値の更新
            for i in range(n):
                x2[i] = x0[i]+x[2]*d[i]
            f1 = f2
            f2 = func(x2)

        beta -= tau      # 新しい[x0,x3]の長さ(現在の長さから[x0,x1](or[x2,x3])の長さtauを引く)
        tau = eta*beta   # 新しい[x0,x1](or[x2,x3])の長さ

        if beta < eps:
            break

    return 0, x[1], k, beta


def quasinewton_bfgs(func, dfunc, x0, max_iter, eps):
    """!
    BFGS法(準ニュートン法の一種, Broyden-Fletcher-Goldfarb-Shanno method)
     - quasinewton.pyとほぼ同じ(各反復での画面出力を削除し，直線探索の範囲を制約範囲内に制限)
    @param[in] func 関数値を与える関数
    @param[in] dfunc 関数勾配値を与える関数
    @param[in] x0 初期探索地点
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 解xout, 実際の反復数, 実際の誤差)
    """
    n = len(x0)
    H = zeros2(n, n)  # ヘッセ行列の逆行列H^-1
    for i in range(n):
        H[i][i] = 1.0

    x = list(x0)
    g = dfunc(x)
    gnorm = 1.0

    k = 0
    for k in range(max_iter):
        # 探索方向dの計算(d = -H∇f)
        p = zeros(n)
        for i in range(n):
            for j in range(n):
                p[i] -= H[i][j]*g[j]

        # 探索範囲[0,tmax]の端が制約範囲外にならないようにtmaxを小さくする(バリア関数法のため)
        # (範囲外ではfuncが大きな値を返すので，両端が範囲外だと黄金分割法がうまく動かない)
        tmax = 3.0
        xt = zeros(n)
        for l in range(60):
            for i in range(n):
                xt[i] = x[i]+tmax*p[i]
            if func(xt) < 1e30:
                break
            tmax *= 0.5

        # 黄金分割法で探索方向d上の最小値までの距離を調べてalphaとして設定
        gmax_iter = 40
        geps = 1e-6*tmax
        ret, alpha, gmax_iter, geps = goldensection(func, x, p, 0.0, tmax, gmax_iter, geps)

        # xの更新
        s = zeros(n)
        for i in range(n):
            s[i] = alpha*p[i]
            x[i] += s[i]

        # yの計算
        gprev = g
        g = dfunc(x)
        y = zeros(n)
        for i in range(n):
            y[i] = g[i]-gprev[i]

        gnorm = 0.0
        for i in range(n):
            gnorm += g[i]*g[i]
        if math.sqrt(gnorm) < eps:
            break  # 勾配ベクトル∇fの大きさで収束判定

        # H*yの計算(行列×縦ベクトル⇒縦ベクトル)
        Hy = zeros(n)
        for i in range(n):
            for j in range(n):
                Hy[i] += H[i][j]*y[j]

        # (s^T)*yと(y^T)*(H*y)の計算(どちらもベクトル同士の内積)
        sy = 0.0
        yHy = 0.0
        for i in range(n):
            sy += s[i]*y[i]
            yHy += y[i]*Hy[i]

        # H^-1の更新
        for i in range(n):
            for j in range(n):
                H[i][j] += s[i]*s[j]/sy + yHy*s[i]*s[j]/(sy*sy) - Hy[i]*s[j]/sy - s[i]*Hy[j]/sy

    return 0, x, k, math.sqrt(gnorm)


# -----------------------------------------------------------------------------
# 制約付きの最適化関数
# -----------------------------------------------------------------------------
def barrierfunction(oac, x0, mu0, r, max_iter, eps, inverse=False):
    """!
    内点ペナルティ関数法(バリア関数法, barrier function method)
     1. 制約範囲内の内点を初期点x^(0)とし，バリアパラメータの初期値μ^(0)>0を設定
     2. k = 1～k_maxまで以下を反復
        a. f(x)+b(x)を目的関数として最小解x^(k)を求める(準ニュートン法．初期点はx^(k-1))
        b. μ^(k) < ε ならばx^(k)を解として反復終了
        c. μ^(k+1) = r μ^(k) (0<r<1)
    @param[in] oac 目的関数と制約関数(不等式制約のみ)
    @param[in] x0 初期探索地点(制約範囲の内点)
    @param[in] mu0 バリアパラメータの初期値
    @param[in] r バリアパラメータを小さくする割合(0<r<1)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @param[in] inverse Falseならログバリア関数，Trueなら逆数バリア関数を用いる
    @return (戻り値(1:成功,0:失敗), 解xout, 実際の反復数, 最後のμ)
    """
    global g_mu, g_oac

    # 等式制約があるかをチェック(等式制約では制約範囲の内点がないので，内点ペナルティ関数法(バリア関数法)は使えない)
    num_equalities = 0
    for i in range(oac.num_constraint):
        num_equalities += (1 if oac.constraintequalities[i] == 1 else 0)
    if num_equalities > 0:
        print("barrier function method cannot handle equality constraints")
        return 0, list(x0), 0, eps
    # 初期点が内点かどうかの確認
    for i in range(oac.num_constraint):
        h = oac.constraintfuncs[i](x0)
        if h >= 0.0:
            print("initial point is not an interior point (h" + str(i) + " = " + fmt(h) + ")")
            return 0, list(x0), 0, eps

    # BFGSに渡す関数(LogBarrierFuncなど)で使う目的関数と制約関数を設定
    g_oac = oac

    # 使用するバリア関数
    func = InverseBarrierFunc if inverse else LogBarrierFunc
    dfunc = DInverseBarrierFunc if inverse else DLogBarrierFunc

    x = list(x0)
    g_mu = mu0

    k = 1
    for k in range(1, max_iter+1):
        # f(x)+b(x)を準ニュートン法で最小化(初期点は前回の解)
        qmax_iter = 100
        qeps = 1e-6
        ret, x, qmax_iter, qeps = quasinewton_bfgs(func, dfunc, x, qmax_iter, qeps)

        print("k = " + str(k) + ", mu = " + fmt(g_mu) + " : x = (" + fmt(x) + "), f = " + fmt(oac.objectivefunc(x))
              + " (BFGS iter = " + str(qmax_iter) + ")")

        if g_mu < eps:
            break       # バリアパラメータが十分小さくなったら終了
        g_mu *= r       # バリアパラメータを小さくする

    return 1, x, k, g_mu


def penaltyfunction(oac, x0, rho0, r, max_iter, eps):
    """!
    外点ペナルティ関数法(exterior penalty function method)
     1. 初期点x^(0)(制約範囲外でもよい)とペナルティパラメータの初期値ρ^(0)>0を設定
     2. k = 1～k_maxまで以下を反復
        a. f(x)+ρ(...)を目的関数として最小解x^(k)を求める(準ニュートン法．初期点はx^(k-1))
        b. 制約の違反量が ε 未満ならばx^(k)を解として反復終了
        c. ρ^(k+1) = r ρ^(k) (r>1)
    @param[in] oac 目的関数と制約関数(等式制約，不等式制約のどちらでもよい)
    @param[in] x0 初期探索地点
    @param[in] rho0 ペナルティパラメータの初期値
    @param[in] r ペナルティパラメータを大きくする割合(r>1)
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値(1:成功,0:失敗), 解xout, 実際の反復数, 最後の制約の違反量)
    """
    global g_rho, g_oac

    # BFGSに渡す関数(PenaltyFunc,DPenaltyFunc)で使う目的関数と制約関数を設定
    g_oac = oac

    x = list(x0)
    g_rho = rho0

    v = Violation(oac, x)
    k = 1
    for k in range(1, max_iter+1):
        # f(x)+ペナルティ関数を準ニュートン法で最小化(初期点は前回の解)
        qmax_iter = 100
        qeps = 1e-6
        ret, x, qmax_iter, qeps = quasinewton_bfgs(PenaltyFunc, DPenaltyFunc, x, qmax_iter, qeps)

        v = Violation(oac, x)
        print("k = " + str(k) + ", rho = " + fmt(g_rho) + " : x = (" + fmt(x) + "), f = " + fmt(oac.objectivefunc(x))
              + ", violation = " + fmt(v) + " (BFGS iter = " + str(qmax_iter) + ")")

        if v < eps:
            break       # 制約の違反量が十分小さくなったら終了
        g_rho *= r      # ペナルティパラメータを大きくする

    return 1, x, k, v


def print_result(oac, x, iter):
    """! 最終計算結果の表示 """
    print("  (x,y) = (" + fmt(x[0]) + "," + fmt(x[1]) + "), f(x,y) = " + fmt(oac.objectivefunc(x)) + ", iter = " + str(iter))
    print("")


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    set_precision(8)

    # 問題1 : 不等式制約
    print("argmin x^2+y^2-xy-x-1  s.t. x+y-0.5 <= 0, -x <= 0")
    print("exact : (x,y) = (" + fmt(5.0/12.0) + "," + fmt(1.0/12.0) + "), f(x,y) = " + fmt(-183.0/144.0))
    print("")

    oac = ObjectiveAndConstraints()
    oac.objectivefunc = Func4
    oac.objectivegrad = DFunc4
    oac.AddConstraint(Constraint2, DConstraint2, 0)  # x+y-0.5 <= 0
    oac.AddConstraint(Constraint3, DConstraint3, 0)  # -x <= 0

    # 内点ペナルティ関数法(バリア関数法) : ログバリア関数
    print("[barrier function method (interior penalty) : log barrier]")
    x0 = zeros(2, 0.1)  # 初期値(制約範囲の内点)
    ret, x, max_iter, eps = barrierfunction(oac, x0, 1.0, 0.1, 100, 1e-6, False)
    print_result(oac, x, max_iter)

    # 内点ペナルティ関数法(バリア関数法) : 逆数バリア関数
    print("[barrier function method (interior penalty) : inverse barrier]")
    x0 = zeros(2, 0.1)  # 初期値(制約範囲の内点)
    ret, x, max_iter, eps = barrierfunction(oac, x0, 1.0, 0.1, 100, 1e-6, True)
    print_result(oac, x, max_iter)

    # 外点ペナルティ関数法
    print("[exterior penalty method]")
    x0 = zeros(2, 1.0)  # 初期値(制約範囲外でもよい)
    ret, x, max_iter, eps = penaltyfunction(oac, x0, 1.0, 10.0, 100, 1e-6)
    print_result(oac, x, max_iter)

    # 問題2 : 等式制約(外点ペナルティ関数法のみ)
    print("argmin x^2+2y^2  s.t. x+y-1 = 0")
    print("exact : (x,y) = (" + fmt(2.0/3.0) + "," + fmt(1.0/3.0) + "), f(x,y) = " + fmt(2.0/3.0))
    print("")

    oac = ObjectiveAndConstraints()
    oac.objectivefunc = Func6
    oac.objectivegrad = DFunc6
    oac.AddConstraint(Constraint1, DConstraint1, 1)  # x+y-1 = 0

    # 外点ペナルティ関数法
    print("[exterior penalty method]")
    x0 = zeros(2, 0.0)  # 初期値
    ret, x, max_iter, eps = penaltyfunction(oac, x0, 1.0, 10.0, 100, 1e-6)
    print_result(oac, x, max_iter)

    return 0


if __name__ == "__main__":
    main()
