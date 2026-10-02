# -*- coding: utf-8 -*-
"""!
@file simplex.py

@brief シンプレックス法

@author Makoto Fujisawa
@date 2019-07 (Python版)
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
# シンプレックス法(線形計画問題の解法)
# -----------------------------------------------------------------------------
def simplex(a, eqn, n_cond, m, ans, max_iter):
    """!
    シンプレックス法(simplex method)
     - 線形計画問題の解法
     - 二段階シンプレックスではないので原点が解に含まれない場合は適応できない
    @param[in] a 係数リスト(条件式の係数を格納した2次元配列,最後の1行は最適化式)
    @param[in] eqn 各式の等号不等号(1:<=, -1:>=)
    @param[in] n_cond 条件式の数
    @param[in] m 変数の数
    @param[out] ans 解(リストなので破壊的に更新される)
    @param[in] max_iter 最大反復数
    @return (戻り値, 実際の反復数)
    """
    n = n_cond+1        # 式の数(条件式の数+最適化式の数)
    m_slack = n_cond    # スラック変数の数
    m_all = m+m_slack   # スラック変数を含む全変数の数

    # 単体表の作成
    #  単体表の行数は条件式の数+1=n，最後の行が目的関数
    #  単体表の列数はスラック変数を含む全変数の数+1，最後の列が基底可能解
    s = zeros2(n, m_all+1)  # 単体表
    xi = [0]*n              # それぞれの行の基底変数(x1=0,x2=1,x3=2,...とインデックス値を格納)
    for i in range(n):
        if i != n-1:  # 条件式の行(0～n-2行)
            xi[i] = i+m  # 初期基底変数はスラック変数
            sgn = -1 if a[i][m] < 0 else 1  # 条件式の右辺項の符号
            for j in range(m):
                s[i][j] = sgn*a[i][j]       # 変数の係数を入れていく(右辺項の符号が正になるようにしている)
            for j in range(m, m_all):
                # スラック変数部分の初期係数は単位行列のような形になる(ただし条件式の符号が>=の場合は-1をセット)
                s[i][j] = sgn*eqn[i] if xi[i] == j else 0
            s[i][m_all] = sgn*a[i][m]       # 右辺項を初期基底可能解としてセット
        else:  # 最終行は最適化式(ここだけ別で設定)
            xi[i] = -1  # 最適化式の変数インデックスには-1を格納しておく
            for j in range(m):
                s[i][j] = -a[i][j]  # 最終行の最適化式は係数の符号を反転
            for j in range(m, m_all+1):
                s[i][j] = 0         # 最後の行の最適化式では初期基底可能解に0をセット

    # チェック用
    for i in range(n):
        line = " z: " if i == n-1 else "x" + str(xi[i]+1) + ": "
        for j in range(m_all+1):
            line += fmt(s[i][j]).rjust(8) + " "
        print(line)
    print("")

    k = 0
    for k in range(max_iter):
        # 非基底変数(初期状態ではスラック変数以外)から負で絶対値最大のものを選択
        b = -1
        xmax = 0.0
        for j in range(m_all):
            if s[n-1][j] < 0 and -s[n-1][j] > xmax:
                b = j
                xmax = -s[n-1][j]
        if b == -1:
            break  # 負の変数が見つからなければ収束したとしてループを抜ける

        # ピボット要素の探索
        p = 0
        ymin = s[0][m_all]/s[0][b]
        for i in range(1, n-1):
            # 基底可能解をb列の値で割った値が最小のものを選択
            y = s[i][m_all]/s[i][b]
            if y < ymin:
                p = i
                ymin = y

        # ピボット要素行を選択された非基底変数(b)で置き換えて，ピボット要素でその行を割る
        xi[p] = b
        xp = s[p][b]  # ループ内で値が変わるので一時変数に値を確保しておく
        for j in range(m_all+1):
            s[p][j] /= xp

        # ピボット行(p)以外の行を選択された非基底変数(b)の値が0になるようにピボット行の係数をx倍して引く
        #  -> s[p][b]は1になっているのでs[p][j]にs[i][b]を掛けて引けば良い
        for i in range(n):
            if i == p:
                continue  # ピボット行は飛ばす
            xp = s[i][b]  # ループ内で値が変わるので一時変数に値を確保しておく
            for j in range(m_all+1):
                s[i][j] -= s[p][j]*xp

    # 解の格納
    for i in range(n-1):
        if xi[i] < m:
            ans[xi[i]] = s[i][m_all]
    ans[m] = s[n-1][m_all]

    # チェック用
    for i in range(n):
        line = " z: " if i == n-1 else "x" + str(xi[i]+1) + ": "
        for j in range(m_all+1):
            line += fmt(s[i][j]).rjust(8) + " "
        print(line)

    return 0, k


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 条件式と目的関数の係数
    # a0 = [[3, 1, 9], [2.5, 2, 12.5], [1, 2, 8], [3, 2, 0]]
    # a0 = [[1, 2, 800], [3, 4, 1800], [3, 1, 1500], [20, 30, 0]]
    a0 = [[2, 1, 8], [1, 3, 9], [1, 1, 0]]
    # eqn0 = [1, 1, 1, 0]
    eqn0 = [1, 1, 0]  # 条件式の符号(1:<=, -1:>=, 0:=), 最後の要素は目的関数

    n_cond = 2  # 条件式の数
    m = 2       # 変数の数
    x = zeros(m+1)
    a = zeros2(n_cond+1, m+1)
    eqn = [0]*(n_cond+1)

    for i in range(n_cond+1):
        eqn[i] = eqn0[i]
        for j in range(m+1):
            a[i][j] = a0[i][j]

    # シンプレックス法で線形計画問題を解く
    max_iter = 100
    ret, max_iter = simplex(a, eqn, n_cond, m, x, max_iter)

    # 結果の画面表示
    s = "f("
    for i in range(m):
        s += fmt(x[i]) + (") = " if i == m-1 else ",")
    s += fmt(x[m]) + ", "
    s += "iter = " + str(max_iter)
    print(s)
    print("")

    return 0


if __name__ == "__main__":
    main()
