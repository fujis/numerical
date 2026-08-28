# -*- coding: utf-8 -*-
"""!
@file sor.py

@brief 逐次加速緩和法(SOR法:Successive Over-Relaxation)

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


# -----------------------------------------------------------------------------
# 緩和法による線形システムソルバ
# -----------------------------------------------------------------------------
def SOR(A, n, w, max_iter, eps):
    """!
    逐次加速緩和法(SOR法:Successive Over-Relaxation)
     - ガウスザイデル反復法に加速係数をかけたもの
    @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
    @param[in] n n元連立一次方程式
    @param[in] w 加速緩和乗数 ([0,2])
    @param[in] max_iter 最大反復数
    @param[in] eps 許容誤差
    @return (戻り値, 実際の反復数, 実際の誤差)
    """
    x = zeros(n)  # 初期値はすべて0とする
    set_precision(8)

    e = 0.0  # 誤差
    k = 0    # 計算反復回数
    for k in range(max_iter):
        # 現在の値を代入して，次の解候補を計算
        l = 0
        e = 0.0
        for i in range(n):
            tmp = x[i]
            x[i] = A[i][n]
            for j in range(n):
                x[i] -= (A[i][j]*x[j] if j != i else 0.0)
            x[i] /= A[i][i]

            x[i] = (1.0-w)*tmp+w*x[i]  # 加速緩和係数wを使って次の解を計算

            # if abs((tmp-x[i])/tmp) > eps:  # 相対誤差の場合
            if abs(tmp-x[i]) > eps:  # 絶対誤差の場合
                e += abs(tmp-x[i])
                l += 1

        # 確認のため現在の解を画面出力
        s = str(k) + " : "
        for i in range(n):
            s += "x" + str(i) + " = " + fmt(x[i]) + ("" if i == n-1 else ", ")
        print(s)

        # 収束判定
        if l == 0:  # すべての解が許容誤差以下なら反復終了
            break

    for i in range(n):
        A[i][n] = x[i]

    return 0, k, e


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # ファイルから行列要素を読み込んで解く
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.txt")
    A = ReadMatrix(fn, ",")

    n = len(A)     # n元連立一次方程式
    m = len(A[0])

    # 繰り返し計算のために右辺項ベクトルを退避しておく
    b = zeros(n)
    for i in range(n):
        b[i] = A[i][n]

    # 読み込んだ行列を確認用に画面表示
    print("A(" + str(n) + " x " + str(m) + ") = ")
    OutputMatrix(A, n, n+1)
    print("")

    # 加速係数の入力
    w = 1.0
    while w >= 0.0:
        try:
            line = input("relaxation factor (enter -1 to stop) : ")
        except EOFError:
            break
        if not IsNumeric(line):
            continue
        w = float(line)
        if w < 0.0:
            break

        # 退避しておいたbを元に戻す
        for i in range(n):
            A[i][n] = b[i]

        # SOR法で線形システムを解く
        max_iter = 100
        eps = 1e-6
        ret, max_iter, eps = SOR(A, n, w, max_iter, eps)

        # 結果の画面表示
        s = ""
        for i in range(n):
            s += "x" + str(i) + " = " + fmt(A[i][n]) + ("" if i == n-1 else ", ")
        print(s)

        print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
        print("")

    return 0


if __name__ == "__main__":
    main()
