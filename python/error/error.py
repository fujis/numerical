# -*- coding: utf-8 -*-
"""!
@file error.py

@brief 数値の精度と誤差

@author Makoto Fujisawa
@date 2019-09 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math
import struct

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *


# -----------------------------------------------------------------------------
# 単精度(float)を模擬するための関数
#  - Pythonのfloatは常にC++のdouble(倍精度,64bit)相当なので，
#    単精度(32bit)の挙動を見るにはstructで32bitに丸めてやる必要がある
# -----------------------------------------------------------------------------
def to_float32(x):
    """!
    倍精度の値を単精度(32bit)に丸める(C++のfloat型への代入に相当)
    @param[in] x 実数値
    @return 単精度に丸めた値
    """
    return struct.unpack("f", struct.pack("f", x))[0]


# -----------------------------------------------------------------------------
# メイン関数
# -----------------------------------------------------------------------------
def main():
    # 整数の範囲
    #  - Pythonのint型は桁数に制限がない(メモリの許す限りいくらでも大きくできる)ので，
    #    ここではC++での各型の範囲を参考値として表示している
    print("char : {0} - {1}  ({2}bit)".format(-128, 127, 8))
    print("short : {0} - {1}  ({2}bit)".format(-32768, 32767, 16))       # 正確にはshort int
    print("int : {0} - {1}  ({2}bit)".format(-2147483648, 2147483647, 32))
    print("long long : {0} - {1}  ({2}bit)".format(-2**63, 2**63-1, 64))
    print("float : {0} - {1}  ({2}bit)".format(fmt(1.17549435e-38), fmt(3.40282347e+38), 32))
    print("double : {0} - {1}  ({2}bit)".format(fmt(sys.float_info.min), fmt(sys.float_info.max), 64))
    print("python int : 桁数の制限なし (必要なだけメモリを使う)")
    print("")

    # float型の有効桁数
    x = to_float32(123456789)
    y = to_float32(123456700)
    z = to_float32(x-y)
    set_precision(10)  # 表示桁数を10桁にする
    print("x = " + fmt(x))
    print("y = " + fmt(y))
    print("x-y = " + fmt(z))
    print("")

    # 丸め誤差
    # print("\n[round-off error]")
    a = to_float32(1.1)
    set_precision(20)  # 表示する桁数を20桁に設定
    print("a = " + fmt(a))  # 1.1ぴったりになるはずだけど...
    print("")

    # 丸め誤差による計算への影響
    print("\n[effect by round-off error]")
    b = to_float32(a*3)
    if b == 3.3:
        print("b is equal to 3.3")
    else:
        print("b is not equal to 3.3")
    print("")

    # 演算で発生する誤差
    print("\n[effect by round-off error 2]")
    a1 = to_float32(1234567)
    a2 = to_float32(0.00123)
    a3 = to_float32(1234567)
    print("(a1+a2)-a3=" + fmt(to_float32(to_float32(a1+a2)-a3)))
    print("(a1-a3)+a2=" + fmt(to_float32(to_float32(a1-a3)+a2)))

    # # 桁落ち誤差を生じないように計算する例
    # t1 = 0.0
    # for j in range(100000000):
    #     t1 += 0.1
    #
    # t = [0.0]*10000
    # for j in range(10000):
    #     t[j] = 0.0
    #     for i in range(10000):
    #         t[j] += 0.1
    # t1 = 0.0
    # for j in range(10000):
    #     t1 += t[j]
    # print(t1)

    # 倍精度と単精度の違い
    print("\n[single/double precision]")
    x1 = to_float32(1.1)  # 単精度
    x2 = 1.1              # 倍精度(Pythonのfloatはこちら)
    print("x1 = " + fmt(x1) + "  (single precision - 32bit)")
    print("x2 = " + fmt(x2) + "  (double precision - 64bit)")

    # 打ち切り誤差
    # ライプニッツの公式でπを計算
    print("\n[trancation error]")
    n = 100
    pi = 0
    sgn = 1
    for i in range(n+1):
        pi += sgn/(2.0*i+1.0)
        sgn *= -1
        if i % (n//10) == 0:
            print("n=" + str(i) + " : " + fmt(4*pi))
    pi0 = 3.141592653589793  # 真値
    set_precision(20)        # 表示する桁数を20桁に設定
    print("error = " + fmt(4*pi-pi0))

    # 真値が分からない場合の収束判定
    print("\n[convergence test]")
    eps = 1.0e-5  # 精度を更に上げる場合は桁落ち誤差に注意
    pi = 0
    sgn = 1
    m = 0
    for i in range(1001):
        pi0 = pi   # 収束判定のために前の反復の値を確保しておく
        pi += sgn/(2.0*i+1.0)*4.0
        sgn *= -1
        if abs(pi-pi0) <= eps:
            break  # 収束判定
        m += 1
    print("pi = " + fmt(pi))
    set_precision(20)  # 表示する桁数を20桁に設定
    print("iterations : " + str(m))

    return 0


if __name__ == "__main__":
    main()
