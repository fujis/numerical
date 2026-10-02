# -*- coding: utf-8 -*-
"""!
@file rx_funcs.py

@brief 数値計算テスト：方程式の記述(C++版 rx_funcs.h のPython移植)

@author Makoto Fujisawa
@date 2012-06 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import math

from rx_utils import RX_PI


# -----------------------------------------------------------------------------
# 関数定義
# -----------------------------------------------------------------------------

def Func1(x):
    """!
    1次元の代数方程式
    @param[in] x 変数
    @return 方程式の値
    """
    # return x*x-2
    return 2*x*x*x*x*x+5*x*x*x+3*x+1
    # return math.cos(0.5*x)


def DFunc1(x):
    # return 2*x
    return 10*x*x*x*x+15*x*x+3


def Func2(x):
    """!
    1次元の凸関数
    @param[in] x 変数
    @return 方程式の値
    """
    return x*x-2*x+2


def DFunc2(x):
    return 2*x-2


def Func3(x):
    """!
    1次元の非凸関数
    @param[in] x 変数
    @return 方程式の値
    """
    return x*math.cos(2*x)


def DFunc3(x):
    return math.cos(2*x)-2*x*math.sin(2*x)


def Func4(x):
    """!
    2次元の凸関数
    @param[in] x 変数(2要素のリスト)
    @return 方程式の値
    """
    return x[0]*x[0]+x[1]*x[1]-x[0]*x[1]-x[0]-1
    # return x[0]*x[0]-4*x[0]*x[1]+x[1]*x[1]


def DFunc4(x):
    g = [0.0]*len(x)
    g[0] = 2*x[0]-x[1]-1
    g[1] = 2*x[1]-x[0]
    # g[0] = 2*x[0]-4*x[1]
    # g[1] = -4*x[0]+2*x[1]
    return g


def Func4a(x):
    """!
    2次元の凸関数2
    @param[in] x 変数(2要素のリスト)
    @return 方程式の値
    """
    return x[0]*x[0]+x[1]*x[1]-2


def DFunc4a(x):
    g = [0.0]*len(x)
    g[0] = 2*x[0]
    g[1] = 2*x[1]
    return g


def Func5(x):
    """!
    2次元関数(鞍点が生じる例)
    @param[in] x 変数(2要素のリスト)
    @return 方程式の値
    """
    return 2*x[0]*x[0]-x[0]*x[1]-6*x[1]*x[1]+2*x[1]


def DFunc5(x):
    g = [0.0]*len(x)
    g[0] = 4*x[0]-x[1]
    g[1] = -x[0]-12*x[1]+2
    return g


def FuncExp(x):
    """!
    1次元の指数関数
    @param[in] x 変数
    @return 方程式の値
    """
    return math.exp(x)


def DFuncExp(x):
    return math.exp(x)


def FuncPi(x):
    """!
    円周率計算用関数
    @param[in] x 変数
    @return 方程式の値
    """
    return math.cos(x/2.0)


def DFuncPi(x):
    return -math.sin(x/2.0)/2.0


def FuncLinear(x):
    """!
    1次元の1次関数
    @param[in] x 変数
    @return 方程式の値
    """
    return 2.0*x+1.0


def DFuncLinear(x):
    return 2.0


# ルンゲ関数の係数
ainv = 25


def FuncRunge(x):
    """!
    1次元のルンゲ関数
    @param[in] x 変数
    @return 方程式の値
    """
    return 1.0/(1.0+ainv*x*x)


def DFuncRunge(x):
    return -2*ainv*x/((1.0+ainv*x*x)*(1.0+ainv*x*x))


def FuncT17(x):
    """!
    2017年度筑波大前期日程入試問題[4]の数式
    @param[in] x 変数
    @return 方程式の値
    """
    return 2*x*x-9*x+14-9/x+2/(x*x)


def DFuncT17(x):
    return 4*x-9+9/(x*x)-4/(x*x*x)


def FuncP2(x, y):
    """!
    2次元関数(積分用2次多項式)
    @param[in] x,y 変数
    @return 方程式の値
    """
    return 8*x*x+4*y


def FuncY1(x):
    return 2-x


def FuncY2(x):
    return x*x


# -----------------------------------------------------------------------------
# 円,球の面積,体積計算用
# -----------------------------------------------------------------------------
sr = 1.0  # 円,球の半径


def FuncCircleTop(x):
    """!
    円の上半分の形状
    @param[in] x 変数
    @return 方程式の値
    """
    y = sr*sr-x*x
    return 0.0 if y <= 0 else math.sqrt(y)


def DFuncCircleTop(x):
    y = sr*sr-x*x
    return 0.0 if y <= 1e-10 else -x/math.sqrt(y)


def FuncCircleBottom(x):
    """!
    円の下半分の形状
    @param[in] x 変数
    @return 方程式の値
    """
    y = sr*sr-x*x
    return 0.0 if y <= 0 else -math.sqrt(y)


def DFuncCircleBottom(x):
    y = sr*sr-x*x
    return 0.0 if y <= 1e-10 else x/math.sqrt(y)


def FuncSphere(x, y):
    """!
    2次元関数(球体の体積計算用)
    @param[in] x,y 変数
    @return 方程式の値
    """
    z = sr*sr-x*x-y*y
    return 0.0 if z <= 0 else math.sqrt(z)


def FuncSphereY1(x):
    z = sr*sr-x*x
    return 0.0 if z <= 0 else -math.sqrt(z)


def FuncSphereY2(x):
    z = sr*sr-x*x
    return 0.0 if z <= 0 else math.sqrt(z)


def FuncCircle(x):
    """!
    与えられた点が円の内なら1，外なら0を返す関数
    @param[in] x 変数(2要素のリスト)
    @return 円の内なら1，外なら0
    """
    return 1 if x[0]*x[0]+x[1]*x[1] <= sr*sr else 0


# -----------------------------------------------------------------------------
# 常微分方程式(ODE:Ordinary Differential Equation)用
# -----------------------------------------------------------------------------
# f(x,y)=-λy の λ (Pythonではlambdaが予約語なのでlmbdaとしている)
lmbda = 1.0


def FuncOdeY(x, y):
    """!
    f(x,y)=-λy
     - dy/dx=f(x,y)の真値 : y = C e^-λx
    @param[in] x,y 変数
    @return f(x,y)の値
    """
    return -lmbda*y


def DyFuncOdeY(x, y):
    return -lmbda


def FuncOdeY_true(x, C):
    """! 微分方程式の真値 """
    return C*math.exp(-lmbda*x)


def FuncOdeXY(x, y):
    """!
    f(x,y)=xy (変数分離型)
     - dy/dx=f(x,y)の真値 : y = C e^(x^2)
    @param[in] x,y 変数
    @return f(x,y)の値
    """
    return 2*x*y


def DyFuncOdeXY(x, y):
    return 2*x


def FuncOdeXY_true(x, C):
    """! 微分方程式の真値 """
    return C*math.exp(x*x)


# ロトカ・ヴォルテラ方程式の係数
lva, lvb = 0.01, 0.0001   # 被食者の増殖係数と被食者による減少係数
lvc, lvd = 0.0001, 0.05   # 捕食者の被食者の数による増殖係数と捕食者の数が増えることによる減少係数


def FuncOdeLV(t, y):
    """!
    ロトカ・ヴォルテラ方程式
     - 多変数の場合の例
     - 被食者と捕食者の生存競争をモデル化
     - 現段階で解析的な解は求められないことが分かっている
    @param[in] t,y 変数 (y[0]:x,y[1]:y)
    @return f(t,x,y)の値
    """
    f = [0.0]*len(y)
    f[0] = lva*y[0]-lvb*y[0]*y[1]
    f[1] = lvc*y[0]*y[1]-lvd*y[1]
    return f


# 単振り子の係数
g = 9.8    # 重力加速度
lp = 1.0   # 振り子のひもの長さ(回転中心からの距離)


def FuncOdePendulum(x, y):
    """!
    単振り子
     - 多変数の場合の例2
    @param[in] x,y 変数 (y[0]:θ,y[1]:ω)
    @return f(x,y)の値
    """
    f = [0.0]*len(y)
    f[0] = y[1]
    f[1] = -g*lp*math.sin(y[0])
    return f


# -----------------------------------------------------------------------------
# 偏微分方程式(PDE:Partial Differential Equation)用
# -----------------------------------------------------------------------------
def FuncPdeX3(x):
    """!
    g(x)=-20x^3
     - d^2 f/dx^2=g(x)の真値 : f(x)=x-x^5
    @param[in] x 変数
    @return g(x)の値
    """
    return -20*x*x*x


def FuncPdeX3T(x):
    """! 微分方程式の真値 """
    return x-x*x*x*x*x
