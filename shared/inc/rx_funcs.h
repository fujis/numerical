/*!
  @file rx_funcs.h

  @brief 数値計算テスト：方程式の記述

  @author Makoto Fujisawa
  @date 2012-06
*/


#ifndef _RX_FUNC_H_
#define _RX_FUNC_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"

using namespace std;


//-----------------------------------------------------------------------------
// 関数定義
//-----------------------------------------------------------------------------

/*!
 * 1次元の代数方程式
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func1(const double x)
{
	//return x*x-2;
	return 2*x*x*x*x*x+5*x*x*x+3*x+1;
}
inline double DFunc1(const double x)
{
	//return 2*x;
	return 10*x*x*x*x+15*x*x+3;
}


/*!
 * 1次元の凸関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func2(const double x)
{
	return x*x-2*x+2;
}
inline double DFunc2(const double x)
{
	return 2*x-2;
}


/*!
 * 1次元の非凸関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func3(const double x)
{
	return x*cos(2*x);
}
inline double DFunc3(const double x)
{
	return cos(2*x)-2*x*sin(2*x);
}



/*!
 * 2次元の凸関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func4(const vector<double> &x)
{
	return x[0]*x[0]+x[1]*x[1]-x[0]*x[1]-x[0]-1;
	//return x[0]*x[0]-4*x[0]*x[1]+x[1]*x[1];
}
inline vector<double> DFunc4(const vector<double> &x)
{
	vector<double> g(x.size(), 0.0);
	g[0] = 2*x[0]-x[1]-1;
	g[1] = 2*x[1]-x[0];
	//g[0] = 2*x[0]-4*x[1];
	//g[1] = -4*x[0]+2*x[1];
	return g;
}

/*!
 * 2次元の凸関数2
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func4a(const vector<double> &x)
{
	return x[0]*x[0]+x[1]*x[1]-2;
}
inline vector<double> DFunc4a(const vector<double> &x)
{
	vector<double> g(x.size(), 0.0);
	g[0] = 2*x[0];
	g[1] = 2*x[1];
	return g;
}

/*!
 * 2次元関数(鞍点が生じる例)
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double Func5(const vector<double> &x)
{
	return 2*x[0]*x[0]-x[0]*x[1]-6*x[1]*x[1]+2*x[1];
}
inline vector<double> DFunc5(const vector<double> &x)
{
	vector<double> g(x.size(), 0.0);
	g[0] = 4*x[0]-x[1];
	g[1] = -x[0]-12*x[1]+2;
	return g;
}


/*!
 * 1次元の指数関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncExp(const double x)
{
	return exp(x);
}
inline double DFuncExp(const double x)
{
	return exp(x);
}

/*!
 * 円周率計算用関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncPi(const double x)
{
	return cos(x/2.0);
}
inline double DFuncPi(const double x)
{
	return -sin(x/2.0)/2.0;
}


/*!
 * 1次元の1次関数
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncLinear(const double x)
{
	return 2.0*x+1.0;
}
inline double DFuncLinear(const double x)
{
	return 2.0;
}

/*!
 * 1次元のルンゲ関数
 * @param[in] x 変数
 * @return 方程式の値
 */
const double ainv = 25;
inline double FuncRunge(const double x)
{
	return 1.0/(1.0+ainv*x*x);
}
inline double DFuncRunge(const double x)
{
	return -2*ainv*x/((1.0+ainv*x*x)*(1.0+ainv*x*x));
}

/*!
 * 2017年度筑波大前期日程入試問題[4]の数式
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncT17(const double x)
{
	return 2*x*x-9*x+14-9/x+2/(x*x);
}
inline double DFuncT17(const double x)
{
	return 4*x-9+9/(x*x)-4/(x*x*x);
}



/*!
 * 2次元関数(積分用2次多項式)
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncP2(const double x, const double y)
{
	return 8*x*x+4*y;
}
inline double FuncY1(const double x)
{
	return 2-x;
}
inline double FuncY2(const double x)
{
	return x*x;
}



//-----------------------------------------------------------------------------
// 円,球の面積,体積計算用
//-----------------------------------------------------------------------------
const double sr = 1.0; // 円,球の半径
/*!
 * 円の上半分の形状
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncCircleTop(const double x)
{
	double y = sr*sr-x*x;
	return (y <= 0 ? 0.0 : sqrt(y));
}
inline double DFuncCircleTop(const double x)
{
	double y = sr*sr-x*x;
	return (y <= 1e-10 ? 0.0 : -x/sqrt(y));
}
/*!
 * 円の下半分の形状
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncCircleBottom(const double x)
{
	double y = sr*sr-x*x;
	return (y <= 0 ? 0.0 : -sqrt(y));
}
inline double DFuncCircleBottom(const double x)
{
	double y = sr*sr-x*x;
	return (y <= 1e-10 ? 0.0 : x/sqrt(y));
}

/*!
 * 2次元関数(球体の体積計算用)
 * @param[in] x 変数
 * @return 方程式の値
 */
inline double FuncSphere(const double x, const double y)
{
	double z = sr*sr-x*x-y*y;
	return (z <= 0 ? 0.0 : sqrt(z));
}
inline double FuncSphereY1(const double x)
{
	double z = sr*sr-x*x;
	return (z <= 0 ? 0.0 : -sqrt(z));
}
inline double FuncSphereY2(const double x)
{
	double z = sr*sr-x*x;
	return (z <= 0 ? 0.0 : sqrt(z));
}

/*!
 * 与えられた点が円の内なら1，外なら0を返す関数
 * @param[in] x 変数
 * @return 円の内なら1，外なら0
 */
inline int FuncCircle(const vector<double> &x)
{
	return (x[0]*x[0]+x[1]*x[1] <= sr*sr ? 1 : 0);
}


//-----------------------------------------------------------------------------
// 常微分方程式(ODE:Ordinary Differential Equation)用
//-----------------------------------------------------------------------------

/*!
 * f(x,y)=-λy
 *  - dy/dx=f(x,y)の真値 : y = C e^-λx
 * @param[in] x,y 変数
 * @return f(x,y)の値
 */
static double lambda = 1.0;
inline double FuncOdeY(double x, double y)
{
	return -lambda*y;
}
inline double DyFuncOdeY(double x, double y)
{
	return -lambda;
}
inline double FuncOdeY_true(double x, double C) // 微分方程式の真値
{
	return C*exp(-lambda*x);
}

/*!
 * f(x,y)=xy (変数分離型)
 *  - dy/dx=f(x,y)の真値 : y = C e^(x^2)
 * @param[in] x,y 変数
 * @return f(x,y)の値
 */
inline double FuncOdeXY(double x, double y)
{
	return 2*x*y;
}
inline double DyFuncOdeXY(double x, double y)
{
	return 2*x;
}
inline double FuncOdeXY_true(double x, double C) // 微分方程式の真値
{
	return C*exp(x*x);
}

/*!
 * ロトカ・ヴォルテラ方程式
 *  - 多変数の場合の例
 *  - 被食者と捕食者の生存競争をモデル化
 *  - 現段階で解析的な解は求められないことが分かっている
 * @param[in] t,y 変数 (y[0]:x,y[1]:y)
 * @return f(t,x,y)の値
 */
static double lva = 0.01, lvb = 0.0005; // 被食者の増殖係数と被食者による減少係数
static double lvc = 0.0005, lvd = 0.05; // 捕食者の被食者の数による増殖係数と捕食者の数が増えることによる減少係数
inline vector<double> FuncOdeLV(double t, const vector<double> &y)
{
	vector<double> f(y.size(), 0.0);
	f[0] = lva*y[0]-lvb*y[0]*y[1];
	f[1] = lvc*y[0]*y[1]-lvd*y[1];
	return f;
}

/*!
 * 単振り子
 *  - 多変数の場合の例2
 * @param[in] t,y 変数 (y[0]:θ,y[1]:ω)
 * @return f(x,y)の値
 */
static double g = 9.8; // 重力加速度
static double lp = 1.0; // 振り子のひもの長さ(回転中心からの距離)
inline vector<double> FuncOdePendulum(double x, const vector<double> &y)
{
	vector<double> f(y.size(), 0.0);
	f[0] = y[1];
	f[1] = -g*lp*sin(y[0]);
	return f;
}


//-----------------------------------------------------------------------------
// 偏微分方程式(PDE:Partial Differential Equation)用
//-----------------------------------------------------------------------------
/*!
 * g(x)=20x^3
 *  - d^2 f/dx^2=g(x)の真値 : f(x)=x-x^5
 * @param[in] x 変数
 * @return g(x)の値
 */
inline double FuncPdeX3(double x)
{
	return -20*x*x*x;
}
inline double FuncPdeX3T(double x) // 微分方程式の真値
{
	return x-x*x*x*x*x;
}


#endif // #ifndef _RX_FUNC_H_