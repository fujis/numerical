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
	return x*x*x-2*x*x-x+2;
}
inline double DFunc1(const double x)
{
	return 3*x*x-4*x-1;
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
	//return x[0]*x[0]-2*x[0]*x[1]+2*x[1]*x[1];
	return x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1]+3*x[0];
}
inline vector<double> DFunc4(const vector<double> &x)
{
	vector<double> g(x.size(), 0.0);
	//g[0] = 2*x[0]-2*x[1];
	//g[1] = -2*x[0]+4*x[1];
	g[0] = 2*x[0]+2*x[1]+3;
	g[1] = 2*x[0]+6*x[1];
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


#endif // #ifndef _RX_FUNC_H_