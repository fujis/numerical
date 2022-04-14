/*! 
  @file gauss.cpp
	
  @brief 数値積分法
         ガウス型積分公式
 
  @author Makoto Fujisawa
  @date 2019-08
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

//-----------------------------------------------------------------------------
// 数値積分法
//-----------------------------------------------------------------------------
/*!
 * ガウス・ルジャンドル公式による積分計算
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[out] S 解
 */
int gauss2(double func(const double), double a, double b, double &S)
{
	double x[2], w[2];

	// 分点と重みの計算(n=2)
	x[0] = -sqrt(3.0);
	x[1] = -x[0];
	w[0] = w[1] = 1;

	// Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
	S = 0.0;
	for(int i = 0; i < 2; ++i){
		S += w[i]*func((b-a)*x[i]/2+(a+b)/2);
	}
	S *= (b-a)/2;

	return 0;
}
int gauss3(double func(const double), double a, double b, double &S)
{

	// 分点と重みの計算(n=3)
	double x[3], w[3];
	x[0] = -sqrt(3.0/5.0);
	x[1] = 0;
	x[2] = -x[0];
	w[0] = w[2] = 5.0/9.0;
	w[1] = 8.0/9.0;

	// Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
	S = 0.0;
	for(int i = 0; i < 3; ++i){
		S += w[i]*func((b-a)*x[i]/2+(a+b)/2);
	}
	S *= (b-a)/2;

	return 0;
}
int gauss4(double func(const double), double a, double b, double &S)
{
	double x[4], w[4];

	// 分点と重みの計算(n=2)
	double tmp = 2.0*sqrt(6.0/5.0);
	x[0] = -sqrt((3.0+tmp)/7.0);
	x[1] = -sqrt((3.0-tmp)/7.0);;
	x[2] = -x[1];
	x[3] = -x[0];
	w[0] = w[3] = (18.0-sqrt(30))/36.0;
	w[1] = w[2] = (18.0+sqrt(30))/36.0;

	// Σ wi f(xi)の計算([a,b]と[-1,1]の変換付き)
	S = 0.0;
	for(int i = 0; i < 4; ++i){
		S += w[i]*func((b-a)*x[i]/2+(a+b)/2);
	}
	S *= (b-a)/2;

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 指数関数の[0,1]での積分 (解析解はe-1=1.718281828459045235360287471352...)
	//double(*func)(double) = FuncExp;
	//double a = 0.0, b = 1.0;
	//double t = (exp(1.0)-1); // 真値

	// 講義で示した例題(2017年度筑波大学前期日程入試問題)
	double(*func)(double) = FuncT17;
	double a = 0.5, b = 2.0;
	double t = -99.0/8.0+18.0*log(2.0); // 真値

	cout.precision(10);
	double s = 0.0;

	cout << endl;
	gauss2(func, a, b, s);
	cout << "gauss2 int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	gauss3(func, a, b, s);
	cout << "gauss3 int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	gauss4(func, a, b, s);
	cout << "gauss4 int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	cout << "ground truth = " << t << endl;

	return 0;
}


