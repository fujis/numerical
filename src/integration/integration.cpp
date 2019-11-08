/*! 
  @file integration.cpp
	
  @brief 数値積分法
		 区分求積法と台形公式
 
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
 * 区分求積法
 *  - 積分区間を矩形で分割する方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数
 * @param[out] ans 解
 * @return 
 */
int segment_integration(double func(const double), double a, double b, int n, double &S)
{
	double h = (b-a)/n; // 分割区間の横幅

	S = 0;
	for(int i = 0; i < n; ++i){
		double f = func(a+(i+1)*h); // xの大きい方の辺の長さを縦幅とする
		S += f*h;
	}

	return 0;
}

/*!
 * 台形法(台形公式)
 *  - 積分区間を台形で分割する方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数
 * @param[out] ans 解
 * @return
 */
int trapezoidal_integration(double func(const double), double a, double b, int n, double &S)
{
	double h = (b-a)/n; // 分割区間の横幅
	double f1, f2;		// 分割区間の縦幅(台形の長辺と短辺)

	f2 = func(a);
	S = 0;
	for(int i = 0; i < n; ++i){
		f1 = f2;
		f2 = func(a+(i+1)*h);
		S += (f1+f2)*h/2;
	}

	return 0;
}


int monte_carlo(double func(const double), double a, double b, double fmin, double fmax, int n, double &S)
{
	int nu = 0;
	for(int i = 0; i < n; ++i){
		double rnd1 = (double)rand()/(double)RAND_MAX; // [0,1]の乱数1
		double rnd2 = (double)rand()/(double)RAND_MAX; // [0,1]の乱数2
		double x = a+(b-a)*rnd1;
		double y = fmin+(fmax-fmin)*rnd2;

		double f = func(x);
		if(y <= f) nu++;
	}
	S = (double)nu/(double)n*(b-a)*(fmax-fmin);

	return 0;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 指数関数の[0,1]での積分 (解析解はe-1=1.718281828459045235360287471352...)
	double(*func)(double) = FuncExp;
	double a = 0.0, b = 1.0;
	double t = (exp(1.0)-1); // 真値

	// 講義で示した例題(2017年度筑波大学前期日程入試問題)
	//double(*func)(double) = FuncT17;
	//double a = 0.5, b = 2.0;
	//double t = -99.0/8.0+18.0*log(2.0); // 真値

	cout.precision(10);
	int n = 8;
	double s = 0.0;

	segment_integration(func, a, b, n, s);
	cout << "segment int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	trapezoidal_integration(func, a, b, n, s);
	cout << "trapezoidal int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	cout << "ground truth = " << t << endl;

	monte_carlo(func, a, b, 0.0, 3.0, 10000, s);
	cout << "monte carlo int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	return 0;
}


