/*! 
  @file integration.chpp
	
  @brief 数値積分法
 
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
int segment_integration(double func(const double), double a, double b, int n, double &ans)
{
	double h = (b-a)/n; // 分割区間の横幅

	ans = 0;
	for(int i = 0; i < n; ++i){
		double f = func(a+(i+1)*h); // xの大きい方の辺の長さを縦幅とする
		ans += f*h;
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
int trapezoidal_integration(double func(const double), double a, double b, int n, double &ans)
{
	double h = (b-a)/n; // 分割区間の横幅
	double f1, f2;		// 分割区間の縦幅(台形の長辺と短辺)

	f2 = func(a);
	ans = 0;
	for(int i = 0; i < n; ++i){
		f1 = f2;
		f2 = func(a+(i+1)*h);
		ans += (f1+f2)*h/2;
	}

	return 0;
}


/*!
 * シンプソン法(シンプソン公式)
 *  - 台形法の誤差を評価してより精度を高めた方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数
 * @param[out] ans 解
 * @return
 */
int simpson_integration(double func(const double), double a, double b, int n, double &ans)
{
	double h = (b-a)/n; // 分割区間の横幅
	double f;			// 分割区間の縦幅

	ans = func(a)+func(b);
	// 奇数項
	for(int i = 0; i < n/2; ++i){
		f = func(a+(2*i+1)*h);
		ans += 4*f;
	}
	// 偶数項
	for(int i = 1; i < n/2; ++i){
		f = func(a+(2*i)*h);
		ans += 2*f;
	}
	ans *= h/3;

	return 0;
}




//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	//double(*func)(double) = FuncExp;
	//double a = 0.0, b = 1.0;
	//double t = (exp(1.0)-1); // 真値

	double(*func)(double) = FuncT17;
	double a = 0.5, b = 2.0;
	double t = -99.0/8.0+18.0*log(2.0); // 真値

	cout.precision(10);
	int n = 100;
	double s = 0.0;

	// 指数関数の[0,1]での積分 (解析解はe-1=1.718281828459045235360287471352...)
	segment_integration(func, a, b, n, s);
	cout << "segment int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	trapezoidal_integration(func, a, b, n, s);
	cout << "trapezoidal int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	simpson_integration(func, a, b, n, s);
	cout << "simpson int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	cout << "ground truth = " << t << endl;

		
	return 0;
}


