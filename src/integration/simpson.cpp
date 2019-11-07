/*! 
  @file integration.cpp
	
  @brief 数値積分法
		 シンプソン公式
 
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
 * シンプソン法(シンプソンの1/3公式)
 *  - 台形法の誤差を評価してより精度を高めた方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数(2区間で1つの多項式なのでここでの分割数はデータ点数/2)
 * @param[out] ans 解
 * @return
 */
int simpson_integration(double func(const double), double a, double b, int n, double &S)
{
	double dx = (b-a)/(2*n); // 分割区間の横幅
	double f;			// 分割区間の縦幅

	S = func(a)+func(b);

	// 奇数項
	for(int i = 1; i <= n; ++i){
		f = func(a+(2*i-1)*dx);
		S += 4*f;
	}
	// 偶数項
	for(int i = 1; i <= n-1; ++i){
		f = func(a+(2*i)*dx);
		S += 2*f;
	}
	S *= dx/3;

	return 0;
}

/*!
 * シンプソン法(シンプソンの3/8公式)
 *  - 台形法の誤差を評価してより精度を高めた方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数(3区間で1つの多項式なのでここでの分割数はデータ点数/3)
 * @param[out] ans 解
 * @return
 */
int simpson38_integration(double func(const double), double a, double b, int n, double &S)
{
	double dx = (b-a)/(3*n); // 分割区間の横幅
	double f;			// 分割区間の縦幅

	S = func(a)+func(b);

	// 3の倍数-2の項
	for(int i = 1; i <= n; ++i){
		f = func(a+(3*i-2)*dx);
		S += 3*f;
	}
	// 3の倍数-1の項
	for(int i = 1; i <= n; ++i){
		f = func(a+(3*i-1)*dx);
		S += 3*f;
	}
	// 3の倍数の項
	for(int i = 1; i <= n-1; ++i){
		f = func(a+(3*i)*dx);
		S += 2*f;
	}
	S *= 3*dx/8;

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
	int n = 120; // シンプソン公式のためにnは2と3の公倍数にしておくこと
	double s = 0.0;

	simpson_integration(func, a, b, n/2, s);
	cout << "simpson int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	simpson38_integration(func, a, b, n/3, s);
	cout << "simpson(3/8) int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	cout << "ground truth = " << t << endl;

		
	return 0;
}


