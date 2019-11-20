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
 * @return 積分値
 */
double segment_integration(double func(const double), double a, double b, int n)
{
	double h = (b-a)/n; // 分割区間の横幅

	double S = 0;
	for(int i = 0; i < n; ++i){
		double f = func(a+(i+1)*h); // xの大きい方の辺の長さを縦幅とする
		S += f*h;
	}

	return S;
}

/*!
 * 台形法(台形公式)
 *  - 積分区間を台形で分割する方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数
 * @return 積分値
 */
double trapezoidal_integration(double func(const double), double a, double b, int n)
{
	double h = (b-a)/n; // 分割区間の横幅
	double f1, f2;		// 分割区間の縦幅(台形の長辺と短辺)

	f2 = func(a);
	double S = 0;
	for(int i = 0; i < n; ++i){
		f1 = f2;
		f2 = func(a+(i+1)*h);
		S += (f1+f2)*h/2;
	}

	return S;
}


/*!
 * 2重積分(台形公式)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b x方向積分範囲
 * @param[in] y1f,y2f y方向積分範囲を与える関数ポインタ
 * @param[in] n,m x,y方向の積分区間の分割数
 * @return 積分値
 */
double trapezoidal_integration2(double func(const double, const double), double a, double b, double y1f(const double), double y2f(const double), int n, int m)
{
	double h1 = (b-a)/n; // x方向刻み幅
	double S = 0.0;
	for(int i = 0; i <= n; ++i){
		double xi = a+i*h1;
		double h2 = (y2f(xi)-y1f(xi))/m; // y方向刻み幅
		double y1 = y1f(xi);
		if(fabs(h2) < 1e-10) continue;

		// 台形公式によるy方向積分(Fiの計算)
		double f1, f2;
		f2 = func(xi, y1);
		double Fi = 0.0;
		for(int j = 0; j < n; ++j){
			f1 = f2;
			f2 = func(xi, y1+(j+1)*h2);
			Fi += (f1+f2)*h2/2;
		}

		if(i == 0 || i == n) S += Fi;
		else S += 2*Fi;
	}
	S *= h1/2;

	return S;
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
	int n = 100;
	double s = 0.0;

	s = segment_integration(func, a, b, n);
	cout << "segment int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	s = trapezoidal_integration(func, a, b, n);
	cout << "trapezoidal int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	cout << "ground truth = " << t << endl;


	// 円の面積の計算(上半分の積分-下半分の積分)
	double r = 1.0;
	a = -r; b = r;
	t = RX_PI*r*r;
	s = trapezoidal_integration(FuncCircleTop, a, b, n)-trapezoidal_integration(FuncCircleBottom, a, b, n);
	cout << "area of circle = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	cout << "ground truth = " << t << endl;


	// 重積分
	//a = 1, b = 2;
	//t = 54; // 真値
	//s = trapezoidal_integration2(FuncP2, a, b, FuncY1, FuncY2, n, n);
	//cout << "trapezoidal int2 f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	//cout << "ground truth = " << t << endl;

	// 球の体積(こちらは上半分の積分x2で計算)
	a = -sr, b = sr;
	t = (4*RX_PI*sr*sr*sr)/3.0; // 真値
	s = 2*trapezoidal_integration2(FuncSphere, a, b, FuncSphereY1, FuncSphereY2, n, n);
	cout << "trapezoidal int2 f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	cout << "ground truth = " << t << endl;


	return 0;
}


