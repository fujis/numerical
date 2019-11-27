/*! 
  @file eular.cpp
	
  @brief オイラー法，ホイン法(改良オイラー法)
 
  @author Makoto Fujisawa
  @date 2019-10
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
//! 常微分方程式の近似解
//-----------------------------------------------------------------------------
/*!
 * オイラー法(1次精度)
 *  - y(n+1)=y(n)+hf(x(n),y(n))
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double eular(double func(double,double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << y0*exp(-25*x) << ", " << fabs(y-y0*exp(-25*x)) << endl; // グラフ描画用に真値も出力
	for(int i = 0; i <= n-1; ++i){
		double fi = func(x, y);
		y = y+h*fi; // yの更新
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << y0*exp(-25*x) << ", " << fabs(y-y0*exp(-25*x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}

/*!
 * ホイン法(2次精度)
 *  - y(n+1)=y(n)+h/2(f(x(n),y(n))+f(x(n+1),y(n+1)))
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double heun(double func(double, double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << y0*exp(-25*x) << ", " << fabs(y-y0*exp(-25*x)) << endl; // グラフ描画用に真値も出力
	for(int i = 0; i <= n-1; ++i){
		// Y(i+1) = f(x(i+1)+y(i+1))をオイラー法で求める
		double fi = func(x, y);
		double Yi = y+h*fi; // Yiの計算

		// f(x,y)とY(i+1)の平均でyを更新
		y = y+h*(fi+Yi)/2;
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << y0*exp(-25*x) << ", " << fabs(y-y0*exp(-25*x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}

//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 常微分方程式dy/dx=ay (解析解はy = C e^ax = y0 e^ax)
	double(*func)(double,double) = FuncOdeY;
	double a = 0.0, b = 1.0; // 範囲[a,b]
	double y0 = 1.0; // 初期値
	double t = y0*exp(-25*b); // 真値

	// 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
	//double(*func)(double,double) = FuncOdeXY;
	//double a = 0.0, b = 1.0; // 範囲[a,b]
	//double y0 = 1.0; // 初期値
	//double t = y0*exp(b*b); // 真値

	cout.precision(10);
	int n = 8;
	double y = 0.0;

	// オイラー法(1次精度)
	cout << "[Eular method]" << endl;
	y = eular(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	// ホイン法(2次精度)
	cout << "[Heun method]" << endl;
	y = heun(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	cout << "ground truth = " << t << endl;





	return 0;
}


