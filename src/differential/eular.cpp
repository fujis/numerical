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
// デバッグ用変数
//-----------------------------------------------------------------------------
std::function<double(double)> TF;

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
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	for(int i = 0; i <= n-1; ++i){
		double fi = func(x, y);
		y = y+h*fi; // yの更新
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
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
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	for(int i = 0; i <= n-1; ++i){
		// Y(i+1) = f(x(i+1)+y(i+1))をオイラー法で求める
		double fi = func(x, y);
		double Yi = y+h*fi; // Yiの計算

		// f(x,y)とY(i+1)の平均でyを更新
		y = y+h*(fi+func(x+h, Yi))/2.0;
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}

/*!
 * 後退オイラー法(1次精度)
 *  - y(n+1)=y(n)+hf(x(n+1),y(n+1))
 *  - ニュートン法を使って計算
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double backward_eular(double func(double, double), double dfunc(double, double), double y0, double a, double b, int n, int max_iter, double eps)
{
	int k_avg = 0;
	double h = (b-a)/n; // 刻み幅
	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	for(int i = 0; i <= n-1; ++i){
		double yk = y; // y_(i+1)計算用の一時変数

		// ニュートン法でy_(i+1)を求める
		int k;
		for(k = 0; k < max_iter; ++k){
			double g = yk-h*func(x+h, yk)-y;
			double dg = 1-h*dfunc(x+h, yk);
			yk = yk-g/dg;
			if(fabs(g/dg) < eps || fabs(g) < eps) break;
		}
		k_avg += k;

		y = yk; // yの更新
		//y = y/(1+25*h);
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}
	cout << "average number of iterations : " << k_avg/(double)n << endl;

	return y;
}




//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 常微分方程式dy/dx=ay (解析解はy = C e^ax = y0 e^ax)
	double(*func)(double, double) = FuncOdeY;
	double(*dfunc)(double, double) = DyFuncOdeY;
	double a = 0.0, b = 1.0; // 範囲[a,b]
	double y0 = 1.0; // 初期値
	lambda = 25;
	TF = std::bind(FuncOdeY_true, placeholders::_1, y0);// 真値

	//// 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
	//double(*func)(double, double) = FuncOdeXY;
	//double(*dfunc)(double, double) = DyFuncOdeXY;
	//double a = 0.0, b = 1.0; // 範囲[a,b]
	//double y0 = 1.0; // 初期値
	//TF = std::bind(FuncOdeXY_true, placeholders::_1, y0);// 真値

	cout.precision(10);
	int n = 10;
	double y = 0.0;
	double t = TF(b);  // 真値

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

	// 後退オイラー法(1次精度)
	cout << "[backward Eular method]" << endl;
	y = backward_eular(func, dfunc, y0, a, b, n, 30, 1e-6);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	cout << "ground truth = " << t << endl;





	return 0;
}


