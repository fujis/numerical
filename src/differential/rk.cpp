/*! 
  @file rk.cpp
	
  @brief ルンゲ・クッタ法
 
  @author Makoto Fujisawa
  @date 2019-11
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
 * 3段3次のルンゲ・クッタ法(クッタの3次公式)(3次精度)
 *  - y(n+1)=y(n)+(h/6)(k1+4*k2+k3)
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double rk3(double func(double,double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double k1, k2, k3;
	for(int i = 0; i <= n-1; ++i){
		k1 = func(x, y);				// k1の算出
		k2 = func(x+h/2, y+(h/2)*k1);	// k2の算出
		k3 = func(x+h, y-h*k1+2*h*k2);	// k3の算出

		y = y+(h/6)*(k1+4*k2+k3); // yの更新
		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}


/*!
 * 4段4次のルンゲ・クッタ法(4次精度)
 *  - y(n+1)=y(n)+(h/6)(k1+2*k2+2*k3+k4)
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double rk4(double func(double, double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double k1, k2, k3, k4;
	for(int i = 0; i <= n-1; ++i){
		k1 = func(x, y);				// k1の算出
		k2 = func(x+h/2, y+(h/2)*k1);	// k2の算出
		k3 = func(x+h/2, y+(h/2)*k2);	// k3の算出
		k4 = func(x+h, y+h*k3);			// k4の算出

		y = y+(h/6)*(k1+2*k2+2*k3+k4); // yの更新
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
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
	lambda = 25;
	TF = std::bind(FuncOdeY_true, placeholders::_1, y0);// 真値

	// 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
	//double(*func)(double, double) = FuncOdeXY;
	//double a = 0.0, b = 1.0; // 範囲[a,b]
	//double y0 = 1.0; // 初期値
	//TF = std::bind(FuncOdeXY_true, placeholders::_1, y0);// 真値

	cout.precision(10);
	int n = 20;
	double y = 0.0;
	double t = TF(b);  // 真値

	// 3段3次のルンゲ・クッタ法(3次精度)
	cout << "[RK3]" << endl;
	y = rk3(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	// 4段4次のルンゲ・クッタ法(3次精度)
	cout << "[RK4]" << endl;
	y = rk4(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	cout << "ground truth = " << t << endl;





	return 0;
}


