/*! 
  @file odesystem.cpp
	
  @brief 連立常微分方程式(system of ordinary differential equations)
 
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
 * 前進オイラー法(多変数版)
 *  - y(n+1)=y(n)+hf(x(n),y(n))
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
vector<double> eular(vector<double> func(double,const vector<double>&), vector<double> y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	vector<double> y(y0), f(y0.size(), 0.0);
	cout << x << ", " << y << endl; // グラフ描画用
	for(int k = 0; k <= n-1; ++k){
		f = func(x, y);
		for(int i = 0; i < y0.size(); ++i){
			y[i] = y[i]+h*f[i]; // yの更新
		}
		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << endl; // グラフ描画用
	}

	return y;
}

/*!
 * ホイン法(多変数版)
 *  - y(n+1)=y(n)+h/2(f(x(n),y(n))+f(x(n+1),y(n+1)))
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
vector<double> heun(vector<double> func(double, const vector<double>&), vector<double> y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	vector<double> y(y0), f(y0.size(), 0.0);
	vector<double> Yi(y0), f1(y0.size(), 0.0);
	cout << x << ", " << y << endl; // グラフ描画用
	for(int k = 0; k <= n-1; ++k){
		f = func(x, y);
		for(int i = 0; i < y0.size(); ++i){
			Yi[i] = y[i]+h*f[i]; // Yiの計算
		}

		// f(x,y)とY(i+1)の平均でyを更新
		f1 = func(x+h, Yi);
		for(int i = 0; i < y0.size(); ++i){
			y[i] = y[i]+h*(f[i]+f1[i])/2.0;
		}

		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << endl; // グラフ描画用
	}

	return y;
}



/*!
 * 4段4次のルンゲ・クッタ法(多変数版)
 *  - y(n+1)=y(n)+(h/6)(k1+2*k2+2*k3+k4)
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
vector<double> rk4(vector<double> func(double, const vector<double>&), vector<double> y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	vector<double> y(y0), f(y0.size(), 0.0);
	vector<double> k1(y0), k2(y0), k3(y0), k4(y0), yk(y0);
	cout << x << ", " << y << endl; // グラフ描画用
	for(int k = 0; k <= n-1; ++k){
		k1 = func(x, y);

		for(int i = 0; i < y0.size(); ++i){
			yk[i] = y[i]+(h/2)*k1[i];
		}
		k2 = func(x+h/2, yk);

		for(int i = 0; i < y0.size(); ++i){
			yk[i] = y[i]+(h/2)*k2[i];
		}
		k3 = func(x+h/2, yk);

		for(int i = 0; i < y0.size(); ++i){
			yk[i] = y[i]+h*k3[i];
		}
		k4 = func(x+h, yk);

		// yの更新
		for(int i = 0; i < y0.size(); ++i){
			y[i] = y[i]+(h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
		}

		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << endl; // グラフ描画用
	}

	return y;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// Lotka-Volterraモデル
	//vector<double>(*func)(double, const vector<double>&) = FuncOdeLV;
	//double a = 0.0, b = 1000.0; // 範囲[a,b] → 時刻(日数)
	//int n = 1000;
	//vector<double> y0(2, 0.0); // 初期値
	//y0[0] = 300; y0[1] = 300;

	// 単振り子モデル
	vector<double>(*func)(double, const vector<double>&) = FuncOdePendulum;
	double a = 0.0, b = 10.0; // 範囲[a,b] → 時刻(秒数)
	int n = 100;
	vector<double> y0(2, 0.0); // 初期値
	y0[0] = RX_PI/4.0; y0[1] = 0;

	cout.precision(10);
	vector<double> y(y0);

	// 多変数版オイラー法(1次精度)
	//cout << "[Eular method]" << endl;
	//y = eular(func, y0, a, b, n);
	//cout << endl;

	// 多変数版ホイン法(2次精度)
	//cout << "[Heun method]" << endl;
	//y = heun(func, y0, a, b, n);
	//cout << endl;

	// 多変数版RK4(4次精度)
	cout << "[RK4]" << endl;
	y = rk4(func, y0, a, b, n);
	cout << endl;


	return 0;
}


