/*! 
  @file adams.cpp
	
  @brief アダムス・バッシュホース法
 
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
 * アダムス・バッシュホース法(3点)
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double adamsbashforth3(double func(double,double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double fi, fi1 = 0, fi2 = 0;  // f_i, f_(i-1), f_(i-2)
	for(int i = 0; i <= n-1; ++i){
		fi = func(x, y);
		if(i < 2){ // 最初の方はオイラー法を使う
			y = y+h*fi; // yの更新
		}
		else{
			y = y+h*(5*fi2-16*fi1+23*fi)/12.0; // yの更新
		}

		fi2 = fi1;
		fi1 = fi;
		x = x+h;    // xの更新

		// 出力用
		cout << "y(" << x << ") = " << y << endl;
		//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}

/*!
 * アダムス・バッシュホース法(4点)
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double adamsbashforth4(double func(double, double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	//cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double fi, fi1 = 0, fi2 = 0, fi3 = 0;  // f_i, f_(i-1), f_(i-2), f_(i-3)
	for(int i = 0; i <= n-1; ++i){
		fi = func(x, y);
		if(i < 3){ // 最初の方はオイラー法を使う
			y = y+h*fi; // yの更新
		} else{
			y = y+h*(-9*fi3+37*fi2-59*fi1+55*fi)/24.0; // yの更新
		}

		fi3 = fi2;
		fi2 = fi1;
		fi1 = fi;
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
	//double(*func)(double,double) = FuncOdeY;
	//double a = 0.0, b = 1.0; // 範囲[a,b]
	//double y0 = 1.0; // 初期値
	//lambda = 25;
	//TF = std::bind(FuncOdeY_true, placeholders::_1, y0);// 真値

	// 常微分方程式dy/dx=2xy (解析解はy = C e^(x^2) = y0 e^(x^2))
	double(*func)(double, double) = FuncOdeXY;
	double a = 0.0, b = 1.0; // 範囲[a,b]
	double y0 = 1.0; // 初期値
	TF = std::bind(FuncOdeXY_true, placeholders::_1, y0);// 真値

	cout.precision(10);
	int n = 10;
	double y = 0.0;
	double t = TF(b);  // 真値

	// アダムス・バッシュホース法(3点,3次精度)
	cout << "[Adams-Bashforth(3)]" << endl;
	y = adamsbashforth3(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	// アダムス・バッシュホース法(4点,4次精度)
	cout << "[Adams-Bashforth(4)]" << endl;
	y = adamsbashforth4(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	cout << "ground truth = " << t << endl;





	return 0;
}


