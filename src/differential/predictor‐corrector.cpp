/*! 
  @file predictor‐corrector.cpp
	
  @brief アダムス・バッシュホース法とアダムス・ムルトン法を追加した予測子修正子法
 
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
 * アダムス・バッシュホース法(3点)+アダムスムルトン法(3点)による予測子修正子法
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double predictorcorrector_abm3(double func(double,double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double f[4] = { 0, 0, 0, 0 };  // f_(i-2), f_(i-1), f_(i), f_(i+1)
	for(int i = 0; i <= n-1; ++i){
		f[2] = func(x, y);
		if(i < 2){ // 最初の方はオイラー法を使う
			// 前進オイラー法でy_(i+1)の予測値を計算
			f[3] = func(x+h, y+h*f[2]);

			// 改良オイラー法で解を修正
			y = y+h*(f[2]+f[3])/2.0;
		}
		else{
			// アダムス・バッシュフォース法でy_(i+1)の予測値を計算
			double y1 = y+h*(5*f[0]-16*f[1]+23*f[2])/12.0;
			f[3] = func(x+h, y1);

			// アダムス・ムルトン法で解を修正
			y = y+h*(-f[1]+8*f[2]+5*f[3])/12.0;
		}

		f[0] = f[1];
		f[1] = f[2];
		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}

/*!
 * アダムス・バッシュホース法(4点)+アダムスムルトン法(4点)による予測子修正子法
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double predictorcorrector_abm4(double func(double, double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double f[5] = { 0, 0, 0, 0, 0 };  // f_(i-3), f_(i-2), f_(i-1), f_(i), f_(i+1)
	for(int i = 0; i <= n-1; ++i){
		f[3] = func(x, y);
		if(i < 3){ // 最初の方はオイラー法を使う
			// 前進オイラー法でy_(i+1)の予測値を計算
			f[4] = func(x+h, y+h*f[3]);

			// 改良オイラー法で解を修正
			y = y+h*(f[3]+f[4])/2.0;
		}
		else{
			// アダムス・バッシュフォース法でy_(i+1)の予測値を計算
			double y1 = y+h*(-9*f[0]+37*f[1]-59*f[2]+55*f[3])/24.0;
			f[4] = func(x+h, y1);

			// アダムス・ムルトン法で解を修正
			y = y+h*(f[1]-5*f[2]+19*f[3]+9*f[4])/24.0;
		}

		f[0] = f[1];
		f[1] = f[2];
		f[2] = f[3];
		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	}

	return y;
}


/*!
 * 前進オイラー+改良オイラーによる予測子修正子法
 * @param[in] func 関数f(x,y)の値を与える関数ポインタ
 * @param[in] y0 初期値(y0=g(a))
 * @param[in] a,b 計算範囲
 * @param[in] n 計算半以内での分割数(h=(b-a)/n)
 * @return x=bでの解
 */
double predictorcorrector_fbe(double func(double, double), double y0, double a, double b, int n)
{
	double h = (b-a)/n; // 刻み幅

	double x = a;  // xの初期値
	double y = y0; // yの初期値
	cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
	double f[2] = { 0, 0 };  // f_(i), f_(i+1)
	for(int i = 0; i <= n-1; ++i){
		f[0] = func(x, y);
		
		// 前進オイラー法でy_(i+1)の予測値を計算
		double y1 = y+h*f[0];
		f[1] = func(x+h, y1);

		// 改良オイラー法で解を修正
		y = y+h*(f[0]+f[1])/2.0;

		x = x+h;    // xの更新

		// 出力用
		//cout << "y(" << x << ") = " << y << endl;
		cout << x << ", " << y << ", " << TF(x) << ", " << fabs(y-TF(x)) << endl; // グラフ描画用に真値も出力
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
	//double(*func)(double,double) = FuncOdeXY;
	//double a = 0.0, b = 1.0; // 範囲[a,b]
	//double y0 = 1.0; // 初期値
	//TF = std::bind(FuncOdeXY_true, placeholders::_1, y0);// 真値

	cout.precision(10);
	int n = 10;
	double y = 0.0;
	double t = TF(b);  // 真値

	// 予測子修正子法(4次精度)
	cout << "[Predictor-Corrector method with Adams-Bashforth(4) + Adams-Moulton(4)]" << endl;
	y = predictorcorrector_abm4(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	// 予測子修正子法(3次精度)
	cout << "[Predictor-Corrector method with Adams-Bashforth(3) + Adams-Moulton(3)]" << endl;
	y = predictorcorrector_abm3(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	// 予測子修正子法(2次精度)
	cout << "[Predictor-Corrector method with Forward Eular + Modified Eular]" << endl;
	y = predictorcorrector_fbe(func, y0, a, b, n);
	cout << "y(" << b << ") = " << y << ",  error = " << fabs(y-t) << endl;
	cout << endl;

	cout << "ground truth = " << t << endl;





	return 0;
}


