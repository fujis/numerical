/*! 
  @file romberg.cpp
	
  @brief 数値積分法
		 ロンバーグ法
 
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
 * 台形法(台形公式)
 *  - 積分区間を台形で分割する方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数
 * @param[out] ans 解
 * @return
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
 * ロンバーグ法
 *  - 台形法の誤差を評価してより精度を高めた方法
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n 積分区間の分割数(2区間で1つの多項式なのでここでの分割数はデータ点数/2)
 * @param[out] ans 解
 * @return
 */
int romberg(double func(const double), double a, double b, int &k_max, double &eps, double &S)
{
	double h = b-a; // 分割幅(初期分割幅はn=1のもの)
	double e = 0.0;		// 誤差
	vector<double> I(k_max+1, 0.0);

	I[0] = (h/2)*(func(a)+func(b)); // I_0,0の計算

	// 計算の途中経過確認用出力
	cout << setw(10) << I[0] << endl;

	int k, n = 1, l;
	for(k = 1; k <= k_max; ++k){
		h = h/2; n *= 2; // 分割幅を1/2にしていく
		I[k] = trapezoidal_integration(func, a, b, n); // 台形公式でI_k,0を計算

		// 計算の途中経過確認用出力
		cout << setw(10) << I[k];

		// 収束判定
		if((e = fabs(I[k]-I[k-1])) < eps){
			l = k; // 結果が格納されている位置
			break;
		}

		// 漸化式の計算
		int m, m4 = 4; // 4^m
		for(m = 1; m <= k; ++m){
			int i = k-m; // I_k,mを格納する配列上の位置
			I[i] = (m4*I[i+1]-I[i])/(m4-1); // I_k,mを計算
			m4 *= 4;

			// 計算の途中経過確認用出力
			cout << "   " << setw(10) << I[i];

			// 収束判定2
			if(i >= 1 && (e = fabs(I[i]-I[i-1])) < eps){
				l = i; // 結果が格納されている位置
				break;
			}
		}
		cout << endl;
		if(m <= k) break;
	}
	S = I[l];
	k_max = k;
	eps = e;

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
	double s = 0.0;

	int kmax = 5; // n=2^mmaxまで分割
	double eps = 1.0e-6;
	romberg(func, a, b, kmax, eps, s);
	cout << "romberg int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	cout << "n = " << (int)pow(2.0, (double)kmax) << ", eps = " << eps << endl;

	cout << "ground truth = " << t << endl;


	cout << endl;
	gauss2(func, a, b, s);
	cout << "gauss2 int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	gauss3(func, a, b, s);
	cout << "gauss3 int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;

	return 0;
}


