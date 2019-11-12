/*! 
  @file monte-carlo.cpp
	
  @brief 数値積分法
		 モンテカルロ積分
 
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
 * モンテカルロ法(1次元)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] a,b 積分範囲
 * @param[in] n サンプル点個数
 * @return 積分値
 */
double monte_carlo(double func(const double), double a, double b, int n)
{
	double s = 0;
	for(int i = 0; i < n; ++i){
		double rnd = (double)rand()/(double)RAND_MAX; // [0,1]の乱数
		double x = a+(b-a)*rnd;
		double f = func(x);
		s += f;
	}
	return (b-a)*s/n;
}

/*!
 * Xorshift(32bit)による乱数生成
 *  - https://ja.wikipedia.org/wiki/Xorshift のコード例そのまま
 *  - オリジナルは G. Marsaglia, "Xorshift RNGs". Journal of Statistical Software, 8(14), 2013.
 */
uint32_t rand_xorshift(void)
{
	//static uint32_t y = 2463534242;
	static uint32_t y = time(0);
	y = y ^ (y << 13); y = y ^ (y >> 17);
	return y = y ^ (y << 5);
}

/*!
 * モンテカルロ法(面積)
 *  - 関数f(x)=1として2重積分で面積を計算する方法
 * @param[in] iofunc サンプル点が範囲内かどうかを判定する関数(範囲内で1,外で0を返す)
 * @param[in] x1,x2,y1,y2 積分範囲
 * @param[in] n 総サンプリング数
 * @return 積分値
 */
double monte_carlo_2d(int iofunc(const vector<double>&), double x1, double x2, double y1, double y2, int n)
{
	vector<double> x(2);
	int n_in = 0;
	for(int i = 0; i < n; ++i){
		double rnd1 = (double)rand()/(double)RAND_MAX; // [0,1]の乱数2
		double rnd2 = (double)rand()/(double)RAND_MAX; // [0,1]の乱数2
		//double rnd1 = (double)(rand_xorshift()%0xFFFFFF)/0xFFFFFF; // [0,1]の乱数1
		//double rnd2 = (double)(rand_xorshift()%0xFFFFFF)/0xFFFFFF; // [0,1]の乱数1
		x[0] = x1+(x2-x1)*rnd1;
		x[1] = y1+(y2-y1)*rnd2;
		if(iofunc(x)) n_in++;  // 範囲内のサンプル点個数を数える
	}
	double V = (x2-x1)*(y2-y1);
	return V*((double)n_in/(double)n);
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

	// モンテカルロ積分
	s = monte_carlo(func, a, b, 100);
	cout << "n=" << setw(4) << n << " : ";
	cout << "monte carlo int f(x) = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
	cout << "ground truth = " << t << endl;
	cout << endl;

	// 円の面積の計算(nを変えて実行)
	double r = sr;
	t = RX_PI*r*r;
	srand((unsigned int)time(0));
	n = 10;
	for(int i = 0; i < 7; ++i){
		s = monte_carlo_2d(FuncCircle, -r, r, -r, r, n);
		cout << "n=" << setw(10) << n << " : ";
		cout << "area of circle = " << fabs(s) << ",  error = " << fabs(fabs(s)-t) << endl;
		n *= 10;
	}
	cout << "ground truth = " << t << endl;

	// 誤差の平均値を求めるためのコード
	//int m = 6;
	//vector<double> avg_error(7, 0.0);
	//for(int j = 0; j < 10; ++j){
	//	srand((unsigned int)time(0));
	//	n = 10;
	//	for(int i = 0; i < m; ++i){
	//		s = monte_carlo_2d(FuncCircle, -r, r, -r, r, n);
	//		avg_error[i] += fabs(fabs(s)-t);
	//		n *= 10;
	//	}
	//}
	//n = 10;
	//for(int i = 0; i < m; ++i){
	//	avg_error[i] /= m;
	//	cout << n << ", " << avg_error[i] << endl;
	//	//cout << "n=" << setw(10) << n << " : " << "avg. error = " << avg_error[i] << endl;
	//	n *= 10;
	//}


	return 0;
}


