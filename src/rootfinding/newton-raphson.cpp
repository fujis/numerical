/*! 
  @file newton-raphson.cpp
	
  @brief ニュートン・ラフソン法(Newton-Raphson method)
         もしくは単にニュートン法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// ニュートン法による求根問題の解法
//-----------------------------------------------------------------------------
/*!
 * ニュートン・ラフソン法(Newton-Raphson method)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] dfunc 導関数値を与える関数ポインタ
 * @param[inout] x 解(探索開始位置を渡す)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 */
int newton(double func(const double), double dfunc(const double), 
		   double &x, int &max_iter, double &eps)
{
	double f, df, dx;

	int k;
	for(k = 0; k < max_iter; ++k){
		// 現在の位置xにおける関数値と導関数の計算
		f = func(x);
		df = dfunc(x);

		// 確認用の画面出力
		cout << k << " : f(" << x << ") = " << f << endl;

		// 導関数の結果から次の位置を計算
		x = x-f/df;

		// 収束判定
		dx = fabs(f/df);
		if(dx < eps || fabs(f) < eps){
			break;
		}
	}

	max_iter = k; eps = dx;
	
	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 探索開始位置
	double x = -1;

	// ニュートン法でf(x)=0を解く
	int max_iter = 100;
	double eps = 1e-6;
	newton(Func1, DFunc1, x, max_iter, eps);

	// 結果の画面表示
	cout << "x = " << x << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


