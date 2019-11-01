/*! 
  @file bisection.cpp
	
  @brief 二分法による求根問題の解法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// 二分法による求根問題の解法
//-----------------------------------------------------------------------------
/*!
 * 2分法(bisection method)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] x1,x2 探索範囲
 * @param[out] x 解
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 */
int bisection(double func(const double), double xl, double xr, double &x, int &max_iter, double &eps)
{
	double f = func(xl);
	double fmid = func(xr);

	// 探索範囲の境界での値の符号が異なる場合のみ探索
	if(f*fmid >= 0.0) return 0.0;

	int k;
	double dx = fabs(xr-xl), xmid;
	for(k = 0; k < max_iter; ++k){
		xmid = 0.5*(xl+xr);	// 中点
		dx *= 0.5;

		// 中点での関数値を求める
		fmid = func(xmid);

		// 確認用の画面出力
		cout << k << " : [" << xl << ", " << xr << "], fmid = " << fmid << endl;

		// 収束判定
		if(dx < eps || fmid == 0.0){
			break;
		}

		// 新しい区間
		if(f*fmid < 0){
			xr = xmid;
		}
		else{
			xl = xmid;
			f = fmid;
		}
	}

	x = xmid;
	max_iter = k;
	eps = dx;
	
	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double x1, x2;
	double x = 0.0;

	// 探索範囲
	x1 = -1.0;
	x2 = 1.0;

	// 二分法でf(x)=0を解く
	int max_iter = 100;
	double eps = 1e-10;
	bisection(Func1, x1, x2, x, max_iter, eps);

	cout.precision(15);
	// 結果の画面表示
	cout << "x = " << x << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	cout << log2(2.0/1.0e-6) << endl;

	return 0;
}


