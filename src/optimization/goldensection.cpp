/*! 
  @file goldensection.cpp
	
  @brief 黄金分割法
 
  @author Makoto Fujisawa
  @date 2019-07
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

const double SQRT5 = 2.23606797749979;

//-----------------------------------------------------------------------------
// 最適化関数(最小値探索)
//-----------------------------------------------------------------------------
/*!
 * 黄金分割探索法(golden-section method)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] xl,xr 初期探索範囲
 * @param[out] x 解
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */
int goldensection(double func(const double), double xl, double xr, double &ans, int &max_iter, double &eps)
{
	const double e = (SQRT5-1.0)/(SQRT5+1.0);
	//const double e = 2.0/(SQRT5+3.0);	// 上の式の変形(√5を事前計算してないならこちらの方が速い)
	double b = xr-xl;
	double d = e*b;
	cout.precision(7);

	double x[4];
	x[0] = xl;
	x[1] = xl+d;
	x[2] = xr-d;
	x[3] = xr;

	// x[1],x[2]における関数値
	double f1 = func(x[1]);
	double f2 = func(x[2]);

	int k;
	for(k = 0; k < max_iter; ++k){

		cout << k << " : " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << ", b = " << b << endl;

		if(f1 < f2){ // 区間[x2,x3]に最小値なし
			double x2 = x[2];
			x[2] = x[1];
			x[1] = x2-d;
			x[3] = x2;
			f2 = f1;
			f1 = func(x[1]);
		}
		else{ // 区間[x0,x1]に最小値なし
			double x1 = x[1];
			x[1] = x[2];
			x[2] = x1+d;
			x[0] = x1;
			f1 = f2;
			f2 = func(x[2]);
		}

		b -= d;	// 新しい[x0,x3]の長さ(現在の長さから[x0,x1](or[x2,x3])の長さdを引く)
		d = e*b; // 新しい[x0,x1](or[x2,x3])の長さ

		if(b < eps) break;
	}

	ans = 0.5*(x[1]+x[2]);
	max_iter = k;
	eps = b;
	
	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double x1, x2;
	double x = 0.0;

	x1 = 0.0;
	x2 = 2.0;

	int max_iter = 100;
	double eps = 1e-6;
	goldensection(Func2, x1, x2, x, max_iter, eps);

	cout << "x = " << x << endl;
	cout << "f(x,y) = " << Func2(x) << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;


	return 0;
}


