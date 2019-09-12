/*! 
  @file gradientdecent.cpp
	
  @brief 最急降下法
 
  @author Makoto Fujisawa
  @date 2019-07
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// 最適化関数(最小値探索)
//-----------------------------------------------------------------------------
/*!
 * 最急降下法(gradient decent method, steepest decent method)
 * @param[in] dfunc 関数勾配値を与える関数ポインタ
 * @param[in] x0 初期探索地点
 * @param[out] xout 解
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */
int gradientdecent(double dfunc(const double), double x0, double alpha, double &xout, int &max_iter, double &eps)
{
	double x = x0, dx = 0.0;

	int k;
	for(k = 0; k < max_iter; ++k){
		cout << "x(" << k << ") = " << x << endl;
		dx = alpha*dfunc(x);	// 関数の勾配に係数αを掛ける
		x -= dx;				// 勾配方向に探索点を移動させる
		if(fabs(dx) < eps) break;	// 勾配の大きさで収束判定
	}

	xout = x;
	max_iter = k;
	eps = fabs(dx);
	
	return 0;
}

/*!
 * 最急降下法(gradient decent method, steepest decent method)
 *  - n次元(n>1)版
 * @param[in] dfunc 関数勾配値を与える関数ポインタ
 * @param[in] x0 初期探索地点
 * @param[out] xout 解
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 1:成功,0:失敗
 */
int gradientdecent(vector<double> dfunc(const vector<double>&), vector<double> x0, double alpha, vector<double> &xout, int &max_iter, double &eps)
{
	int n = x0.size();
	vector<double> x = x0, dx(n, 0.0);

	double norm_dx; // 勾配ベクトルのノルム(収束判定用)
	int k;
	for(k = 0; k < max_iter; ++k){
		
		dx = dfunc(x);
		norm_dx = 0.0;
		for(int i = 0; i < n; ++i){
			dx[i] *= alpha;		// 関数の勾配ベクトルに係数αを掛ける
			x[i] -= dx[i];		// 勾配方向に探索点を移動させる
			norm_dx += dx[i]*dx[i];
		}

		// 勾配ベクトルのノルムで収束判定
		norm_dx = sqrt(norm_dx);
		if(norm_dx < eps) break;
	}

	xout = x;
	max_iter = k;
	eps = norm_dx;

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double alpha = 0.2;
	int max_iter = 100;
	double eps = 1e-6;

#if 0
	double x0 = 1.0;
	double x = 0.0;
	gradientdecent(DFunc3, x0, alpha, x, max_iter, eps);
	cout << "x = " << x << endl;
#else
	vector<double> x0(2, 2.0);
	vector<double> x(2);
	gradientdecent(DFunc4, x0, alpha, x, max_iter, eps);

	cout << "x,y = " << x[0] << ", " << x[1] << endl;

#endif


	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;


	return 0;
}


