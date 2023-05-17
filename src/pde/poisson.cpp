/*! 
  @file poisson.cpp
	
  @brief ポアソン方程式のソルバー
 
  @author Makoto Fujisawa
  @date 2019-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

// 境界値
double alpha = 0.0, beta = 0.0;

//-----------------------------------------------------------------------------
//! ポアソン方程式のソルバー
//-----------------------------------------------------------------------------

/*!
 * 共役勾配法によりA・x=bを解く
 * @param[in] A n×n正値対称行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 0:成功,1:失敗
 */
int cg_solver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;

	vector<double> r(n), p(n), y(n);
	x.assign(n, 0.0);

	// 第0近似解に対する残差の計算
	for(int i = 0; i < n; ++i){
		double ax = 0.0;
		for(int j = 0; j < n; ++j){
			ax += A[i][j]*x[j];
		}
		r[i] = b[i]-ax;
		p[i] = r[i];
	}

	double rr0 = dot(r, r, n), rr1;
	double alpha, beta;

	double e = 0.0;
	int k;
	for(k = 0; k < max_iter; ++k){
		//cout << k;

		// y = AP の計算
		for(int i = 0; i < n; ++i){
			y[i] = dot(A[i], p, n);
		}

		// alpha = r*r/(P*AP)の計算
		alpha = rr0/dot(p, y, n);

		// 解x、残差rの更新
		for(int i = 0; i < n; ++i){
			x[i] += alpha*p[i];
			r[i] -= alpha*y[i];
		}

		// (r*r)_(k+1)の計算
		rr1 = dot(r, r, n);

		// 収束判定 (||r||<=eps)
		e = sqrt(rr1);
		if(e < eps){
			k++;
			break;
		}

		// βの計算とPの更新
		beta = rr1/rr0;
		for(int i = 0; i < n; ++i){
			p[i] = r[i]+beta*p[i];
		}

		// (r*r)_(k+1)を次のステップのために確保しておく
		rr0 = rr1;
	}

	max_iter = k+1;
	eps = e;

	return 0;
}

/*!
 * 境界条件設定(ディリクレ境界条件)
 * @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 */
void setbc_d(vector<double> &f, int n)
{
	f[0] = alpha;
	f[n] = beta;
}

/*!
 * 中心差分で1次元ポアソン方程式を解く
 * @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
 * @param[in] x0,x1 計算範囲
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 * @param[in] g 右辺項を与える関数ポインタ
 * @return 
 */
int poisson1d_central(vector<double> &f, double x0, double x1, int n, double g(double))
{
	double h = (x1-x0)/n; // 空間刻み幅

	// 線形システムの係数行列と右辺項ベクトルの計算
	vector<double> b(n-1, 0.0), x(n-1, 0.0);
	vector< vector<double> > A(n-1, b); // 各要素をbで初期化することでn-1×n-1の2次元配列になる
	for(int i = 1; i < n-2; ++i){
		double xi = x0+h*(i+1);
		A[i][i-1] = -1; A[i][i] = 2; A[i][i+1] = -1;
		b[i] = -h*h*g(xi);
	}
	A[0][0] = 2; A[0][1] = -1; b[0] = -h*h*g(x0+h)+f[0];
	A[n-2][n-3] = -1; A[n-2][n-2] = 2; b[n-2] = -h*h*g(x1-h)+f[n];
	for(int i = 0; i < n-1; ++i) x[i] = f[i+1];


	// CG法で線形システムを解く
	int max_iter = 100;
	double eps = 1e-6;
	cg_solver(A, b, x, n-1, max_iter, eps);
	cout << "cg solver : max_iter = " << max_iter << ", eps = " << eps << endl;

	// 結果を配列fに戻す
	for(int i = 0; i < n-1; ++i) f[i+1] = x[i];
	setbc_d(f, n); // 境界条件

	return 1;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncPdeX3;
	double a = 0.0, b = 1.0;
	double alpha = 0.0, beta = 0.0;
	int n = 10;
	double(*func_t)(double) = FuncPdeX3T;

	vector<double> f(n+1, 0.0);
	setbc_d(f, n);

	// 中心差分+CG法でポアソン方程式を解く
	poisson1d_central(f, a, b, n, func);

	// 結果の出力
	double h = (b-a)/n;
	double avg_error = 0.0;
	for(int i = 0; i < n+1; ++i){
		double x = a+h*i;

		double e = fabs(f[i]-func_t(x));
		cout << "f(" << x << ") = " << f[i] << ", error = " << e << endl;
		//cout << x << ", " << f[i] << ", " << func_t(x) << endl;
		avg_error += e;
	}

	cout << "avg. error = " << avg_error/(n-1) << endl;



	return 0;
}


