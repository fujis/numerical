/*! 
  @file spline.cpp
	
  @brief スプライン補間
 
  @author Makoto Fujisawa
  @date 2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

#include <iomanip>



/*!
 * 共役勾配法によりA・x=bを解く
 * @param[in] A n×n正値対称行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 1:成功,0:失敗
 */
int cg_solver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 0;

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

	return 1;
}

//-----------------------------------------------------------------------------
// 補間法
//-----------------------------------------------------------------------------
/*!
 * スプライン補間
 *  - 
 * @param[in] yi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
 * @param[in] n データ数
 * @param[in] x 補間した値が必要な位置x
 * @param[out] ans 解
 * @return
 */
int spline_interpolation(const vector<double> &yi, const vector<double> &xi, int n, double x, double &ans)
{
	int m = n-1; // 補間区間の数
	vector<double> a(m+1, 0.0), b(m, 0.0), c(m+1, 0.0), d(m, 0.0); // 各区間における補間係数(計算上aとcだけサイズがm+1になっていることに注意)
	vector<double> h(m); // 各補間区間の幅
	for(int i = 0; i < m; ++i){
		h[i] = xi[i+1]-xi[i];
	}

	// 係数aの計算(0次項の係数=定数項)
	for(int i = 0; i < m+1; ++i){
		a[i] = yi[i];
	}

	// 係数cの計算(2次項の係数) -> 線形システムを解くことで求める
	// 線形システムの係数行列と右辺項ベクトルの計算
	vector< vector<double> > Ac(m+1);
	vector<double> bc(m+1, 0.0);
	for(int i = 0; i < m+1; ++i){
		Ac[i].resize(m+1, 0.0); // 係数行列の各行を0で初期化
		if(i == 0){
			Ac[i][0] = 1.0; bc[0] = 0.0;
		}
		else if(i == m){
			Ac[i][m] = 1.0; bc[m] = 0.0;
		}
		else{
			Ac[i][i-1] = h[i-1];
			Ac[i][i]   = 2*(h[i-1]+h[i]);
			Ac[i][i+1] = h[i];
			bc[i] = 3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1];
		}
	}
	// CG法で線形システムを解く(Acは対称疎行列)
	int max_iter = 100;
	double eps = 1e-6;
	cg_solver(Ac, bc, c, m+1, max_iter, eps);

	// 係数b,dの計算(1次項,3次項の係数)
	for(int i = 0; i < m-1; ++i){
		b[i] = (a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2*c[i])/3;
		d[i] = (c[i+1]-c[i])/(3*h[i]);
	}
	b[m-1] = d[m-1] = 0.0;

	// xが含まれる区間の探索(区間幅がすべて同じならint(x/h)で求められる)
	int k = 0;
	for(int i = 0; i < m; ++i){
		if(x >= xi[i] && x < xi[i+1]){
			k = i; break;
		}
	}

	// 計算済みの係数を使って3次多項式で補間値を計算
	ans = a[k] + b[k]*(x-xi[k]) + c[k]*(x-xi[k])*(x-xi[k]) + d[k]*(x-xi[k])*(x-xi[k])*(x-xi[k]);

	return 0;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncExp;
	//double(*func)(double) = FuncLinear;

	double x = 0.5;
	double fx;
	double gt = func(x); // 真値
	cout.precision(10);

	// 補間用の点
	vector<double> xi, yi;
	xi.push_back(0.0);
	yi.push_back(func(xi.back()));
	xi.push_back(0.2);
	yi.push_back(func(xi.back()));
	xi.push_back(0.4);
	yi.push_back(func(xi.back()));
	xi.push_back(0.6);
	yi.push_back(func(xi.back()));
	xi.push_back(0.8);
	yi.push_back(func(xi.back()));
	xi.push_back(1.0);
	yi.push_back(func(xi.back()));

	// 線形補間
	spline_interpolation(yi, xi, xi.size(), x, fx);
	cout << "f_spline(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	cout << "ground truth = " << gt << endl;

		
	return 0;
}


