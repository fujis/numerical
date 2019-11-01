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
 * 3次スプライン補間
 * @param[in] f 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
 * @param[in] m データ数(=n+1)
 * @param[in] x 補間した値が必要な位置x
 * @param[out] ans 解
 * @param[out] a,b,c,d 各区間における補間係数
 * @return
 */
int spline_interpolation(const vector<double> &f, const vector<double> &xi, int m, double x, double &ans,
						 vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d)
{
	int n = m-1; // 補間区間の数
	// 補間係数配列のメモリ確保
	a.resize(n); b.resize(n); c.resize(n); d.resize(n);
	vector<double> h(n); // 各補間区間の幅
	for(int i = 0; i < n; ++i){
		h[i] = xi[i+1]-xi[i];
	}

	// 位置x_iでの多項式の2階微分u_iの計算
	// 線形システムの係数行列Hと右辺項ベクトルbの計算
	vector< vector<double> > H(n-1);
	vector<double> tmp(n-1, 0.0), v(n-1, 0.0);
	for(int i = 0; i < n-1; ++i){
		H[i].resize(n-1, 0.0); // 係数行列の各行を0で初期化
		// 係数行列と右辺項ベクトルの要素の計算
		if(i != 0) H[i][i-1] = h[i];
		H[i][i] = 2*(h[i]+h[i+1]);
		if(i != n-2) H[i][i+1] = h[i+1];
		v[i] = 6*((f[i+2]-f[i+1])/h[i+1]-(f[i+1]-f[i])/h[i]);
	}
	// CG法で線形システムを解く(Hは対称疎行列)
	int max_iter = 100;
	double eps = 1e-6;
	cg_solver(H, v, tmp, n-1, max_iter, eps);

	// u_0とu_nを追加
	vector<double> u(n+1);
	u[0] = 0.0; u[n] = 0.0;
	for(int i = 0; i < n-1; ++i) u[i+1] = tmp[i];

	// 係数a,b,c,dの計算
	for(int i = 0; i < n; ++i){
		a[i] = (u[i+1]-u[i])/(6.0*h[i]);
		b[i] = u[i]/2.0;
		c[i] = (f[i+1]-f[i])/h[i]-h[i]*(2*u[i]+u[i+1])/6.0;
		d[i] = f[i];
	}

	// xが含まれる区間の探索(区間幅がすべて同じならint(x/h)で求められる)
	int k = 0;
	for(int i = 0; i < n; ++i){
		if(x >= xi[i] && x < xi[i+1]){
			k = i; break;
		}
	}

	// 計算済みの係数を使って3次多項式で補間値を計算
	double dx = x-xi[k];
	ans = a[k]*dx*dx*dx + b[k]*dx*dx + c[k]*dx + d[k];

	return 0;
}


/*!
 * 3次スプライン多項式で関数値を計算
 * @param[in] n 区間数
 * @param[in] x 補間した値が必要な位置
 * @param[out] a,b,c,d 各区間における補間係数
 * @return スプライン多項式を計算した結果
 */
double cubic_spline(double x, const vector<double> &xi, int n, const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d)
{
	// xが含まれる区間の探索(区間幅がすべて同じならint(x/h)で求められる)
	int k = 0;
	for(int i = 0; i < n; ++i){
		if(x >= xi[i] && x < xi[i+1]){
			k = i; break;
		}
	}

	// 計算済みの係数を使って3次多項式で補間値を計算
	double dx = x-xi[k];
	return a[k]*dx*dx*dx + b[k]*dx*dx + c[k]*dx + d[k];
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncExp;
	double x0 = 0, x1 = 1;

	//double(*func)(double) = FuncRunge;
	//double x0 = -1, x1 = 1;

	double x = 0.5, fx;
	double gt = func(x); // 真値
	cout.precision(6);

	// 補間用の点(サンプリング点数6)
	vector<double> xi, yi;
	MakeSamplingPoints(x0, x1, (x1-x0)/5, func, xi, yi);
	OutputSamplingPoints(xi, yi); // サンプリング点の画面出力

	// スプライン補間
	vector<double> a, b, c, d;
	spline_interpolation(yi, xi, xi.size(), x, fx, a, b, c, d);
	cout << "f_spline(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	cout << "ground truth = " << gt << endl;



	// グラフ描画用にデータ出力
	int m = 11; // データ点数
	MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func, xi, yi);
	spline_interpolation(yi, xi, xi.size(), x, fx, a, b, c, d);
	OutputSamplingPoints(xi, yi, "dat/spline"+TOSTR(m)+"_data.txt"); // サンプリング点のファイル出力
	OutputFunction(x0, x1, (x1-x0)/200, std::bind(cubic_spline, std::placeholders::_1, xi, xi.size()-1, a, b, c, d), "dat/spline"+TOSTR(m)+".txt");

	// 真値のグラフ作成用ファイル出力
	OutputFunction(x0, x1, (x1-x0)/200, func, "dat/spline_ground_truth.txt");
	
	return 0;
}


