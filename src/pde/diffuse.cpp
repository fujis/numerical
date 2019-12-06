/*! 
  @file diffuse.cpp
	
  @brief 拡散方程式のソルバー
 
  @author Makoto Fujisawa
  @date 2019-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


// 差分法用関数ポインタ
typedef int(*FUNC)(vector<double> &, vector<double> &, double, double, int, double, double);
typedef int(*FUNC2)(vector< vector<double> > &, vector< vector<double> > &, double, double, int, double, double, double, double);
// 境界値設定用関数ポインタ
typedef void(*BC)(vector<double> &, int);
typedef void(*BC2)(vector< vector<double> > &, int);


// 境界値
double alpha = 0.0, beta = 0.0;
double cn_lambda = 0.5;

//-----------------------------------------------------------------------------
//! 拡散方程式のソルバー
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
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 */
void setbc_d(vector<double> &f, int n)
{
	f[0] = alpha;
	f[n] = beta;
}

/*!
 * 境界条件設定(ディリクレ境界条件)
 * @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 */
void setbc2_d(vector< vector<double> > &f, int n)
{
	for(int i = 0; i < n; ++i){ f[i][0] = alpha; f[i][n] = beta; }
	for(int j = 0; j < n; ++j){ f[0][j] = alpha; f[n][j] = beta; }
}

/*!
 * FTCS法(前進オイラー法+中心差分)による拡散方程式のソルバー
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 * @return 問題なければ0を返す
 */
int diffuse1d_ftcs(vector<double> &f1, vector<double> &f0, double a, double dt, int n, double x0, double xn)
{
	if(n < 0) return 1;
	double dx = (xn-x0)/n; // 空間刻み幅
	double eta = a*dt/(dx*dx);
	for(int i = 1; i <= n-1; ++i){
		f1[i] = f0[i]+eta*(f0[i+1]-2*f0[i]+f0[i-1]);
	}

	return 0;
}

/*!
 * クランク・ニコルソン法による拡散方程式のソルバー
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 * @return 問題なければ0を返す
 */
int diffuse1d_cn(vector<double> &f1, vector<double> &f0, double a, double dt, int n, double x0, double xn)
{
	if(n < 0) return 1;
	double dx = (xn-x0)/n; // 空間刻み幅

	double eta = a*dt/(dx*dx);

	// 線形システムの係数行列と右辺項ベクトルの計算
	vector<double> b(n-1, 0.0), x(n-1, 0.0);
	vector< vector<double> > A(n-1, b); // 各要素をbで初期化することでn-1×n-1の2次元配列になる
	for(int i = 0; i < n-1; ++i){
		int j = i+1; // 行列としてはi行目でf[i+1]での式となる
		double xi = a+dx*j;
		if(i != 0) A[i][i-1] = -cn_lambda*eta;
		A[i][i] = 1+2*cn_lambda*eta;
		if(i != n-2) A[i][i+1] = -cn_lambda*eta;
		b[i] = f0[j]+(1-cn_lambda)*eta*(f0[j+1]-2*f0[j]+f0[j-1]);
	}
	b[0] += cn_lambda*eta*f0[0];
	b[n-2] += cn_lambda*eta*f0[n];
	for(int i = 0; i < n-1; ++i) x[i] = f0[i+1];

	// CG法で線形システムを解く
	int max_iter = 100;
	double eps = 1e-6;
	cg_solver(A, b, x, n-1, max_iter, eps);
	//cout << "cg solver : max_iter = " << max_iter << ", eps = " << eps << endl;

	// 結果を配列fに戻す
	for(int i = 0; i < n-1; ++i) f1[i+1] = x[i];

	return 0;
}


/*!
 * FTCS法(前進オイラー法+中心差分)による拡散方程式のソルバー
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
 * @return 問題なければ0を返す
 */
int diffuse2d_ftcs(vector< vector<double> > &f1, vector< vector<double> > &f0, double a, double dt, int n, double x0, double xn, double y0, double yn)
{
	if(n < 0) return 1;
	double dx = (xn-x0)/n; // 空間刻み幅
	double dy = (yn-y0)/n; // 空間刻み幅

	double eta1 = a*dt/(dx*dx);
	double eta2 = a*dt/(dy*dy);
	for(int i = 1; i <= n-1; ++i){
		for(int j = 1; j <= n-1; ++j){
			f1[i][j] = f0[i][j]+eta1*(f0[i+1][j]-2*f0[i][j]+f0[i-1][j])+eta2*(f0[i][j+1]-2*f0[i][j]+f0[i][j-1]);
		}
	}

	return 0;
}

/*!
 * クランク・ニコルソン法による拡散方程式のソルバー
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
 * @return 問題なければ0を返す
 */
int diffuse2d_cn(vector< vector<double> > &f1, vector< vector<double> > &f0, double a, double dt, int n, double x0, double xn, double y0, double yn)
{
	if(n < 0) return 1;
	double dx = (xn-x0)/n; // 空間刻み幅
	double dy = (yn-y0)/n; // 空間刻み幅

	double eta1 = a*dt/(dx*dx);
	double eta2 = a*dt/(dy*dy);

	// 線形システムの係数行列と右辺項ベクトルの計算
	vector<double> b((n-1)*(n-1), 0.0), x((n-1)*(n-1), 0.0);
	vector< vector<double> > A((n-1)*(n-1), b); // 各要素をbで初期化することで(n-1)^2×(n-1)^2の2次元配列になる
	for(int i = 0; i < n-1; ++i){
		int i1 = i+1;
		for(int j = 0; j < n-1; ++j){
			int j1 = j+1;
			int k = i+j*(n-1);

			A[k][k] = 1+2*cn_lambda*eta1+2*cn_lambda*eta2;

			if(i != 0  ) A[k][(i-1)+(j)*(n-1)] = -cn_lambda*eta1;
			if(i != n-2) A[k][(i+1)+(j)*(n-1)] = -cn_lambda*eta1;
			if(j != 0  ) A[k][(i)+(j-1)*(n-1)] = -cn_lambda*eta2;
			if(j != n-2) A[k][(i)+(j+1)*(n-1)] = -cn_lambda*eta2;

			b[k] = f0[i1][j1]+(1-cn_lambda)*eta1*(f0[i1+1][j1]-2*f0[i1][j1]+f0[i1-1][j1])+(1-cn_lambda)*eta2*(f0[i1][j1+1]-2*f0[i1][j1]+f0[i1][j1-1]);
		}
	}
	for(int j = 0; j < n-1; ++j){
		b[(0)+j*(n-1)]   += cn_lambda*eta1*f0[0][j+1];
		b[(n-2)+j*(n-1)] += cn_lambda*eta1*f0[n][j+1];
	}
	for(int i = 0; i < n-1; ++i){
		b[i+(0)*(n-1)]   += cn_lambda*eta2*f0[i+1][0];
		b[i+(n-2)*(n-1)] += cn_lambda*eta2*f0[i+1][n];
	}
	for(int i = 0; i < n-1; ++i){
		for(int j = 0; j < n-1; ++j){
			x[i+j*(n-1)] = f0[i+1][j+1];
		}
	}

	// CG法で線形システムを解く
	int max_iter = 100;
	double eps = 1e-6;
	cg_solver(A, b, x, (n-1)*(n-1), max_iter, eps);
	//cout << "cg solver : max_iter = " << max_iter << ", eps = " << eps << endl;

	// 結果を配列fに戻す
	for(int i = 0; i < n-1; ++i){
		for(int j = 0; j < n-1; ++j){
			f1[i+1][j+1] = x[i+j*(n-1)];
		}
	}

	return 0;
}

/*!
 * 1ステップ進める
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[inout] f 未知関数fの各グリッドでの値
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 */
void diffuse1d_step(FUNC sdfunc, BC bcfunc, vector<double> &f, double a, double dt, int n, double x0, double xn)
{
	vector<double> f0(f);
	bcfunc(f0, n);
	sdfunc(f, f0, a, dt, n, x0, xn);
	bcfunc(f, n);
}

/*!
 * 時間ステップを進めていったときの結果をファイル出力するための関数
 * @param[in] tdfunc 時間差分を行う関数
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピー渡しにしている)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 * @param[in] filename 出力ファイル名
 */
void diffuse1d(FUNC sdfunc, BC bcfunc, vector<double> f, double a, double dt, int n, double x0, double xn, const string &filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	fo << "#1d,0.0,1.0" << endl;

	// 初期値の出力
	OutputValueToFile(f, n, x0, xn, 0, fo);

	// タイムステップを進める
	for(int k = 1; k < 1000; ++k){
		diffuse1d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn);

		// 結果の出力
		OutputValueToFile(f, n, x0, xn, k*dt, fo);
	}
	fo.close();

	cout << "diffuse1d : " << filename << endl;
}


/*!
 * 1ステップ進める
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピー渡しにしている)
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
 * @param[in] x0,xn,y0,yn x方向,y方向の計算範囲
 */
void diffuse2d_step(FUNC2 sdfunc, BC2 bcfunc, vector< vector<double> > &f, double a, double dt, int n, double x0, double xn, double y0, double yn)
{
	vector< vector<double> > f0(f);
	bcfunc(f0, n);
	sdfunc(f, f0, a, dt, n, x0, xn, y0, yn);
	bcfunc(f, n);
}

/*!
 * 時間ステップを進めていったときの結果をファイル出力するための関数
 * @param[in] tdfunc 時間差分を行う関数
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] a 拡散係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(dx=(xn-x0)/n,dy=(yn-y0)/n)
 * @param[in] filename 出力ファイル名
 */
void diffuse2d(FUNC2 sdfunc, BC2 bcfunc, vector< vector<double> > f, double a, double dt, int n, double x0, double xn, double y0, double yn, const string &filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	fo << "#2d,0.0,1.0" << endl;

	// 初期値の出力
	OutputValueToFile(f, n, x0, xn, y0, yn, 0, fo);

	// タイムステップを進める
	for(int k = 1; k < 200; ++k){
		diffuse2d_step(sdfunc, bcfunc, f, a, dt, n, x0, xn, y0, yn);

		// 結果の出力
		OutputValueToFile(f, n, x0, xn, y0, yn, k*dt, fo);
	}
	fo.close();

	cout << "diffuse2d : " << filename << endl;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double dt = 0.01;
	double x0 = 0.0, xn = 1.0;
	double a = 0.02;
	int n = 50;
	BC bcfunc = setbc_d;

	vector<double> f(n+1, 0.0);

	// 初期値設定
	for(int i = 0; i < n+1; ++i) f[i] = 1.0;

	// データ保存パス(ファイル名の最後の_ignはGitで除外ファイルにするためのもの)
	string path = "../../bin/data";


	// 前進オイラー+中心差分
	diffuse1d(diffuse1d_ftcs, bcfunc, f, a, dt, n, x0, xn, path+"diffuse1d_ftcs_ign.txt");
	// クランク・ニコルソン法
	cn_lambda = 0.5;
	diffuse1d(diffuse1d_cn, bcfunc, f, a, dt, n, x0, xn, path+"diffuse1d_cn_ign.txt");

	dt = 0.011;
	// 前進オイラー+中心差分
	diffuse1d(diffuse1d_ftcs, bcfunc, f, a, dt, n, x0, xn, path+"diffuse1d_ftcs_dt11_ign.txt");
	// クランク・ニコルソン法
	cn_lambda = 0.5;
	diffuse1d(diffuse1d_cn, bcfunc, f, a, dt, n, x0, xn, path+"diffuse1d_cn_dt11_ign.txt");


	// 2次元
	double y0 = 0.0, yn = 1.0;
	dt = 0.005;
	BC2 bcfunc2 = setbc2_d;
	for(int i = 0; i < n+1; ++i) f[i] = 1.0;
	vector< vector<double> > f2(n+1, f);

	// 前進オイラー+中心差分
	diffuse2d(diffuse2d_ftcs, bcfunc2, f2, a, dt, n, x0, xn, y0, yn, path+"diffuse2d_ftcs_ign.txt");
	// クランク・ニコルソン法
	cn_lambda = 0.5;
	diffuse2d(diffuse2d_cn, bcfunc2, f2, a, dt, n, x0, xn, y0, yn, path+"diffuse2d_cn_ign.txt");






	return 0;
}


