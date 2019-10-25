/*! 
  @file conjugate-gradient.cpp
	
  @brief 共役勾配法(conjugate gradient method)
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// (前処理付き)共役勾配法(CG法)による線形システムソルバ
//-----------------------------------------------------------------------------
/*!
 * 修正コレスキー分解(modified Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii * d_i = 1とした場合 <- この部分がcholesky.cppと異なるので注意
 *  - L: i > jの要素が非ゼロで対角成分は1
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int ModifiedCholeskyDecomp2(const vector< vector<double> > &A, vector< vector<double> > &L, vector<double> &d, int n)
{
	if(n <= 0) return 1;

	L[0][0] = A[0][0];
	d[0] = 1.0/L[0][0];

	for(int i = 1; i < n; ++i){
		for(int j = 0; j <= i; ++j){
			double lld = A[i][j];
			for(int k = 0; k < j; ++k){
				lld -= L[i][k]*L[j][k]*d[k];
			}
			L[i][j] = lld;
		}
		d[i] = 1.0/L[i][i];
	}

	return 0;
}

/*!
 * 不完全コレスキー分解(incomplete Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii * d_i = 1とした場合 <- この部分がcholesky.cppと異なるので注意
 *  - L: i > jの要素が非ゼロで対角成分は1
 *  - 行列Aの値が0である要素に対応する部分を飛ばす
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int IncompleteCholeskyDecomp2(const vector< vector<double> > &A, vector< vector<double> > &L, vector<double> &d, int n)
{
	if(n <= 0) return 1;

	L[0][0] = A[0][0];
	d[0] = 1.0/L[0][0];

	for(int i = 1; i < n; ++i){
		for(int j = 0; j <= i; ++j){
			if(fabs(A[i][j]) < 1.0e-10) continue;

			double lld = A[i][j];
			for(int k = 0; k < j; ++k){
				lld -= L[i][k]*L[j][k]*d[k];
			}
			L[i][j] = lld;
		}

		d[i] = 1.0/L[i][i];
	}

	return 0;
}




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
int CGSolver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;
	cout.precision(8);

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

		// 確認のため現在の解を画面出力
		cout << k << " : ";
		for(int i = 0; i < n; ++i) cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");

		// (r*r)_(k+1)の計算
		rr1 = dot(r, r, n);

		// 収束判定 (||r||<=eps)
		e = sqrt(rr1);
		//cout << ", " << e/n;
		cout << endl;
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
 * (LDL^T)^-1 r の計算
 * @param[in] L,d IC分解で得られた下三角行列と対角行列(対角成分のみのベクトル)
 * @param[in] r 残差ベクトル
 * @param[in] u (LDL^T)^-1 rを計算した結果
 */
inline void ICRes(const vector< vector<double> > &L, const vector<double> &d, const vector<double> &r, vector<double> &u, int n)
{
	vector<double> y(n);
	for(int i = 0; i < n; ++i){
		double rly = r[i];
		for(int j = 0; j < i; ++j){
			rly -= L[i][j]*y[j];
		}
		y[i] = rly/L[i][i];
	}

	for(int i = n-1; i >= 0; --i){
		double lu = 0.0;
		for(int j = i+1; j < n; ++j){
			lu += L[j][i]*u[j];
		}
		u[i] = y[i]-d[i]*lu;
	}
}


/*!
 * 不完全コレスキー分解による前処理付共役勾配法によりA・x=bを解く
 * @param[in] A n×n正値対称行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 0:成功,1:失敗
 */
int ICCGSolver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;
	cout.precision(8);

	vector<double> r(n), p(n), y(n), r2(n);

	// 初期値を設定
	x.assign(n, 0.0);

	// コレスキー分解
	vector<double> d(n);
	vector< vector<double> > L(n, vector<double>(n, 0.0));
	//ModifiedCholeskyDecomp2(A, L, d, n);
	IncompleteCholeskyDecomp2(A, L, d, n);

	// 第0近似解に対する残差の計算
	for(int i = 0; i < n; ++i){
		double ax = 0.0;
		for(int j = 0; j < n; ++j){
			ax += A[i][j]*x[j];
		}
		r[i] = b[i]-ax;
	}

	// p_0 = (LDL^T)^-1 r_0 の計算
	ICRes(L, d, r, p, n);

	float rr0 = dot(r, p, n), rr1;
	float alpha, beta;


	double e = 0.0;
	int k;
	for(k = 0; k < max_iter; ++k){
		cout << k;

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

		// 確認のため現在の解を画面出力
		//cout << k << " : ";
		//for(int i = 0; i < n; ++i) cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");
		//cout << endl;

		// (r*r)_(k+1)の計算
		ICRes(L, d, r, r2, n);
		rr1 = dot(r, r2, n);

		// 収束判定 (||r||<=eps)
		e = sqrt(rr1);
		cout << ", " << e/n;
		cout << endl;
		if(e < eps){
			k++;
			break;
		}

		// βの計算とPの更新
		beta = rr1/rr0;
		for(int i = 0; i < n; ++i){
			p[i] = r2[i]+beta*p[i];
		}

		// (r*r)_(k+1)を次のステップのために確保しておく
		rr0 = rr1;
	}

	max_iter = k;
	eps = e;

	return 0;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込んで解く
	vector< vector<double> > Ab, A;
	ReadMatrix("matrix.txt", ",", Ab);

	int n = (int)Ab.size();	// n元連立一次方程式

	// 読み込んだ行列を確認用に画面表示
	cout << "A(" << n << " x " << n+1 << ") = " << endl;
	OutputMatrix(Ab, n, n+1);
	cout << endl;

	// 左辺の行列の抽出
	A.resize(n, vector<double>(n));
	for(int i = 0; i < n; ++i){
		A[i].assign(Ab[i].begin(), Ab[i].begin()+n);
	}

	// 右辺項の抽出
	vector<double> b(n);
	for(int i = 0; i < n; ++i){
		b[i] = Ab[i][n];
	}
	
	// 共役勾配法を用いて連立1次方程式を解く
	int max_iter = 100;
	double eps = 1e-6;
	vector<double> x(n);
	//CGSolver(A, b, x, n, max_iter, eps);
	ICCGSolver(A, b, x, n, max_iter, eps);
	
	// 結果の画面表示
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");
	}
	cout << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;



	return 0;
}


