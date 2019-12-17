/*!
  @file qr.cpp

  @brief 固有値・固有ベクトル
		 QR法

  @author Makoto Fujisawa
  @date 2012-06
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "qr.h"


//-----------------------------------------------------------------------------
// QR法による固有値計算
//-----------------------------------------------------------------------------
/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
 *  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
 *  - LとUを一つの行列にまとめた形で結果を返す
 * @param[inout] A n×nの係数行列．LU分解した結果を格納する．
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int lu_decomp(vector< vector<double> > &A, int n)
{
	if(n <= 0) return 1;

	for(int i = 0; i < n; ++i){
		// l_ijの計算(i >= j)
		for(int j = 0; j <= i; ++j){
			double lu = A[i][j];
			for(int k = 0; k < j; ++k){
				lu -= A[i][k]*A[k][j];	// l_ik * u_kj
			}
			A[i][j] = lu;
		}

		// u_ijの計算(i < j)
		for(int j = i+1; j < n; ++j){
			double lu = A[i][j];
			for(int k = 0; k < i; ++k){
				lu -= A[i][k]*A[k][j];	// l_ik * u_kj
			}
			A[i][j] = lu/A[i][i];
		}
	}

	return 0;
}



/*!
 * LU分解した行列a(n×n)から前進代入・後退代入によりA・x=bを解く
 * @param[in] A LU分解された行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int lu_solver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n)
{
	if(n <= 0) return 1;

	// 前進代入(forward substitution)
	//  LY=bからYを計算
	for(int i = 0; i < n; ++i){
		double bly = b[i];
		for(int j = 0; j < i; ++j){
			bly -= A[i][j]*x[j];
		}
		x[i] = bly/A[i][i];
	}

	// 後退代入(back substitution)
	//  UX=YからXを計算
	for(int i = n-1; i >= 0; --i){
		double yux = x[i];
		for(int j = i+1; j < n; ++j){
			yux -= A[i][j]*x[j];
		}
		x[i] = yux;
	}

	return 0;
}


/*!
 * ハウスホルダー変換でQR分解 (A=QR)
 * @param[in] A 元の行列(n×n)
 * @param[in] Q 分解された直行行列(n×n)
 * @param[in] R 分解された上三角行列(n×n)
 * @param[in] n 行列のサイズ
 */
void qr_decomposition(const vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R, int n)
{
	vector<double> u(n, 0.0);
	vector< vector<double> > H(n, u); // 行列H
	R = A; // A^(0)
	Q.resize(n, u); unit(Q, n); // Qを単位行列で初期化

	for(int k = 0; k < n-1; ++k){
		// a_kk^(k+1)の計算 (ベクトルxの大きさ)
		double akk = 0.0;
		for(int i = k; i < n; ++i) akk += R[i][k]*R[i][k];
		akk = sqrt(akk);

		// u=(x1-y1)/|x1-y1|の計算
		double l = 0.0;
		for(int i = k; i < n; ++i){
			u[i] = R[i][k]-(i == k ? akk : 0.0);
			l += u[i]*u[i];
		}
		if(fabs(l) < 1e-10) break;
		l = sqrt(l);
		for(int i = k; i < n; ++i) u[i] /= l;

		// ハウスホルダー行列H^(k)の計算(H=I-2uu^T)
		unit(H, n); // Hを単位行列で初期化
		for(int i = k; i < n; ++i){
			for(int j = k; j < n; ++j){
				H[i][j] -= 2*u[i]*u[j];
			}
		}

		// A^(k+1) = H^(k) A^(k)
		R = mul_mm(H, R, n);

		// Q^(k+1) = Q^(k) (H^(k))^T
		Q = mul_mm(Q, transpose(H, n), n);
	}
}

/*!
 * 逆反復法による固有ベクトルの算出
 * @param[inout] A 元の行列
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 反復回数
 */
int inverse_iteration(const vector< vector<double> > &A, const vector<double> lambda, vector< vector<double> > &v, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;

	vector< vector<double> > LU(A); // LU分解した行列を格納するための配列
	vector<double> y(n, 1.0), y1(n); // y^(k), y^(k+1)格納用
	double sum_e = 0.0;
	int sum_k = 0;
	v.resize(n);

	for(int l = 0; l < n; ++l){
		// (A-λI)を計算
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				LU[i][j] = A[i][j]-(i == j ? lambda[l] : 0.0);
			}
		}

		// (A-λI)をLU分解
		lu_decomp(LU, n);

		for(int i = 0; i < n; ++i) y[i] = 1.0;
		double l2 = dot(y, y, n);
		double mu, mu0 = 0.0;
		double e;
		int k;
		for(k = 0; k < max_iter; ++k){
			// |y^(k)|=1となるように正規化
			double ly = sqrt(l2); // |y^(k)|の計算
			y = mul_sv(1.0/ly, y, n); // |y^(k)|で割る(1/|y|を掛ける)

			// y^(k+1) = B y^(k)の計算(By^(k+1)=y^(k)をLU分解で解く)
			lu_solver(LU, y, y1, n);

			// 固有値1/|λ-λi|の計算(|y^(k)|=1で正規化されていることが前提)
			mu = dot(y, y1, n);

			// 収束判定
			if((e = fabs((mu-mu0)/mu)) < eps) break;

			y = y1;
			l2 = dot(y, y, n);
			mu0 = mu;
		}
		v[l] = mul_sv(1.0/sqrt(dot(y, y, n)), y, n); // 最後に正規化しておく
		sum_e += e;
		sum_k += k;
	}


	max_iter = sum_k/n;
	eps = sum_e/n;


	return 0;
}

/*!
 * QR法による固有値の算出
 * @param[inout] A 元の行列．計算後，対角要素に固有値が入る
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 反復回数
 */
int qr(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;
	vector< vector<double> > R(A), Q(A); // R,Q用配列を確保(Aと同じ大きさに初期化)
	vector<double> lambda(n, 0.0); // 収束判定用
	double e = eps*n;
	int k;
	for(k = 0; k < max_iter; ++k){
		// A^(k) ⇐ Q^(k) R^(k)
		qr_decomposition(A, Q, R, n);

		// A^(k+1) ⇐ R^(k) Q^(k)
		A = mul_mm(R, Q, n);

		// 確認用表示
		//cout << k << " : lambda = ";
		//for(int i = 0; i < n; ++i) cout << A[i][i] << (i == n-1 ? "" : ", ");
		//cout << endl;

		// 収束判定
		e = 0;
		for(int i = 0; i < n; ++i) e += fabs(lambda[i]-A[i][i]);
		if(e/n < eps) break;

		for(int i = 0; i < n; ++i) lambda[i] = A[i][i];
	}
	eps = e/n;
	max_iter = k;

	return 0;
}


/*!
 * 2次元座標データから共分散行列を作成
 * @param[in] data 2次元座標データ
 * @param[out] A 2x2の共分散行列
 */
void covariance_matrix2(const vector<rxPoint2> &data, vector< vector<double> > &A, rxPoint2 &c)
{
	int m = data.size();
	int n = 2;

	// データの平均値(重心)の計算
	c.x = c.y = 0.0;
	for(rxPoint2 p : data){
		c.x += p.x; c.y += p.y;
	}
	c.x /= m;
	c.y /= m;

	// 共分散行列の計算
	A.resize(n);
	for(int i = 0; i < A.size(); ++i) A[i].resize(n, 0.0);
	for(rxPoint2 p : data){
		A[0][0] += (p.x-c.x)*(p.x-c.x);
		A[1][1] += (p.y-c.y)*(p.y-c.y);
		A[0][1] += (p.x-c.x)*(p.y-c.y);
		A[1][0] += (p.y-c.y)*(p.x-c.x);
	}
	A[0][0] /= m;
	A[1][1] /= m;
	A[0][1] /= m;
	A[1][0] /= m;
}
