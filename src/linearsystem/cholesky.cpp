/*! 
  @file cholesky.cpp
	
  @brief 修正/不完全コレスキー分解
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// コレスキー分解
//-----------------------------------------------------------------------------
/*!
 * コレスキー分解
 *  - 正定値対称行列A(n×n)を下三角行列(L:Lower triangular matrix)とその転置(L^T)に分解する
 *  - L: i >= j (i>jとi==jに分けて処理)
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int CholeskyDecomp(vector< vector<double> > &A, vector< vector<double> > &L, int n)
{
	if(n <= 0) return 1;

	for(int j = 0; j < n; ++j){
		// i == jについて解く
		double ll = A[j][j];
		for(int k = 0; k < j; ++k){
			ll -= L[j][k]*L[j][k];
		}
		L[j][j] = sqrt(ll);

		for(int i = j+1; i < n; ++i){
			ll = A[i][j];
			for(int k = 0; k < j; ++k){
				ll -= L[i][k]*L[j][k];
			}
			L[i][j] = ll/L[j][j];
		}

	}

	return 0;
}


/*!
 * 修正コレスキー分解(modified Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii = 1とした場合
 *  - L: i > jの要素が非ゼロで対角成分は1
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int ModifiedCholeskyDecomp(const vector< vector<double> > &A, vector< vector<double> > &L, vector<double> &d, int n)
{
	if(n <= 0) return 1;

	d[0] = A[0][0];
	L[0][0] = 1.0;

	for(int i = 1; i < n; ++i){
		// i < k の場合
		for(int j = 0; j < i; ++j){
			double lld = A[i][j];
			for(int k = 0; k < j; ++k){
				lld -= L[i][k]*L[j][k]*d[k];
			}
			L[i][j] = (1.0/d[j])*lld;
		}

		// i == k の場合
		double ld = A[i][i];
		for(int k = 0; k < i; ++k){
			ld -= L[i][k]*L[i][k]*d[k];
		}
		d[i] = ld;
		L[i][i] = 1.0;
	}

	return 0;
}


/*!
 * 不完全コレスキー分解(incomplete Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii = 1とした場合
 *  - L: i > jの要素が非ゼロで対角成分は1
 *  - 行列Aの値が0である要素に対応する部分を飛ばす
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int IncompleteCholeskyDecomp(const vector< vector<double> > &A, vector< vector<double> > &L, vector<double> &d, int n)
{
	if(n <= 0) return 1;

	d[0] = A[0][0];
	L[0][0] = 1.0;

	for(int i = 1; i < n; ++i){
		// i < k の場合
		for(int j = 0; j < i; ++j){
			if(fabs(A[i][j]) < 1.0e-10) continue;

			double lld = A[i][j];
			for(int k = 0; k < j; ++k){
				lld -= L[i][k]*L[j][k]*d[k];
			}
			L[i][j] = (1.0/d[j])*lld;
		}

		// i == k の場合
		double ld = A[i][i];
		for(int k = 0; k < i; ++k){
			ld -= L[i][k]*L[i][k]*d[k];
		}
		d[i] = ld;
		L[i][i] = 1.0;
	}

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

	// 左辺の行列の抽出
	A.resize(n, vector<double>(n));
	for(int i = 0; i < n; ++i){
		A[i].assign(Ab[i].begin(), Ab[i].begin()+n);
	}

	// 読み込んだ行列を確認用に画面表示
	cout << "A(" << n << " x " << n+1 << ") = " << endl;
	OutputMatrix(A, n, n);
	cout << endl;

	// コレスキー分解
	vector< vector<double> > L(n, vector<double>(n, 0.0));
	vector<double> d;
	//CholeskyDecomp(A, L, n);

	// 修正コレスキー分解，不完全コレスキー分解
	d.resize(n, 0.0);
	ModifiedCholeskyDecomp(A, L, d, n);
	//IncompleteCholeskyDecomp(A, L, d, n);

	cout << "L = " << endl;
	OutputMatrix(L, n, n);

	if(d.empty()){
		// 分解結果のチェック用
		vector< vector<double> > LL(n, vector<double>(n, 0.0));
		MulMatrix(L, transpose(L, n), LL, n);

		// 分解した結果を掛け合わせた行列を画面表示(これが元の行列Aと同じならOK)
		cout << "LL^T = " << endl;
		OutputMatrix(LL, n, n);
		cout << endl;
	}
	else{
		vector< vector<double> > D(n, vector<double>(n, 0.0));
		for(int i = 0; i < n; ++i) D[i][i] = d[i];
		cout << "D = " << endl;
		OutputMatrix(D, n, n);
		cout << endl;

		// 分解結果のチェック用
		vector< vector<double> > LD(n, vector<double>(n, 0.0));
		vector< vector<double> > LDL(n, vector<double>(n, 0.0));
		MulMatrix(L, D, LD, n);
		MulMatrix(LD, transpose(L, n), LDL, n);

		// 分解した結果を掛け合わせた行列を画面表示(これが元の行列Aと同じならOK)
		cout << "LDL^T = " << endl;
		OutputMatrix(LDL, n, n);
		cout << endl;
	}

	return 0;
}


