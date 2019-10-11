/*! 
  @file gauss-seidel.cpp
	
  @brief ガウス-ザイデル(Gauss Seidel)反復法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// ガウス・ザイデル法による線形システムソルバ
//-----------------------------------------------------------------------------
/*!
 * ガウス-ザイデル反復法(Gauss Seidel iterative method)
 *  - 解が収束するのは
 *      ・対角有利(diagonal dominant, 対角要素の絶対値>その行の他の要素の絶対値の和)
 *      ・係数行列が対称(symmetric)かつ正定(positive definite)
 *      ・Σ_j |a_ij/a_ii| < 1 (i = 1～n, j != i) 
 * @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 */
int GaussSeidel(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
	vector<double> x(n, 0.0);	// 初期値はすべて0とする
	cout.precision(8);

	double e = 0.0;	// 誤差
	int k;	// 計算反復回数
	for(k = 0; k < max_iter; ++k){
		//cout << k;

		// 現在の値を代入して，次の解候補を計算
		int l = 0;
		e = 0;
		for(int i = 0; i < n; ++i){
			double tmp = x[i];
			x[i] = A[i][n];
			for(int j = 0; j < n; ++j){
				x[i] -= (j != i ? A[i][j]*x[j] : 0.0);
			}
			x[i] /= A[i][i];

			if(fabs(tmp-x[i]) > eps){	// 絶対誤差の場合
			//if(fabs((tmp-x[i])/tmp)){	// 相対誤差の場合
				e += fabs(tmp-x[i]);
				l++;
			}
			// 確認のため現在の誤差を画面出力
			//cout << ", " << fabs(tmp-x[i]);
		}
		//cout << ", " << e/n;

		//// 確認のため現在の解を画面出力
		cout << k << " : ";
		for(int i = 0; i < n; ++i) cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");

		cout << endl;

		// 収束判定
		if(l == 0){ // すべての解が許容誤差以下なら反復終了
			break;
		}
	}

	max_iter = k;
	eps = e/n; // 平均誤差

	for(int i = 0; i < n; ++i){
		A[i][n] = x[i];
	}

	return 0;
}





//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込んで解く
	vector< vector<double> > A;
	ReadMatrix("matrix.txt", ",", A);

	int n = (int)A.size();	// n元連立一次方程式
	int m = (int)A[0].size();

	// 読み込んだ行列を確認用に画面表示 (n+1列目は右辺項b)
	cout << "A(" << n << " x " << m << ") = " << endl;
	OutputMatrix(A, n, n+1);
	cout << endl;

	// ガウス・ザイデル法で線形システムを解く
	int max_iter = 100;
	double eps = 1e-6;
	GaussSeidel(A, n, max_iter, eps);

	// 結果の画面表示
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << A[i][n] << (i == n-1 ? "" : ", ");
	}
	cout << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


