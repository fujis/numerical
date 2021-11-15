/*! 
  @file jacobi.cpp
	
  @brief ヤコビ(Jacobi)反復法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// 反復法による線形システムソルバ
//-----------------------------------------------------------------------------
/*!
 * ヤコビ反復法(Jacobi iterative method)
 *  - 解が収束するのは
 *      ・対角有利(diagonal dominant, 対角要素の絶対値>その行の他の要素の絶対値の和)
 *      ・係数行列が対称(symmetric)かつ正定(positive definite)
 *    のどちらかの場合
 * @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の拡大行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 */
int JacobiIteration(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
	vector<double> x(n, 0.0);	// 初期値はすべて0とする
	vector<double> y(n, 0.0);

	double e = 0.0; // 誤差平均
	int k;
	for(k = 0; k < max_iter; ++k){
		e = 0.0;
		// 現在の値を代入して，次の解候補を計算
		for(int i = 0; i < n; ++i){
			y[i] = A[i][n];
			for(int j = 0; j < n; ++j){
				y[i] -= (j != i ? A[i][j]*x[j] : 0.0);
			}
			y[i] /= A[i][i];

			e += fabs(y[i]-x[i]);
		}

		// 確認のため現在の解を画面出力
		cout.precision(8);
		cout << k << " : ";
		for(int i = 0; i < n; ++i) cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");
		//for(int i = 0; i < n; ++i) cout << ", " << fabs(y[i]-x[i]);
		//cout << ", " << e/n;
		cout << endl;


		// 収束判定
		int l;
		for(l = 0; l < n; ++l){
			if(fabs(y[l]-x[l]) > eps){	// 絶対誤差の場合
			//if(fabs((y[l]-x[l])/y[l])){	// 相対誤差の場合
				break;
			}
		}
		if(l == n){ // すべての解が許容誤差以下なら反復終了
			break;
		}

		swap(x, y);
	}

	max_iter = k;
	eps = e;

	// 解をA_(i,n)に格納
	for(int i = 0; i < n; ++i){
		A[i][n] = y[i];
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

	// 読み込んだ行列を確認用に画面表示
	cout << "A(" << n << " x " << m << ") = " << endl;
	OutputMatrix(A, n, n+1);
	cout << endl;

	// ヤコビ反復法で線形システムを解く
	int max_iter = 100;
	double eps = 1e-6;
	JacobiIteration(A, n, max_iter, eps);

	// 結果の画面表示
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << A[i][n] << (i == n-1 ? "" : ", ");
	}
	cout << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


