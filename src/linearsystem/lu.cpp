/*! 
  @file lu.cpp
	
  @brief LU分解
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// 三角分解
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
int LUDecomp(vector< vector<double> > &A, int n)
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
int LUSolver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n)
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
 * LU分解で逆行列を計算
 * @param[in]  mat 計算したい行列
 * @param[out] inv 逆行列
 * @return 0:成功,1:失敗(n <= 0)
 */
int LUInverse(vector< vector<double> > &mat, vector< vector<double> > &inv)
{
	int n = (int)mat.size();
	if(n <= 0) return 1;

	// 行列を1回だけLU分解
	LUDecomp(mat, n);

	vector<double> b(n), x(n);

	// 行ごとに逆行列を算出
	for(int j = 0; j < n; ++j){
		for(int i = 0; i < n; ++i){
			b[i] = (i == j ? 1.0 : 0.0);
		}

		LUSolver(mat, b, x, n);
		for(int i = 0; i < n; ++i){
			inv[i][j] = x[i];
		}
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

	// LU分解
	LUDecomp(A, n);

	// LU分解結果を画面表示(L行列とU行列をひとつにまとめた行列として表示)
	cout << "LU = " << endl;
	OutputMatrix(A, n, n);
	cout << endl;


	// LU分解を用いて連立1次方程式を解く
	vector<double> x(n);
	LUSolver(A, b, x, n);
	
	// 結果の画面表示
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << x[i] << (i == n-1 ? "" : ", ");
	}
	cout << endl;


	// LU分解を用いて逆行列を算出する
	vector< vector<double> > invA(n, vector<double>(n));
	for(int i = 0; i < n; ++i){
		A[i].assign(Ab[i].begin(), Ab[i].begin()+n);
	}
	LUInverse(A, invA);

	// 逆行列を画面表示
	cout << "inv(A) = " << endl;
	OutputMatrix(invA, n, n);
	cout << endl;

	// 逆行列チェック
	vector< vector<double> > C(n, vector<double>(n));
	MulMatrix(Ab, invA, C, n);

	// 元の行列と計算した逆行列を掛けた結果の行列を画面表示(これがほぼ単位行列ならOK)
	cout << "I = " << endl;
	OutputMatrix(C, n, n);
	cout << endl;


	return 0;
}


