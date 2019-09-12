/*! 
  @file gauss-jordan.cpp
	
  @brief ガウス-ジョルダン(Gauss-Jordan)の消去法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// ガウス・ジョルダン法による逆行列計算
//  - ガウスの消去法を拡張して逆行列計算に用いる
//-----------------------------------------------------------------------------
/*!
 * ガウス・ジョルダン法(ピボット選択なし)
 * @param[inout] A n×2nの拡張行列
 * @param[in] n n元連立一次方程式
 */
int GaussJordan(vector< vector<double> > &A, int n)
{
	// 拡張行列の右半分を単位行列にする
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			A[i][j+n] = (i == j ? 1.0 : 0.0);
		}
	}

	// ガウス・ジョルダン法(Gauss-Jordan method)で逆行列計算
	for(int k = 0; k < n; ++k){
		double akk = A[k][k];
		// 対角要素を1にするために，k行目のすべての要素をa_kkで割る
		for(int j = 0; j < 2*n; ++j){
			A[k][j] /= akk;
		}

		// k列目の非対角要素を0にする
		for(int i = 0; i < n; ++i){
			if(i == k) continue;
			double aik = A[i][k];
			for(int j = 0; j < 2*n; ++j){
				A[i][j] -= A[k][j]*aik;
			}
		}
	}

	return 0;
}

/*!
 * ピボット選択(Pivoting)
 *  - 行入れ替えだけの部分的ピボッティング
 * @param[inout] A n×2nの拡張行列
 * @param[in] n n元連立一次方程式
 * @param[in] k 対象行
 * @return 交換した行
 */
inline int Pivoting(vector< vector<double> > &A, int n, int k)
{
	// k行目以降でk列目の絶対値が最も大きい要素を持つ行を検索
	int p = k; // 絶対値が最大の行
	double am = fabs(A[k][k]);// 最大値
	for(int i = k+1; i < n; ++i){
		if(fabs(A[i][k]) > am){
			p = i;
			am = fabs(A[i][k]);
		}
	}
	// k != pならば行を交換
	if(k != p) swap(A[k], A[p]);

	return p;
}


/*!
 * ガウス・ジョルダン法(ピボット選択あり)
 * @param[inout] A n×2nの拡張行列
 * @param[in] n n元連立一次方程式
 */
int GaussJordanWithPivoting(vector< vector<double> > &A, int n)
{
	// 拡張行列の右半分を単位行列にする
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			A[i][j+n] = (i == j ? 1.0 : 0.0);
		}
	}

	// ガウス・ジョルダン法(Gauss-Jordan method)で逆行列計算
	vector<int> ptable(n);
	for(int k = 0; k < n; ++k) ptable[k] = k;
	for(int k = 0; k < n; ++k){
		// ピボット選択
		int p = Pivoting(A, n, k);
		if(p != k){
			swap(ptable[k], ptable[p]);
		}

		double akk = A[k][k];
		// 対角要素を1にするために，k行目のすべての要素をa_kkで割る
		for(int j = 0; j < 2*n; ++j){
			A[k][j] /= akk;
		}

		// k列目の非対角要素を0にする
		for(int i = 0; i < n; ++i){
			if(i == k) continue;
			double aik = A[i][k];
			for(int j = 0; j < 2*n; ++j){
				A[i][j] -= A[k][j]*aik;
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込む
	vector< vector<double> > A, A0;
	ReadMatrix("matrix.txt", ",", A);
	A0 = A;	// チェック用に元の行列を確保

	int n = (int)A.size();

	// 読み込んだ行列を確認用に画面表示
	cout << "A(" << n << " x " << n << ") = " << endl;
	OutputMatrix(A, n, n);
	cout << endl;

	// 拡張行列を作成
	for(int i = 0; i < n; ++i){
		A[i].resize(2*n);
		for(int j = 0; j < n; ++j){
			A[i][j+n] = (i == j ? 1.0 : 0.0);
		}
	}

	// ガウスジョルダンで逆行列を求める
	//GaussJordan(A, n);				// ピボッティングなし
	GaussJordanWithPivoting(A, n);	// ピボッティングあり

	// 拡張行列A(nx2n)全体の画面表示
	OutputMatrix(A, n, 2*n);
	cout << endl;

	// 拡張行列から逆行列部分だけを抽出
	vector< vector<double> > invA(n, vector<double>(n));
	for(int i = 0; i < n; ++i){
		invA[i].assign(A[i].begin()+n, A[i].end());
	}

	// 逆行列部分のみの画面表示
	cout << "A^-1 = " << endl;
	OutputMatrix(invA, n, n);
	cout << endl;

	// 逆行列チェック
	vector< vector<double> > C(n, vector<double>(n));
	MulMatrix(A0, invA, C, n);

	// 元の行列と計算した逆行列を掛けた結果の行列を画面表示(これがほぼ単位行列ならOK)
	cout << "A A^-1 = " << endl;
	OutputMatrix(C, n, n);
	cout << endl;

	return 0;
}


