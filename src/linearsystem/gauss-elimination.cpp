/*! 
  @file gauss-elimination.chpp
	
  @brief ガウスの消去法(Gaussian elimination)
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


//-----------------------------------------------------------------------------
// 関数
//-----------------------------------------------------------------------------
/*!
 * ガウスの消去法(ピボット交換なし)
 * @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 * @return 
 */
int GaussElimination(vector< vector<double> > &A, int n)
{
	// 前進消去(forward elimination)
	//  - 対角要素をのぞいた左下要素をすべて0にする(上三角行列にする)
	for(int k = 0; k < n-1; ++k){
		double akk = A[k][k];
		for(int i = k+1; i < n; ++i){
			double aik = A[i][k];
			for(int j = k; j < n+1; ++j){ // 確認のため左下要素が0になるようにj=kとしたが，実際にはj=k+1でよい
				A[i][j] = A[i][j]-aik*(A[k][j]/akk);
			}
		}
		cout << "k = " << k << endl;
		OutputMatrix(A, n, n+1);
		cout << endl;
	}

	// 後退代入(back substitution)
	//  - x_nの解はb_n/a_nn，x_nをさらにn-1行の式に代入することでx_(n-1)を求める．
	//  - この作業を1行目まで続けることですべての解を得る．
	A[n-1][n] = A[n-1][n]/A[n-1][n-1];
	for(int i = n-2; i >= 0; --i){
		double ax = 0.0;
		for(int j = i+1; j < n; ++j){
			ax += A[i][j]*A[j][n];
		}
		A[i][n] = (A[i][n]-ax)/A[i][i];

		cout << "i = " << i << endl;
		OutputMatrix(A, n, n+1);
		cout << endl;
	}

	return 0;
}

/*!
 * ピボット選択(Pivoting)
 *  - 行入れ替えだけの部分的ピボッティング
 * @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列
 * @param[in] n n元連立一次方程式
 * @param[in] k 対象行
 */
inline void Pivoting(vector< vector<double> > &A, int n, int k)
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
	// k != pならば行を交換(ピボット選択)
	if(k != p) swap(A[k], A[p]);
}

/*!
 * ガウスの消去法(行に関するピボット選択(部分ピボッティング)あり)
 * @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 */
int GaussEliminationWithPivoting(vector< vector<double> > &A, int n)
{
	// 前進消去(forward elimination)
	//  - 対角要素をのぞいた左下要素をすべて0にする(上三角行列にする)
	for(int k = 0; k < n-1; ++k){
		// ピボット選択
		Pivoting(A, n, k);

		double akk = A[k][k];
		for(int i = k+1; i < n; ++i){
			double aik = A[i][k];
			for(int j = k; j < n+1; ++j){ // 確認のため左下要素が0になるようにj=kとしたが，実際にはj=k+1でよい
				A[i][j] = A[i][j]-aik*(A[k][j]/akk);
			}
		}
		cout << "k = " << k << endl;
		OutputMatrix(A, n, n+1);
		cout << endl;
	}

	// 後退代入(back substitution)
	//  - x_nの解はb_n/a_nn，x_nをさらにn-1行の式に代入することでx_(n-1)を求める．
	//  - この作業を1行目まで続けることですべての解を得る．
	A[n-1][n] = A[n-1][n]/A[n-1][n-1];
	for(int i = n-2; i >= 0; --i){
		double ax = 0.0;
		for(int j = i+1; j < n; ++j){
			ax += A[i][j]*A[j][n];
		}
		A[i][n] = (A[i][n]-ax)/A[i][i];
	}

	return 0;
}


/*!
 * スケーリング処理
 *  - ガウスの消去法の前処理として使用
 * @param[inout] A n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列
 * @param[in] n n元連立一次方程式
 * @return 1:成功
 */
int ScalingForGauss(vector< vector<double> > &A, int n)
{
	// ガウスの消去法において桁落ちを防ぐために，
	// 各係数を各行の最大値で割ることで正規化する．
	for(int i = 0; i < n; ++i){
		// 行要素で絶対値最大のものを検索
		double am = fabs(A[i][0]);
		for(int j = 1; j < n; ++j){	// 右辺項b_i=a_i,nは含まない
			if(fabs(A[i][j]) > am){
				am = fabs(A[i][j]);
			}
		}

		// すべての行要素を最大値で割る
		for(int j = 0; j < n+1; ++j){
			A[i][j] /= am;
		}
	}

	return 1;
}

//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込む場合
	vector< vector<double> > A;
	ReadMatrix("matrix.txt", ",", A);

	// 2次元配列に初期値として直接値を設定する場合
	//vector< vector<double> > A{ {2, 1, 3, 9},
 //                               {1, 3, 2, 1},
 //                               {3, 4, 3, 4} };

	// 行列のサイズの取得と確認のための画面表示
	int n = (int)A.size();	// n元連立一次方程式
	int m = (int)A[0].size();
	cout << "A(" << n << " x " << m << ") = " << endl;
	OutputMatrix(A, n, n+1);
	cout << endl;

	// ガウスの消去法で線形システムを解く
	//ScalingForGauss(A, n);			// スケーリング処理
	//GaussElimination(A, n);			// ピボッティングなし
	GaussEliminationWithPivoting(A, n);	// ピボッティングあり

	// 結果の表示
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << A[i][n] << (i == n-1 ? "" : ", ");
	}
	cout << endl;

	return 0;
}


