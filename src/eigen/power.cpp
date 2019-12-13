/*! 
  @file power.cpp
	
  @brief 固有値・固有ベクトル
         べき乗法(power method)
 
  @author Makoto Fujisawa
  @date 2019-12
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"




/*!
 * べき乗法による固有値の算出
 * @param[in] A 正方行列(n×n)
 * @param[out] lambda 絶対値最大固有値
 * @param[inout] v 絶対値最大固有値に対応する固有ベクトル(大きさnの配列) - vの初期値を入れておく
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 反復回数
 */
int eigen_power(vector< vector<double> > &A, double &lambda, vector<double> &x, int n, int &max_iter, double &eps)
{
	vector<double> x1(n); // v^(k+1)格納用
	double l2 = dot(x, x, n), e;
    int k;
    for(k = 0; k < max_iter; ++k){
		// |x^(k)|=1となるように正規化
		double l = sqrt(l2); // |x^(k)|の計算
		x = mul_sv(1.0/l, x, n); // |x^(k)|で割る(1/|x|を掛ける)

		// x^(k+1) = A x^(k)の計算
		x1 = mul_mv(A, x, n);

		// 固有値λの計算(|x^(k)|=1で正規化されていることが前提)
		lambda = dot(x, x1, n);

		cout << k << " : lambda = " << lambda << ",  v = " << x << endl;

		// 収束判定
		l2 = dot(x1, x1, n);
		if((e = fabs(l2-lambda*lambda)) < eps*eps) break;

		x = x1;
    }
	x = x1;
	x = mul_sv(1.0/sqrt(l2), x, n); // 最後に回を正規化しておく

	max_iter = k;
	eps = sqrt(e);

    return 1;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込んで解く
	vector< vector<double> > A;

	ReadMatrix("matrix.txt", ",", A);

	int n = (int)A.size();	// n×n行列
	cout << "A(" << n << " x " << n << ") = " << endl;
	OutputMatrix(A, n, n);
	cout << endl;

	double lambda;
	vector<double> v(n, 1.0); // vを1で初期化
	cout.precision(8);

	int max_iter = 100;
	double eps = 1e-6;
	eigen_power(A, lambda, v, n, max_iter, eps);

	// 絶対値最大固有値と固有ベクトルの表示
	cout << "lambda_max = " << lambda << endl;
	cout << "v = " << v << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


