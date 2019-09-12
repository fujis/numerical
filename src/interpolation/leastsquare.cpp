/*! 
  @file spline.cpp
	
  @brief スプライン補間
 
  @author Makoto Fujisawa
  @date 2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// LU分解(最小自乗法での係数計算に用いる)
//-----------------------------------------------------------------------------
/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
 *  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
 *  - LとUを一つの行列にまとめた形で結果を返す
 * @param[inout] A n×nの係数行列．LU分解した結果を格納する．
 * @param[in] n 行列の大きさ
 * @return 1:成功,0:失敗
 */
int luDecomp(vector<double> &A, int n)
{
	if(n <= 0) return 0;

	for(int i = 0; i < n; ++i){
		// l_ijの計算(i >= j)
		for(int j = 0; j <= i; ++j){
			double lu = A[i*n+j];
			for(int k = 0; k < j; ++k){
				lu -= A[i*n+k]*A[k*n+j];    // l_ik * u_kj
			}
			A[i*n+j] = lu;
		}

		// u_ijの計算(i < j)
		for(int j = i+1; j < n; ++j){
			double lu = A[i*n+j];
			for(int k = 0; k < i; ++k){
				lu -= A[i*n+k]*A[k*n+j];    // l_ik * u_kj
			}
			A[i*n+j] = lu/A[i*n+i];
		}
	}

	return 1;
}

/*!
 * LU分解した行列a(n×n)から前進代入・後退代入によりA・x=bを解く
 * @param[in] A LU分解された行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @return 1:成功,0:失敗
 */
int luSolver(const vector<double> &A, const vector<double> &b, vector<double> &x, int n)
{
	if(n <= 0) return 0;

	// 前進代入(forward substitution)
	//  LY=bからYを計算
	for(int i = 0; i < n; ++i){
		double bly = b[i];
		for(int j = 0; j < i; ++j){
			bly -= A[i*n+j]*x[j];
		}
		x[i] = bly/A[i*n+i];
	}

	// 後退代入(back substitution)
	//  UX=YからXを計算
	for(int i = n-1; i >= 0; --i){
		double yux = x[i];
		for(int j = i+1; j < n; ++j){
			yux -= A[i*n+j]*x[j];
		}
		x[i] = yux;
	}

	return 1;
}

//! 階乗計算
inline int factorial(int n)
{
	int f = 1;
	for(int i = 1; i <= n; ++i) f *= i;
	return f;
}

//-----------------------------------------------------------------------------
// 補間法
//-----------------------------------------------------------------------------
/*!
 * 最小二乗法による補間
 *  - 
 * @param[in] yi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
 * @param[in] n データ数
 * @param[in] m フィッティングする多項式の次数
 * @param[in] x 補間した値が必要な位置x
 * @param[out] ans 解
 * @return
 */
int leastsquare_interpolation(const vector<double> &yi, const vector<double> &xi, int n, int m, double x, double &ans)
{
	int dim = 1;	// 空間次元数(今回は1次元データなので1)
	int dim_basis = factorial(dim+m)/(factorial(m)*factorial(dim)); // 多項式の次数から行列のサイズ(=基底ベクトルの次元数)を計算

	vector<double> A, b, c; // Ac=b -> c=A^-1 b
	A.resize(dim_basis*dim_basis, 0.0);		// 左辺係数行列
	b.resize(dim_basis, 0.0);				// 右辺項
	c.resize(dim_basis, 0.0);				// 計算結果の多項式係数

	// 係数行列Aと右辺項bの計算
	vector<double> basis(dim_basis, 0.0);
	for(int k = 0; k < n; ++k){
		// 多項式の基底ベクトルの計算(多項式以外でフィッティングする場合はこの部分を変えれば良い)
		double xb = 1;
		for(int i = 0; i < dim_basis; ++i){
			basis[i] = xb; xb *= xi[k];
		}

		// 係数行列Aと右辺項bの計算
		for(int i = 0; i < dim_basis; ++i){
			for(int j = 0; j < dim_basis; ++j){
				A[i*dim_basis+j] += basis[i]*basis[j];
			}
			b[i] += basis[i]*yi[k];
		}
	}

	// Ax=bをLU分解で解く(疎行列とは限らないのでCGソルバは使わない)
	luDecomp(A, dim_basis);
	luSolver(A, b, c, dim_basis);

	// 基底ベクトルに求めた係数を掛けて行くことで位置xにおける値yを計算
	ans = 0;
	double xb = 1;
	for(int i = 0; i < dim_basis; ++i){
		ans += c[i]*xb; xb *= x;
	}

	return 0;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncExp;
	//double(*func)(double) = FuncLinear;

	double x = 0.5;
	double fx;
	double gt = func(x); // 真値
	cout.precision(10);

	// 補間用の点
	vector<double> xi, yi;
	xi.push_back(0.0);
	yi.push_back(func(xi.back()));
	xi.push_back(0.2);
	yi.push_back(func(xi.back()));
	xi.push_back(0.4);
	yi.push_back(func(xi.back()));
	xi.push_back(0.6);
	yi.push_back(func(xi.back()));
	xi.push_back(0.8);
	yi.push_back(func(xi.back()));
	xi.push_back(1.0);
	yi.push_back(func(xi.back()));

	// n次多項式による補間
	leastsquare_interpolation(yi, xi, xi.size(), 2, x, fx);
	cout << "f_ls2(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	leastsquare_interpolation(yi, xi, xi.size(), 3, x, fx);
	cout << "f_ls3(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	leastsquare_interpolation(yi, xi, xi.size(), 4, x, fx);
	cout << "f_ls4(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	leastsquare_interpolation(yi, xi, xi.size(), 5, x, fx);
	cout << "f_ls5(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	cout << "ground truth = " << gt << endl;

		
	return 0;
}


