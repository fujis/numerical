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

//-----------------------------------------------------------------------------
// 補間法
//-----------------------------------------------------------------------------
/*!
 * 最小二乗法による補間
 *  - 1次元データ(xi,fi)に対して，1次元n次多項式をフィッティング
 * @param[in] fi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列(位置は昇順でソートされている必要がある)
 * @param[in] m データ数
 * @param[in] n フィッティングする多項式の次数
 * @param[in] x 補間した値が必要な位置x
 * @param[out] c 係数ベクトル
 * @return 補間値
 */
double leastsquare_interpolation(const vector<double> &fi, const vector<double> &xi, int m, int n, double x, vector<double> &c)
{
	// 多項式の次数から行列のサイズ(=基底ベクトルの次元数)を計算
	int dim_b = n+1; // 1次元の場合

	// Ac=y, bは基底ベクトル
	c.resize(dim_b, 0.0);	// 結果の係数ベクトル
	vector<double> b(dim_b, 0.0), y(dim_b, 0.0); // 基底ベクトルと右辺項
	vector< vector<double> > A(dim_b, b);

	// 係数行列Aと右辺項yの計算
	for(int k = 0; k < m; ++k){ // データ分だけ反復
		// 多項式の基底ベクトルの計算(多項式以外でフィッティングする場合はこの部分を変えれば良い)
		double xb = 1;
		for(int i = 0; i < dim_b; ++i){
			b[i] = xb; xb *= xi[k];
		}

		// 係数行列Aと右辺項bの計算
		for(int i = 0; i < dim_b; ++i){
			for(int j = 0; j < dim_b; ++j){
				A[i][j] += b[i]*b[j];
			}
			y[i] += b[i]*fi[k];
		}
	}

	// Ac=yをLU分解で解く(疎行列とは限らないのでCGソルバは使わない)
	LUDecomp(A, dim_b);
	LUSolver(A, y, c, dim_b);

	// 基底ベクトルに求めた係数を掛けて行くことで位置xにおける値yを計算
	double fx = 0;
	double xb = 1;
	for(int i = 0; i < dim_b; ++i){
		fx += c[i]*xb; xb *= x;
	}
	return fx;
}

/*!
 * 1次元n次多項式を計算
 * @param[in] x 値が必要な位置
 * @param[in] c 多項式の係数
 * @param[in] n 多項式の次数
 * @return 多項式を計算した結果
 */
double polynominal1d(double x, const vector<double> &c, int n)
{
	double fx = 0;
	double xb = 1;
	for(int i = 0; i <= n; ++i){
		fx += c[i]*xb; xb *= x;
	}
	return fx;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncExp;
	double x0 = 0, x1 = 1;

	//double(*func)(double) = FuncRunge;
	//double x0 = -1, x1 = 1;

	double x = 0.5, fx;
	double gt = func(x); // 真値
	cout.precision(6);

	// 補間用の点(サンプリング点の数=6)
	vector<double> xi, yi;
	MakeSamplingPoints(x0, x1, (x1-x0)/5, func, xi, yi);
	OutputSamplingPoints(xi, yi); // サンプリング点の画面出力

	// 3次多項式による補間
	vector<double> c; // 係数ベクトル
	fx = leastsquare_interpolation(yi, xi, xi.size(), 3, x, c);
	cout << "f_ls3(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	cout << "ground truth = " << gt << endl;



	//// グラフ描画用にデータ出力
	//int d = 3;  // 多項式の次数
	//int m = 11; // サンプリング点数
	//MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func, xi, yi);
	//leastsquare_interpolation(yi, xi, xi.size(), d, x, c);
	//OutputSamplingPoints(xi, yi, "dat/leastsquare"+TOSTR(m)+"_data.txt"); // サンプリング点のファイル出力
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, c.size()), "dat/leastsquare"+TOSTR(m)+"_d"+TOSTR(d)+".txt");
	//// チェビシェフ節点を用いる場合
	//MakeChebyshevNodes(x0, x1, (x1-x0)/(m-1.0), func, xi, yi);
	//leastsquare_interpolation(yi, xi, xi.size(), d, x, c);
	//OutputSamplingPoints(xi, yi, "dat/leastsquare"+TOSTR(m)+"c_data.txt");
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, c.size()), "dat/leastsquare"+TOSTR(m)+"_d"+TOSTR(d)+"c.txt");

	//// 真値のグラフ作成用ファイル出力
	//OutputFunction(x0, x1, (x1-x0)/200, func, "dat/leastsquare_ground_truth.txt");


	//// ノイズ付きデータに対する最小２乗法
	//MakeSamplingPointsWithWhiteNoise(x0, x1, (x1-x0)/(m-1.0), func, 0.15, xi, yi);
	//leastsquare_interpolation(yi, xi, xi.size(), d, x, c);
	//OutputSamplingPoints(xi, yi, "dat/leastsquare"+TOSTR(m)+"n_data.txt"); // サンプリング点のファイル出力
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, c.size()), "dat/leastsquare"+TOSTR(m)+"n_d"+TOSTR(d)+".txt");

	//// 外れ値付きデータに対する最小２乗法
	//MakeSamplingPointsWithWhiteNoise(x0, x1, (x1-x0)/(m-1.0), func, 0.01, xi, yi);
	//yi[2] = func(xi[2])+1.2; // 外れ値を意図的に追加
	//leastsquare_interpolation(yi, xi, xi.size(), d, x, c);
	//OutputSamplingPoints(xi, yi, "dat/leastsquare"+TOSTR(m)+"o_data.txt"); // サンプリング点のファイル出力
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, c.size()), "dat/leastsquare"+TOSTR(m)+"o_d"+TOSTR(d)+".txt");

	return 0;
}


