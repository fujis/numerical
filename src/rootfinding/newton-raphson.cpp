/*! 
  @file newton-raphson.cpp
	
  @brief ニュートン・ラフソン法(Newton-Raphson method)
         もしくは単にニュートン法
 
  @author Makoto Fujisawa
  @date 2012-06,2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// ニュートン法による求根問題の解法
//-----------------------------------------------------------------------------
/*!
 * ニュートン・ラフソン法(Newton-Raphson method)
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] dfunc 導関数値を与える関数ポインタ
 * @param[inout] x 解(探索開始位置を渡す)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 */
int newton(double func(const double), double dfunc(const double), 
		   double &x, int &max_iter, double &eps)
{
	double f, df, dx;

	int k;
	for(k = 0; k < max_iter; ++k){
		// 現在の位置xにおける関数値と導関数の計算
		f = func(x);
		df = dfunc(x);

		// 確認用の画面出力
		cout << k << " : f(" << x << ") = " << f << endl;

		// 導関数の結果から次の位置を計算
		x = x-f/df;

		// 収束判定
		dx = fabs(f/df);
		if(dx < eps || fabs(f) < eps){
			break;
		}
	}

	max_iter = k; eps = dx;
	
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

struct FUNCTION
{
	double (*func)(const vector<double>&);
	vector<double> (*dfunc)(const vector<double>&);
};

/*!
 * ニュートン・ラフソン法(n次元版)
 * @param[in] funcs 関数群
 * @param[inout] x 解(探索開始位置を渡す)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 */
int newton(vector<FUNCTION> funcs, vector<double> &x, int n, int &max_iter, double &eps)
{
	vector<double> f(n, 0.0);
	vector< vector<double> > J(n);
	for(int i = 0; i < n; ++i) J[i].resize(n+1, 0.0); // 最後の列に右辺項ベクトルを入れるためにn+1にしている
	vector<double> df(n);

	double d = 0.0;
	int k;
	for(k = 0; k < max_iter; ++k){
		// 現在の位置xにおける関数値の計算
		for(int i = 0; i < n; ++i) f[i] = -funcs[i].func(x);

		// ヤコビ行列の計算
		for(int i = 0; i < n; ++i){
			df = funcs[i].dfunc(x);
			for(int j = 0; j < n; ++j) J[i][j] = df[j];
		}

		// 確認用の画面出力
		cout << k << " : ";
		for(int j = 0; j < n; ++j){
			cout << "f" << j << "(";
			for(int i = 0; i < n; ++i) cout << x[i] << (i == n-1 ? ") = " : ",");
			cout << f[j] << (j == n-1 ? "" : ", ");
		}
		cout << endl;
		
		

		// 線形システムを解いてδを計算
		for(int i = 0; i < n; ++i) J[i][n] = f[i];
		GaussEliminationWithPivoting(J, n);

		// xを更新
		for(int i = 0; i < n; ++i) x[i] += J[i][n];

		// 収束判定
		d = 0.0;
		for(int i = 0; i < n; ++i) d += fabs(f[i]);
		if(d < eps){
			break;
		}
	}

	max_iter = k; eps = d;

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	//// 探索開始位置
	//double x = 1;

	//// ニュートン法でf(x)=0を解く
	//int max_iter = 100;
	//double eps = 1e-6;
	//newton(FuncPi, DFuncPi, x, max_iter, eps);

	//// 結果の画面表示
	//cout.precision(12);
	//cout << "x = " << x << endl;
	//cout << "iter = " << max_iter << ", eps = " << eps << endl;
	//cout << endl;

	// 多次元のニュートン法
	vector<double> xv(2);
	xv[0] = 0.0; xv[1] = -1.0;
	FUNCTION f;
	vector<FUNCTION> funcs;
	f.func = Func4;	f.dfunc = DFunc4;
	funcs.push_back(f);
	f.func = Func4a;	f.dfunc = DFunc4a;
	funcs.push_back(f);

	int max_iter = 100;
	double eps = 1e-6;
	newton(funcs, xv, 2, max_iter, eps);

	// 結果の画面表示
	cout << "(x1,x2) = (" << xv[0] << ", " << xv[1] << ")" << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}

