/*! 
  @file dka.cpp
	
  @brief DKA法(ワイヤストラス法+アバースの初期値)
		 代数方程式(多項式からなる方程式)の求根問題
 
  @author Makoto Fujisawa
  @date 2012-06
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

#include <complex>	// 複素数
typedef complex<double> complexf;	// double型の複素数を定義


//-----------------------------------------------------------------------------
// 代数方程式
//-----------------------------------------------------------------------------
/*!
 * 係数から代数方程式の値を計算
 * @param[in] x 変数
 * @param[in] c 係数
 * @param[in] n 方程式の次数
 * @return 代数方程式の値
 */
template<class T>
inline T func(T x, const vector<T> &b, int n)
{
	T f = b[n];
	for(int i = 0; i < n; ++i){
		f += b[i]*pow(x, (T)(n-i));
	}
	return f;
}

/*!
 * ホーナー法で代数方程式の値を計算
 * @param[in] x 変数
 * @param[in] b 係数
 * @param[in] n 方程式の次数
 * @return 代数方程式の値
 */
template<class T>
inline T func_h(double x, const vector<T> &b, int n)
{
	T f = b[0];
	for(int i = 1; i <= n; ++i){
		f = b[i]+f*x;
	}
	return f;
}
/*!
 * ホーナー法で代数方程式の導関数値を計算
 * @param[in] x 変数
 * @param[in] b 係数
 * @param[in] n 方程式の次数
 * @return 代数方程式の導関数値
 */
template<class T>
inline T dfunc_h(double x, const vector<T> &b, int n)
{
	T df = n*b[0];
	for(int i = 1; i <= n-1; ++i){
		df = (n-i)*b[i]+df*x;
	}
	return df;
}



//-----------------------------------------------------------------------------
// DKA法による求根問題の解法
//-----------------------------------------------------------------------------
/*!
 * ニュートン・ラフソン法(Newton-Raphson method)
 *  - Aberthの初期値を求めるために用いる
 * @param[in] b 多項式の係数
 * @param[in] n 方程式の次数
 * @param[inout] x 探索開始位置を受け取り，解を返す
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */
int newton(const vector<double> &b, int n, double &x, int &max_iter, double &eps)
{
	double f, df, dx;

	int k;
	for(k = 0; k < max_iter; ++k){
		// 現在の位置xにおける関数値と導関数の計算
		f = func_h(x, b, n);
		df = dfunc_h(x, b, n);

		// 導関数の結果から次の位置を計算
		x = x-f/df;

		// 収束判定

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
 * ホーナー法(組立除法)の係数の値を求める
 * @param[in] a 代数方程式の係数
 * @param[out] b x=x0のときの組立除法の係数
 * @param[in] x0 係数を計算するxの値
 */
template<class T>
inline void horner(const vector<T> &a, vector<T> &b, int n, T x0)
{
	if(n <= 2) return;

	b = a;
	for(int i = 0; i <= n; ++i){
		for(int j = 1; j <= n-i; ++j){
			b[j] += x0*b[j-1];
		}
	}
}

/*!
 * Aberthの方法で初期値を算出
 * @param[out] z 解探索のための初期値
 * @param[in] c 多項式の係数
 * @param[in] n 方程式の次数
 * @param[in] max_iter 半径計算のための最大反復数
 * @param[in] eps 半径計算のための許容誤差
 * @return 1:成功,0:失敗
 */
int aberth(vector<complexf> &z, const vector<complexf> &c, int n, int max_iter, double eps)
{
	// 半径算出のための方程式の係数
	vector<complexf> a;
	complexf zc = -c[1]/(c[0]*(double)n);
	horner(c, a, n, zc);

	vector<double> b(a.size());
	b[0] = abs(a[0]);
	for(int i = 1; i <= n; ++i){
		b[i] = -abs(a[i]);
	}

	// Aberthの初期値の半径をニュートン法で算出
	double r = 100.0;
	newton(b, n, r, max_iter, eps);
	//cout << "r = " << r << endl;

	// Aberthの初期値
	for(int j = 0; j < n; ++j){
		double theta = (2*RX_PI/(double)n)*j+RX_PI/(2.0*n);
		z[j] = zc+r*complexf(cos(theta), sin(theta));
	}

	return 1;
}

/*!
 * ワイヤストラス法(DK公式)
 * @param[inout] z 初期値位置を受け取り，解を返す
 * @param[in] c 多項式の係数
 * @param[in] n 方程式の次数
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */
int weierstrass(vector<complexf> &z, const vector<complexf> &c, int n, int &max_iter, double &eps)
{
	double e = 0.0, ej;

	vector<complexf> zp;
	complexf f, df;
	int k;
	for(k = 0; k < max_iter; ++k){
		zp = z;

		// DK式の計算
		for(int j = 0; j < n; ++j){
			f = func(z[j], c, n);
			df = c[0];
			for(int i = 0; i < n; ++i){
				if(i != j){
					df *= zp[j]-zp[i];
				}
			}

			z[j] = zp[j]-f/df;
		}

		// 誤差の算出
		e = 0.0;
		for(int j = 0; j < n; ++j){
			if((ej = abs(func(z[j], c, n))) > e){
				e = ej;
			}
		}

		// 収束判定
		if(e < eps){
			max_iter = k; eps = e;
			return 1;
		}
	}

	eps = e;
	return 0;
}




//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	vector<complexf> c;
	vector<double> c0;
	vector<complexf>::iterator itr;

	// 代数方程式の係数を読み込む
	if(ReadAlgebra("algebra.txt", ",", c0)){
		for(vector<double>::iterator i = c0.begin(); i != c0.end(); ++i){
			c.push_back(complex<double>(*i));
		}
	}

	int n = (int)c.size()-1;		// 代数方程式の次数
	int max_iter = 100;
	double eps = 1e-8;

	// Aberthの初期値
	vector<complexf> z(n);
	aberth(z, c, n, max_iter, eps);

	// ワイヤストラス法(DK公式)で解を求める
	weierstrass(z, c, n, max_iter, eps);

	// 結果の画面表示
	cout << "solutions : " << endl;
	for(itr = z.begin(); itr != z.end(); ++itr){
		cout << "  " << *itr << endl;
	}
	cout << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;



	return 0;
}


