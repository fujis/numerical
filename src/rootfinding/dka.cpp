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


//! complex型に<<オペレータを設定(演算子のオーバーロード)
template<typename T>
inline std::ostream &operator<<(std::ostream &out, const complex<T> &x)
{
	out << x.real();
	if(fabs(x.imag()) > 1e-6) out << " + " << x.imag() << "i";
	return out;
}


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
inline double func_h(double x, const vector<double> &b, int n)
{
	double f = b[0];
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
inline double dfunc_h(double x, const vector<double> &b, int n)
{
	double df = n*b[0];
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
  * ホーナー法(組立除法)
   - P(x) = a0 x^n + a1 x^(n-1) + ... + a_(n-1) x + a_n を (x-b) で割ったときの商と余りを返す
   - 商はn-1次の多項式の係数として返す
  * @param[in] a 代数方程式の係数
  * @param[in] b 割る1次式の係数(x-b)
  * @param[out] c 商であるn-1次の多項式の係数
  * @param[out] rm あまり
  * @param[in] n 方程式の次数(配列aの大きさはn+1,配列cはn)
  */
inline void horner(const vector<double> &a, double b, vector<double> &c, double &rm, int n)
{
	if(n <= 1) return;
	rm = a[0]; // 最終的に余りになる
	c.resize(n);
	for(int i = 1; i < n+1; ++i){
		c[i-1] = rm;
		rm *= b;
		rm += a[i];
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
	vector<complexf> cd(n+1, 1.0); // 係数c'
	double c1n = -c[1].real()/n; 
	cd[0] = c[0];
	double rm; // 組立除法での余り
	vector<double> a(n+1), tmp(n); // 係数格納用の一時的な変数
	for(int i = 0; i <= n; ++i) a[i] = c[i].real();
	// zの多項式をz+c1/nで割っていくことでwの多項式の係数を求める
	for(int i = n; i > 1; --i){
		horner(a, c1n, tmp, rm, i);
		cd[i] = rm;
		a = tmp;
	}
	cd[1] = a[1]+c1n;

	// 多項式S(w)の係数
	vector<double> b(cd.size());
	b[0] = abs(cd[0]);
	for(int i = 1; i <= n; ++i){
		b[i] = -abs(cd[i]);
	}

	// Aberthの初期値の半径をニュートン法で算出
	int m = 0; // 係数cの中で0でないものの個数
	for(int i = 0; i <= n; ++i) m += (fabs(c[i].real()) > 1e-6 ? 1 : 0);
	double rmax = 0.0; // 半径の最大値
	for(int i = 1; i <= n; ++i){
		double ri = pow(m*fabs(c[i].real())/fabs(c[0].real()), 1.0/(i+1.0));
		if(ri > rmax) rmax = ri;
	}
	cout << "r_max = " << rmax << endl;
	double r = rmax;
	newton(b, n, r, max_iter, eps);
	cout << "r = " << r << endl;
	r = 10;

	// Aberthの初期値
	complexf zc = -c[1]/(c[0]*(double)n);
	for(int j = 0; j < n; ++j){
		double theta = (2*RX_PI/(double)n)*j+RX_PI/(2.0*n);
		z[j] = zc+r*complexf(cos(theta), sin(theta));
		cout << "x" << j << "(0) = " << z[j] << (j == n-1 ? "" : ", ");
	}
	cout << endl;
	

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

	cout << "0, ";
	for(int i = 0; i < n; ++i) cout << z[i].real() << ", " << z[i].imag() << ", ";
	cout << endl;

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

		cout << k+1 << ", ";
		for(int i = 0; i < n; ++i) cout << z[i].real() << ", " << z[i].imag() << (i == n-1 ? "" : ", ");
		cout << endl;

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
	int n;	// 代数方程式の次数

	// 代数方程式の係数を読み込む
	if(ReadAlgebra("algebra.txt", ",", c0)){
		n = (int)c0.size()-1;
		c.resize(n+1);
		for(int i = 0; i <= n; ++i){
			c[i] = complex<double>(c0[i]);
		}
	}
	else{
		cout << "file read error!" << endl;
		return 1;
	}

	// (確認用)元の方程式の表示
	cout << "f(x) = ";
	for(int i = 0; i < n; ++i){
		cout << c0[i] << "x^" << n-i << " + ";
	}
	cout << c0[n] << " = 0" << endl;
	cout << endl;


	int max_iter = 100;
	double eps = 1e-6;

	// Aberthの初期値
	vector<complexf> z(n);
	aberth(z, c, n, max_iter, eps);

	// ワイヤストラス法(DK公式)で解を求める
	weierstrass(z, c, n, max_iter, eps);

	// 結果の画面表示
	cout << "solutions : " << endl;
	for(int i = 0; i < n; ++i){
		cout << "x" << i << " = " << z[i];

		// 解の精度のチェック(f(x)を計算してみる)
		complexf f = c0[n], zi = z[i];
		for(int j = n-1; j >= 0; --j){
			f += c0[j]*zi; zi *= z[i];
		}
		cout << "  --> f = " << f;
		if(sqrt(f.real()*f.real()+f.imag()*f.imag()) < 1e-6) cout << " (OK)";
		cout << endl;
	}
	cout << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;



	return 0;
}


