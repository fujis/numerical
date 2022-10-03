/*!
@file quasinewton.cpp

@brief 準ニュートン法

@author Makoto Fujisawa
@date 2019-07
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


const double SQRT5 = 2.23606797749979;


//-----------------------------------------------------------------------------
// 最適化関数(最小値探索)
//-----------------------------------------------------------------------------
/*!
* 黄金分割探索法(golden-section method)
- 多変数で探索方向が決まっている場合
- BGFSでの探索方向に対する係数α決定のために用いる
* @param[in] func 関数値を与える関数ポインタ
* @param[in] x0 探索初期地点ベクトル
* @param[in] d 探索方向ベクトル
* @param[in] xl,xr 初期探索範囲
* @param[out] x 解
* @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
* @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
* @return 1:成功,0:失敗
*/
int goldensection(double func(const vector<double>&), const vector<double> x0, const vector<double> d, double xl, double xr, double &ans, int &max_iter, double &eps)
{
	const double eta = 2.0/(SQRT5+3.0);
	double beta = xr-xl;
	double tau = eta*beta;
	int n = x0.size();

	double x[4];
	x[0] = xl;
	x[1] = xl+tau;
	x[2] = xr-tau;
	x[3] = xr;

	// 探索方向上の2点x1,x2での関数値を計算
	vector<double> x1(n, 0.0), x2(n, 0.0);
	for(int i = 0; i < n; ++i){
		x1[i] = x0[i]+x[1]*d[i];
		x2[i] = x0[i]+x[2]*d[i];
	}
	double f1 = func(x1);
	double f2 = func(x2);

	int k;
	for(k = 0; k < max_iter; ++k){

		if(f1 < f2){ // 区間[x2,x3]に最小値なし
					 // 分割区間更新
			double x2 = x[2];
			x[2] = x[1];
			x[1] = x2-tau;
			x[3] = x2;

			// 関数値の更新
			for(int i = 0; i < n; ++i) x1[i] = x0[i]+x[1]*d[i];
			f2 = f1;
			f1 = func(x1);
		} else{ // 区間[x0,x1]に最小値なし
				// 分割区間更新
			double x1 = x[1];
			x[1] = x[2];
			x[2] = x1+tau;
			x[0] = x1;

			// 関数値の更新
			for(int i = 0; i < n; ++i) x2[i] = x0[i]+x[2]*d[i];
			f1 = f2;
			f2 = func(x2);
		}

		beta -= tau;	// 新しい[x0,x3]の長さ(現在の長さから[x0,x1](or[x2,x3])の長さtauを引く)
		tau = eta*beta; // 新しい[x0,x1](or[x2,x3])の長さ

		if(beta < eps) break;
	}

	ans = x[1];
	max_iter = k;
	eps = beta;

	return 0;
}


/*!
* BGFS法(準ニュートン法の一種, Broyden–Fletcher–Goldfarb–Shanno method)
* @param[in] func 関数値を与える関数ポインタ
* @param[in] dfunc 関数勾配値を与える関数ポインタ
* @param[in] x0 初期探索地点
* @param[out] xout 解
* @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
* @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
* @return 1:成功,0:失敗
*/
int quasinewton_bgfs(double func(const vector<double>&), vector<double> dfunc(const vector<double>&), vector<double> x0, vector<double> &xout, int &max_iter, double &eps)
{
	int n = x0.size();
	vector< vector<double> > H; // ヘッセ行列の逆行列H^-1
	H.resize(n);
	for(int i = 0; i < n; ++i) H[i].resize(n, 0.0);
	for(int i = 0; i < n; ++i) H[i][i] = 1.0;

	vector<double> x = x0;
	vector<double> g = dfunc(x);
	double alpha = 0.1, gnorm = 1.0;

	double f;
	int k;
	for(k = 0; k < max_iter; ++k){
		cout << "f(" << x << ") = " << func(x) << endl;
		//cout << x << ", " << gnorm << endl;

		// 探索方向dの計算(d = -H∇f)
		vector<double> p(n, 0.0);
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				p[i] -= H[i][j]*g[j];
			}
		}

		// 黄金分割法で探索方向d上の最小値までの距離を調べてalphaとして設定
		int gmax_iter = 20;
		double geps = 1e-4;
		goldensection(func, x, p, 0.0, 3.0, alpha, gmax_iter, geps);

		// xの更新
		vector<double> s(n, 0.0);
		for(int i = 0; i < n; ++i){
			s[i] = alpha*p[i];
			x[i] += s[i];
		}

		// yの計算
		vector<double> gprev = g;
		g = dfunc(x);
		vector<double> y(n, 0.0);
		for(int i = 0; i < n; ++i){
			y[i] = g[i]-gprev[i];
		}

		gnorm = 0.0;
		for(int i = 0; i < n; ++i) gnorm += g[i]*g[i];
		if(sqrt(gnorm) < eps) break; // 勾配ベクトル∇fの大きさで収束判定

		// H*yの計算(行列×縦ベクトル⇒縦ベクトル)
		vector<double> Hy(n, 0.0);
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				Hy[i] += H[i][j]*y[j];	
			}
		}
		// (s^T)*yと(y^T)*(H*y)の計算(どちらもベクトル同士の内積)
		double sy = 0.0, yHy = 0.0;
		for(int i = 0; i < n; ++i){
			sy += s[i]*y[i];
			yHy += y[i]*Hy[i];
		}

		// H^-1の更新
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				H[i][j] += s[i]*s[j]/sy + yHy*s[i]*s[j]/(sy*sy) - Hy[i]*s[j]/sy - s[i]*Hy[j]/sy;
			}
		}

	}

	xout = x;
	max_iter = k;
	eps = gnorm;

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	vector<double> x0(2, -1.0);
	vector<double> x(2, 0.0);

	int max_iter = 100;
	double eps = 1e-6;
	quasinewton_bgfs(Func4, DFunc4, x0, x, max_iter, eps);

	cout << "(x,y) = (" << x[0] << "," << x[1] << "), ";
	cout << " f = " << Func4(x) << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;


	return 0;
}


