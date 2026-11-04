/*!
@file penalty.cpp

@brief ペナルティ関数法による制約付き最小化問題
		- 内点ペナルティ関数法(バリア関数法) -> コード中ではLogBarrierFunc, InverseBarrierFunc
		- 外点ペナルティ関数法 -> コード中ではPenaltyFunc

@author Makoto Fujisawa
@date 2026-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


const double SQRT5 = 2.23606797749979;


//-----------------------------------------------------------------------------
// 制約付き最小化問題における制約を与える関数の定義
//  - 等式制約の場合は h(x) = 0
//  - 不等式制約の場合は h(x) <= 0
//  - h(x) >= 0 の不等式制約を掛けたい場合は返値にマイナスを付けて，-h(x) <= 0とする
//-----------------------------------------------------------------------------
/*!
 * 制約条件h(x)とその勾配∇h(x) : h(x) = x+y-1 (等式制約 h(x) = 0 として使用)
 * @param[in] x 変数
 * @return 値or勾配
 */
inline double Constraint1(const vector<double> &x)
{
	return x[0]+x[1]-1.0;
}
inline vector<double> DConstraint1(const vector<double> &x)
{
	vector<double> dg(x.size(), 0.0);
	dg[0] = 1.0; dg[1] = 1.0;
	return dg;
}

/*!
 * 制約条件h(x)とその勾配∇h(x) : h(x) = x+y-0.5 (不等式制約 h(x) <= 0 として使用)
 * @param[in] x 変数
 * @return 値or勾配
 */
inline double Constraint2(const vector<double> &x)
{
	return x[0]+x[1]-0.5;
}
inline vector<double> DConstraint2(const vector<double> &x)
{
	vector<double> dg(x.size(), 0.0);
	dg[0] = 1.0; dg[1] = 1.0;
	return dg;
}

/*!
 * 制約条件h(x)とその勾配∇h(x) : h(x) = -x (不等式制約 h(x) <= 0 すなわち x >= 0 として使用)
 * @param[in] x 変数
 * @return 値or勾配
 */
inline double Constraint3(const vector<double>& x)
{
	return -x[0];
}
inline vector<double> DConstraint3(const vector<double>& x)
{
	vector<double> dg(x.size(), 0.0);
	dg[0] = -1.0; dg[1] = 0.0;
	return dg;
}


// 目的関数と制約条件のための関数オブジェクトの型(C++ならstd::functionを使う方が良い)
typedef double (*FUNCTION)(const vector<double>&);
typedef vector<double> (*DFUNCTION)(const vector<double>&);

//! 目的関数と制約をまとめた構造体
struct ObjectiveAndConstraints
{
	FUNCTION objectivefunc;				//!< 目的関数の値を返す関数ポインタ
	DFUNCTION objectivegrad;			//!< 目的関数の勾配値を返す関数ポインタ
	vector<FUNCTION> constraintfuncs;	//!< 制約関数の値を返す関数ポインタ
	vector<DFUNCTION> constraintgrads;	//!< 制約関数の勾配値を返す関数ポインタ
	vector<int> constraintequalities;	//!< 制約が等式(1)か不等式(0)を格納する配列
	int num_constraint;					//!< 制約の数
	ObjectiveAndConstraints() : objectivefunc(0), objectivegrad(0), num_constraint(0){}

	//! 制約の追加(equality = 1:等式制約 h(x) = 0, 0:不等式制約 h(x) <= 0)
	void AddConstraint(FUNCTION func, DFUNCTION grad, int equality)
	{
		constraintfuncs.push_back(func);
		constraintgrads.push_back(grad);
		constraintequalities.push_back(equality);
		num_constraint = (int)constraintfuncs.size();
	}
};

/*!
 * 制約の違反量 Σ|h_i(x)| (等式制約) + Σmax(0,h_j(x)) (不等式制約)
 * @param[in] oac 目的関数と制約関数
 * @param[in] x 変数
 */
double Violation(const ObjectiveAndConstraints& oac, const vector<double>& x)
{
	double v = 0.0;
	for(int j = 0; j < oac.num_constraint; ++j){
		double h = oac.constraintfuncs[j](x);
		if(oac.constraintequalities[j] == 1) v += fabs(h);	// 等式制約
		else if(h > 0.0) v += h;							// 不等式制約
	}
	return v;
}


//-----------------------------------------------------------------------------
// 内点ペナルティ関数法 : バリア関数を加えた目的関数 F(x) = f(x)+b(x)
//  - LogBarrierFuncやDLogBarrierFuncを既存の最適化問題の関数に関数ポインタとして渡すために，
//    バリアパラメータや制約関数などはグローバル変数で与える形にしている
//    (外点ペナルティ関数法のPenaltyFunc,DPenaltyFuncも同様)
//  - g_oacはbarrierfunction()やpenaltyfunction()の中で引数のoacから設定する
//-----------------------------------------------------------------------------
double g_mu = 1.0;				//!< バリアパラメータμ
ObjectiveAndConstraints g_oac;	//!< 目的関数と制約関数

/*!
* ログバリア関数を加えた目的関数
*  - ログバリア関数 : b(x) = -μ Σ log(-h_j(x))
*  - 制約範囲外(h_j(x) >= 0)では大きな値を返す(探索点が範囲外に出ないように)
* @param[in] x 変数
*/
double LogBarrierFunc(const vector<double> &x)
{
	// 目的関数の値
	double f = g_oac.objectivefunc(x);

	// バリア関数の値を目的関数の値に加える
	for(int j = 0; j < g_oac.num_constraint; ++j){
		double h = g_oac.constraintfuncs[j](x);
		if(h >= 0.0) return 1e30; // ゼロ割やlog(0)を防ぐためのエラー処理
		f -= g_mu*log(-h);
	}
	return f;
}
/*!
* ログバリア関数を加えた目的関数の勾配値
*  - ログバリア関数の勾配:　∇b(x) = -μ Σ ∇h_j(x)/h_j(x)
* @param[in] x 変数
*/
vector<double> DLogBarrierFunc(const vector<double>& x)
{
	// 目的関数の勾配値
	vector<double> df = g_oac.objectivegrad(x);

	// バリア関数の勾配値を目的関数の勾配値に加える
	for (int j = 0; j < g_oac.num_constraint; ++j) {
		double h = g_oac.constraintfuncs[j](x);
		vector<double> dh = g_oac.constraintgrads[j](x);
		double c = -g_mu / h;
		for (size_t k = 0; k < df.size(); ++k) df[k] += c * dh[k];
	}
	return df;
}

/*!
* 逆数バリア関数を加えた目的関数
*  - 逆数バリア関数 : b(x) = μ Σ 1/h_j(x)^2
*  - 制約範囲外(h_j(x) >= 0)では大きな値を返す(探索点が範囲外に出ないように)
* @param[in] x 変数
*/
double InverseBarrierFunc(const vector<double>& x)
{
	// 目的関数の値
	double f = g_oac.objectivefunc(x);

	// バリア関数の値を目的関数の値に加える
	for (int j = 0; j < g_oac.num_constraint; ++j) {
		double h = g_oac.constraintfuncs[j](x);
		if (h >= 0.0) return 1e30; // ゼロ割を防ぐためのエラー処理
		f += g_mu / (h * h);
	}
	return f;
}
/*!
* 逆数バリア関数を加えた目的関数の勾配値
*  - 逆数バリア関数の勾配:　∇b(x) = -2μ Σ ∇h_j(x)/h_j(x)^3
* @param[in] x 変数
*/
vector<double> DInverseBarrierFunc(const vector<double>& x)
{
	// 目的関数の勾配値
	vector<double> df = g_oac.objectivegrad(x);

	// バリア関数の勾配値を目的関数の勾配値に加える
	for (int j = 0; j < g_oac.num_constraint; ++j) {
		double h = g_oac.constraintfuncs[j](x);
		vector<double> dh = g_oac.constraintgrads[j](x);
		double c = -2.0 * g_mu / (h * h * h);
		for (size_t k = 0; k < df.size(); ++k) df[k] += c * dh[k];
	}
	return df;
}


//-----------------------------------------------------------------------------
// 外点ペナルティ関数法 : ペナルティ関数を加えた目的関数
//  F(x) = f(x)+ρ(Σ|h_i(x)|^α + Σ(max(0,h_j(x)))^β)  (h_i:等式制約, h_j:不等式制約)
//  - α = β = 2 とする(微分可能にするため．α = β = 1 では制約の境界で微分できない)
//-----------------------------------------------------------------------------
double g_rho = 1.0;		//!< ペナルティパラメータρ

/*!
* ペナルティ関数を加えた目的関数
*  - 制約範囲内では元の目的関数f(x)と同じ
* @param[in] x 変数
*/
double PenaltyFunc(const vector<double> &x)
{
	// 目的関数の値
	double f = g_oac.objectivefunc(x);

	// ペナルティ関数の値を目的関数の値に加える
	for (int j = 0; j < g_oac.num_constraint; ++j) {
		if (g_oac.constraintequalities[j] == 1) { // 等式制約
			double h = g_oac.constraintfuncs[j](x);
			f += g_rho * h * h;
		}
		else if (g_oac.constraintequalities[j] == 0) { // 不等式制約
			double h = g_oac.constraintfuncs[j](x);
			if (h > 0.0) f += g_rho * h * h;
		}
	}
	return f;
}
/*!
* ペナルティ関数を加えた目的関数の勾配
*  - ∇F(x) = ∇f(x) + 2ρ Σ h_i(x)∇h_i(x) + 2ρ Σ max(0,h_j(x))∇h_j(x)
* @param[in] x 変数
*/
vector<double> DPenaltyFunc(const vector<double> &x)
{
	// 目的関数の勾配値
	vector<double> df = g_oac.objectivegrad(x);

	// ペナルティ関数の勾配値を目的関数の勾配値に加える
	for (int j = 0; j < g_oac.num_constraint; ++j) {
		if (g_oac.constraintequalities[j] == 1) { // 等式制約
			double h = g_oac.constraintfuncs[j](x);
			vector<double> dh = g_oac.constraintgrads[j](x);
			for (size_t k = 0; k < df.size(); ++k) df[k] += 2.0 * g_rho * h * dh[k];
		}
		else if (g_oac.constraintequalities[j] == 0) { // 不等式制約
			double h = g_oac.constraintfuncs[j](x);
			if (h <= 0.0) continue;
			vector<double> dh = g_oac.constraintgrads[j](x);
			for (size_t k = 0; k < df.size(); ++k) df[k] += 2.0 * g_rho * h * dh[k];
		}
	}
	return df;
}


//-----------------------------------------------------------------------------
// 制約なしの最適化関数(最小値探索) : quasinewton.cpp と同じ
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
* BFGS法(準ニュートン法の一種, Broyden–Fletcher–Goldfarb–Shanno method)
*  - quasinewton.cppとほぼ同じ(各反復での画面出力を削除し，直線探索の範囲を制約範囲内に制限)
* @param[in] func 関数値を与える関数ポインタ
* @param[in] dfunc 関数勾配値を与える関数ポインタ
* @param[in] x0 初期探索地点
* @param[out] xout 解
* @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
* @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
* @return 1:成功,0:失敗
*/
int quasinewton_bfgs(double func(const vector<double>&), vector<double> dfunc(const vector<double>&), vector<double> x0, vector<double> &xout, int &max_iter, double &eps)
{
	int n = x0.size();
	vector< vector<double> > H; // ヘッセ行列の逆行列H^-1
	H.resize(n);
	for(int i = 0; i < n; ++i) H[i].resize(n, 0.0);
	for(int i = 0; i < n; ++i) H[i][i] = 1.0;

	vector<double> x = x0;
	vector<double> g = dfunc(x);
	double alpha = 0.1, gnorm = 1.0;

	int k;
	for(k = 0; k < max_iter; ++k){
		// 探索方向dの計算(d = -H∇f)
		vector<double> p(n, 0.0);
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				p[i] -= H[i][j]*g[j];
			}
		}

		// 探索範囲[0,tmax]の端が制約範囲外にならないようにtmaxを小さくする(バリア関数法のため)
		// (範囲外ではfuncが大きな値を返すので，両端が範囲外だと黄金分割法がうまく動かない)
		double tmax = 3.0;
		vector<double> xt(n, 0.0);
		for(int l = 0; l < 60; ++l){
			for(int i = 0; i < n; ++i) xt[i] = x[i]+tmax*p[i];
			if(func(xt) < 1e30) break;
			tmax *= 0.5;
		}

		// 黄金分割法で探索方向d上の最小値までの距離を調べてalphaとして設定
		int gmax_iter = 40;
		double geps = 1e-6*tmax;
		goldensection(func, x, p, 0.0, tmax, alpha, gmax_iter, geps);

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
	eps = sqrt(gnorm);

	return 0;
}


//-----------------------------------------------------------------------------
// 制約付きの最適化関数
//-----------------------------------------------------------------------------
/*!
* 内点ペナルティ関数法(バリア関数法, barrier function method)
*  1. 制約範囲内の内点を初期点x^(0)とし，バリアパラメータの初期値μ^(0)>0を設定
*  2. k = 1～k_maxまで以下を反復
*     a. f(x)+b(x)を目的関数として最小解x^(k)を求める(準ニュートン法．初期点はx^(k-1))
*     b. μ^(k) < ε ならばx^(k)を解として反復終了
*     c. μ^(k+1) = r μ^(k) (0<r<1)
* @param[in] oac 目的関数と制約関数(不等式制約のみ)
* @param[in] x0 初期探索地点(制約範囲の内点)
* @param[in] mu0 バリアパラメータの初期値
* @param[in] r バリアパラメータを小さくする割合(0<r<1)
* @param[out] xout 解
* @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
* @param[inout] eps 許容誤差(反復終了後,最後のμを返す)
* @param[in] inverse falseならログバリア関数，trueなら逆数バリア関数を用いる
* @return 1:成功,0:失敗
*/
int barrierfunction(const ObjectiveAndConstraints& oac, vector<double> x0, double mu0, double r, vector<double> &xout, int &max_iter, double &eps, bool inverse = false)
{
	// 等式制約があるかをチェック(等式制約では制約範囲の内点がないので，内点ペナルティ関数法(バリア関数法)は使えない)
	int num_equalities = 0;
	for (int i = 0; i < oac.num_constraint; ++i) num_equalities += (oac.constraintequalities[i] == 1 ? 1 : 0);
	if(num_equalities > 0){
		cout << "barrier function method cannot handle equality constraints" << endl;
		return 0;
	}
	// 初期点が内点かどうかの確認
	for (int i = 0; i < oac.num_constraint; ++i){
		double h = oac.constraintfuncs[i](x0);
		if(h >= 0.0){
			cout << "initial point is not an interior point (h" << i << " = " << h << ")" << endl;
			return 0;
		}
	}

	// BFGSに渡す関数(LogBarrierFuncなど)で使う目的関数と制約関数を設定
	g_oac = oac;

	// 使用するバリア関数
	FUNCTION func = (inverse ? InverseBarrierFunc : LogBarrierFunc);
	DFUNCTION dfunc = (inverse ? DInverseBarrierFunc : DLogBarrierFunc);

	vector<double> x = x0;
	g_mu = mu0;

	int k;
	for(k = 1; k <= max_iter; ++k){
		// f(x)+b(x)を準ニュートン法で最小化(初期点は前回の解)
		int qmax_iter = 100;
		double qeps = 1e-6;
		quasinewton_bfgs(func, dfunc, x, x, qmax_iter, qeps);

		cout << "k = " << k << ", mu = " << g_mu << " : x = (" << x << "), f = " << oac.objectivefunc(x);
		cout << " (BFGS iter = " << qmax_iter << ")" << endl;

		if(g_mu < eps) break;	// バリアパラメータが十分小さくなったら終了
		g_mu *= r;				// バリアパラメータを小さくする
	}

	xout = x;
	max_iter = k;
	eps = g_mu;

	return 1;
}

/*!
* 外点ペナルティ関数法(exterior penalty function method)
*  1. 初期点x^(0)(制約範囲外でもよい)とペナルティパラメータの初期値ρ^(0)>0を設定
*  2. k = 1～k_maxまで以下を反復
*     a. f(x)+ρ(...)を目的関数として最小解x^(k)を求める(準ニュートン法．初期点はx^(k-1))
*     b. 制約の違反量が ε 未満ならばx^(k)を解として反復終了
*     c. ρ^(k+1) = r ρ^(k) (r>1)
* @param[in] oac 目的関数と制約関数(等式制約，不等式制約のどちらでもよい)
* @param[in] x0 初期探索地点
* @param[in] rho0 ペナルティパラメータの初期値
* @param[in] r ペナルティパラメータを大きくする割合(r>1)
* @param[out] xout 解
* @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
* @param[inout] eps 許容誤差(反復終了後,最後の制約の違反量を返す)
* @return 1:成功,0:失敗
*/
int penaltyfunction(const ObjectiveAndConstraints& oac, vector<double> x0, double rho0, double r, vector<double> &xout, int &max_iter, double &eps)
{
	// BFGSに渡す関数(PenaltyFunc,DPenaltyFunc)で使う目的関数と制約関数を設定
	g_oac = oac;

	vector<double> x = x0;
	g_rho = rho0;

	double v = Violation(oac, x);
	int k;
	for(k = 1; k <= max_iter; ++k){
		// f(x)+ペナルティ関数を準ニュートン法で最小化(初期点は前回の解)
		int qmax_iter = 100;
		double qeps = 1e-6;
		quasinewton_bfgs(PenaltyFunc, DPenaltyFunc, x, x, qmax_iter, qeps);

		v = Violation(oac, x);
		cout << "k = " << k << ", rho = " << g_rho << " : x = (" << x << "), f = " << oac.objectivefunc(x);
		cout << ", violation = " << v << " (BFGS iter = " << qmax_iter << ")" << endl;

		if(v < eps) break;	// 制約の違反量が十分小さくなったら終了
		g_rho *= r;			// ペナルティパラメータを大きくする
	}

	xout = x;
	max_iter = k;
	eps = v;

	return 1;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	int max_iter = 100;
	double eps = 1e-6;

	vector<double> x(2, 0.0);
	cout.precision(8);

	// 問題1 : 不等式制約
	{
		cout << "argmin x^2+y^2-xy-x-1  s.t. x+y-0.5 <= 0, -x <= 0" << endl;
		cout << "exact : (x,y) = (" << 5.0 / 12.0 << "," << 1.0 / 12.0 << "), f(x,y) = " << -183.0 / 144.0 << endl << endl;

		ObjectiveAndConstraints oac;
		oac.objectivefunc = Func4;
		oac.objectivegrad = DFunc4;
		oac.AddConstraint(Constraint2, DConstraint2, 0);	// x+y-0.5 <= 0
		oac.AddConstraint(Constraint3, DConstraint3, 0);	// -x <= 0

		// 内点ペナルティ関数法(バリア関数法) : ログバリア関数
		cout << "[barrier function method (interior penalty) : log barrier]" << endl;
		vector<double> x0(2, 0.1);	// 初期値(制約範囲の内点)
		max_iter = 100; eps = 1e-6;
		barrierfunction(oac, x0, 1.0, 0.1, x, max_iter, eps, false);
		cout << "  (x,y) = (" << x[0] << "," << x[1] << "), f(x,y) = " << oac.objectivefunc(x) << ", iter = " << max_iter << endl << endl;

		// 内点ペナルティ関数法(バリア関数法) : 逆数バリア関数
		cout << "[barrier function method (interior penalty) : inverse barrier]" << endl;
		x0 = vector<double>(2, 0.1);	// 初期値(制約範囲の内点)
		max_iter = 100; eps = 1e-6;
		barrierfunction(oac, x0, 1.0, 0.1, x, max_iter, eps, true);
		cout << "  (x,y) = (" << x[0] << "," << x[1] << "), f(x,y) = " << oac.objectivefunc(x) << ", iter = " << max_iter << endl << endl;

		// 外点ペナルティ関数法
		cout << "[exterior penalty method]" << endl;
		x0 = vector<double>(2, 1.0);	// 初期値(制約範囲外でもよい)
		max_iter = 100; eps = 1e-6;
		penaltyfunction(oac, x0, 1.0, 10.0, x, max_iter, eps);
		cout << "  (x,y) = (" << x[0] << "," << x[1] << "), f(x,y) = " << oac.objectivefunc(x) << ", iter = " << max_iter << endl << endl;
	}

	// 問題2 : 等式制約(外点ペナルティ関数法のみ)
	{
		cout << "argmin x^2+2y^2  s.t. x+y-1 = 0" << endl;
		cout << "exact : (x,y) = (" << 2.0 / 3.0 << "," << 1.0 / 3.0 << "), f(x,y) = " << 2.0 / 3.0 << endl << endl;

		ObjectiveAndConstraints oac;
		oac.objectivefunc = Func6;
		oac.objectivegrad = DFunc6;
		oac.AddConstraint(Constraint1, DConstraint1, 1);	// x+y-1 = 0

		// 外点ペナルティ関数法
		cout << "[exterior penalty method]" << endl;
		vector<double> x0(2, 0.0);	// 初期値
		max_iter = 100; eps = 1e-6;
		penaltyfunction(oac, x0, 1.0, 10.0, x, max_iter, eps);
		cout << "  (x,y) = (" << x[0] << "," << x[1] << "), f(x,y) = " << oac.objectivefunc(x) << ", iter = " << max_iter << endl << endl;
	}

	return 0;
}
