/*! 
  @file advection.cpp
	
  @brief 移流方程式のソルバー
 
  @author Makoto Fujisawa
  @date 2019-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

// 差分法用関数ポインタ
typedef int(*FUNC)(vector<double> &, vector<double> &, double, double, int, double, double);
// 境界値設定用関数ポインタ
typedef void(*BC)(vector<double> &, int);
// 時間ステップ用関数ポインタ
typedef void(*STEP)(FUNC, BC, vector<double> &, double, double, int, double, double);

// 境界値
double alpha = 0.0, beta = 0.0;


//-----------------------------------------------------------------------------
//! 移流法のソルバー
//-----------------------------------------------------------------------------

/*!
 * 中心差分でdt分だけ移流を進める
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
int advect1d_central(vector<double> &f1, vector<double> &f0, double u, double dt, int n, double x0, double xn)
{
	if(n <= 0) return 1;
	double h = (xn-x0)/n; // 空間刻み幅
	double nu = u*dt/h; // クーラン数
	for(int i = 1; i <= n-1; ++i){
		// ∂f/∂xの差分による近似
		double gx = (f0[i+1]-f0[i-1])/2.0;
		// fの値を更新
		f1[i] += -nu*gx;
	}
	return 0;
}

/*!
 * 風上差分でdt分だけ移流を進める
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
int advect1d_upwind(vector<double> &f1, vector<double> &f0, double u, double dt, int n, double x0, double xn)
{
	if(n <= 0) return 1;
	double h = (xn-x0)/n; // 空間刻み幅
	double nu = u*dt/h; // クーラン数
	for(int i = 1; i <= n-1; ++i){
		// ∂f/∂xの差分による近似
		double gx;
		if(u > 0){
			gx = f0[i]-f0[i-1];
		}
		else{
			gx = f0[i+1]-f0[i];
		}
		
		// fの値を更新
		f1[i] += -nu*gx;
	}
	return 0;
}

/*!
 * Lax-Wendroff法でdt分だけ移流を進める
 *  -  P. Lax and B.Wendroff, "Systems of conservation laws", Commun. Pure Appl Math. 13, pp.217-237, 1960.
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
int advect1d_lw(vector<double> &f1, vector<double> &f0, double u, double dt, int n, double x0, double xn)
{
	if(n <= 0) return 1;
	double h = (xn-x0)/n; // 空間刻み幅
	double nu = u*dt/h; // クーラン数
	for(int i = 1; i <= n-1; ++i){
		// ∂f/∂xの差分による近似
		double gx = (f0[i+1]-f0[i-1])/2.0-nu*(f0[i+1]-2*f0[i]+f0[i-1])/2.0;
		// fの値を更新
		f1[i] += -nu*gx;
	}
	return 0;
}

/*!
 * セミラグランジュ法でdt分だけ移流を進める
 *  -  P. Lax and B.Wendroff, "Systems of conservation laws", Commun. Pure Appl Math. 13, pp.217-237, 1960.
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
int advect1d_sl(vector<double> &f1, vector<double> &f0, double u, double dt, int n, double x0, double xn)
{
	if(n <= 0) return 1;
	double h = (xn-x0)/n; // 空間刻み幅
	for(int i = 1; i <= n-1; ++i){
		// 現在の位置からバックトレースした位置を求める
		double x = x0+i*h;
		x -= u*dt; // バックトレース

		// バックトレースした位置のグリッド情報計算
		int ib = (int)((x-x0)/h); // xが含まれるグリッド番号
		if(ib < 0) ib = 0;
		if(ib >= n) ib = n-1;
		double dx = (x-ib*h)/h; // グリッド位置からの距離を[0,1]で正規化
		
		// fの値を線形補間で求める
		f1[i] = (1-dx)*f0[ib]+dx*f0[ib+1];
	}
	return 0;
}

/*!
 * 境界条件設定(ディリクレ境界条件)
 * @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 */
void setbc_d(vector<double> &f, int n)
{
	f[0] = alpha;
	f[n] = beta;
}

/*!
 * 前進オイラー法で1ステップ進める
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 */
void advect1d_step_forwardeular(FUNC sdfunc, BC bcfunc, vector<double> &f, double u, double dt, int n, double x0, double xn)
{
	vector<double> f0(f);
	sdfunc(f, f0, u, dt, n, x0, xn);
	bcfunc(f, n);
}

/*!
 * ホイン法で1ステップ進める
 *  - 空間差分部分の計算切り替えのために(主にセミラグランジュ法のために)，
 *    f(k+1)=f(k)+(dt/2)(F(k)+F(k+1)) を f(k+1)=(f(k)+dt*F(k))/2+(f(k)+dt*F(k+1))/2として計算している
 *  - 手法比較のために無駄な処理やメモリを使っているので実際にRK4で計算する場合は元の式で計算した方が良い
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 */
void advect1d_step_heun(FUNC sdfunc, BC bcfunc, vector<double> &f, double u, double dt, int n, double x0, double xn)
{
	vector<double> f1(f), f2(f);

	// 前進オイラー法で仮のf(k+1)を求める
	sdfunc(f1, f, u, dt, n, x0, xn); bcfunc(f, n);

	// f(k)+dt*F(t(k+1),f(k+1))を計算
	sdfunc(f2, f1, u, dt, n, x0, xn); bcfunc(f2, n);

	// f(k+1)=(f(k)+dt*F(k))/2+(f(k)+dt*F(k+1))/2でfを更新
	for(int i = 1; i <= n-1; ++i){
		f[i] = f1[i]/2+f2[i]/2;
	}
	bcfunc(f, n);
}

/*!
 * RK4で1ステップ進める
 *  - 空間差分部分の計算切り替えのために(主にセミラグランジュ法のために)，
 *    f(i+1)=f(i)+(dt/6)(k1+k2+k3+k4) を f(i+1)=(fi+dt*k1)/6+(fi+dt*k2)/3+(fi+dt*k3)/3+(fi+dt*k4)/6 として計算している
 *  - 手法比較のために無駄な処理やメモリを使っているので実際にRK4で計算する場合は元の式で計算した方が良い
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 */
void advect1d_step_rk4(FUNC sdfunc, BC bcfunc, vector<double> &f, double u, double dt, int n, double x0, double xn)
{
	vector<double> k1(f), k2(f), k3(f), k4(f), fk(f), fk2(f), fk3(f);

	// fi+dt*k1の計算
	sdfunc(k1, f, u, dt, n, x0, xn); bcfunc(k1, n);

	// fi+(dt/2)*k1=fi+(dt/2)*Fの計算
	sdfunc(fk, f, u, dt/2, n, x0, xn); bcfunc(fk, n);
	// fi+dt*k2の計算
	sdfunc(k2, fk, u, dt, n, x0, xn); bcfunc(k2, n);

	// fi+(dt/2)*k2=fi+(dt/2)*F(t+dt/2,fi+(dt/2)*F)の計算
	sdfunc(fk2, fk, u, dt/2, n, x0, xn); bcfunc(fk2, n);
	// fi+dt*k3の計算
	sdfunc(k3, fk2, u, dt, n, x0, xn); bcfunc(k3, n);

	// fi+dt*k3=fi+dt*F(t+dt,fi+(dt/2)*F(t+dt/2,fi+(dt/2)*F))の計算
	sdfunc(fk3, fk2, u, dt, n, x0, xn); bcfunc(fk3, n);
	// fi+dt*k4の計算
	sdfunc(k4, fk3, u, dt, n, x0, xn); bcfunc(k4, n);

	// f(i+1)=(fi+dt*k1)/6+(fi+dt*k2)/3+(fi+dt*k3)/3+(fi+dt*k4)/6
	for(int i = 1; i <= n-1; ++i){
		f[i] = k1[i]/6+k2[i]/3+k3[i]/3+k4[i]/6;
	}
	bcfunc(f, n);
}


/*!
 * 時間ステップを進めていったときの結果をファイル出力するための関数
 * @param[in] tdfunc 時間差分を行う関数
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] u 速度(全体で一定,場所によって変える場合はuを配列にする)
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] a,b 計算範囲
 * @param[in] 出力ファイル名
 */
void advect1d(STEP tdfunc, FUNC sdfunc, BC bcfunc, vector<double> f, double u, double dt, int n, double x0, double xn, const string &filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	fo << "#1d,0.0,1.0" << endl;

	// 初期値の出力
	OutputValueToFile(f, n, x0, xn, 0, fo);

	// タイムステップを進める
	for(int k = 1; k < 750; ++k){
		tdfunc(sdfunc, bcfunc, f, u, dt, n, x0, xn);

		// 結果の出力
		OutputValueToFile(f, n, x0, xn, k*dt, fo);
	}
	fo.close();

	cout << "advect1d : " << filename << endl;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double x0 = 0.0, xn = 5.0; // 計算範囲(空間)
	int n = 128;	// 分割数(空間)
	double u = 0.75; // 移流速度
	double dx = (xn-x0)/n;
	double dt = 0.1*dx/u;	// 時間ステップ幅(時間方向の刻み幅)

	vector<double> f(n+1, 0.0);

	// 初期値設定
	SetValueRectangle(f, n, x0, xn);
	//SetValueSin(f, n, x0, xn);

	// データ保存パス
	string path = "../../bin/data/";

	// 前進オイラー+中心差分
	advect1d(advect1d_step_forwardeular, advect1d_central, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_fe+cen_ign.txt");
	//// 前進オイラー+風上差分
	//advect1d(advect1d_step_forwardeular, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_fe+upwind_ign.txt");
	//// 前進オイラー+Lax-Wendroff
	//advect1d(advect1d_step_forwardeular, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_fe+lw_ign.txt");
	//// 前進オイラー+セミラグランジュ法
	//advect1d(advect1d_step_forwardeular, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_fe+sl_ign.txt");
	//// ホイン+中心差分
	//advect1d(advect1d_step_heun, advect1d_central, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_heun+cen_ign.txt");
	//// ホイン+風上差分
	//advect1d(advect1d_step_heun, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_heun+upwind_ign.txt");
	//// ホイン+Lax-Wendroff
	//advect1d(advect1d_step_heun, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_heun+lw_ign.txt");
	//// ホイン+セミラグランジュ法
	//advect1d(advect1d_step_heun, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_heun+sl_ign.txt");
	//// RK4+中心差分
	//advect1d(advect1d_step_rk4, advect1d_central, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_rk4+cen_ign.txt");
	//// RK4+風上差分
	//advect1d(advect1d_step_rk4, advect1d_upwind, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_rk4+upwind_ign.txt");
	//// RK4+Lax-Wendroff
	//advect1d(advect1d_step_rk4, advect1d_lw, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_rk4+lw_ign.txt");
	//// RK4+セミラグランジュ法
	//advect1d(advect1d_step_rk4, advect1d_sl, setbc_d, f, u, dt, n, x0, xn, path+"advect1d_rk4+sl_ign.txt");

	//ofstream fo;
	//fo.open("../../bin/advect1d_heun+cen.txt", ios::out);
	//fo << "#1d,0.0,1.0" << endl;

	//// 初期値の出力
	//OutputValueToFile(f, n, x0, xn, 0, fo);

	//// タイムステップを進める
	//for(int k = 1; k < 1000; ++k){
	//	advect1d_step_heun(advect1d_central, setbc_d, f, u, dt, n, x0, xn);

	//	// 結果の出力
	//	OutputValueToFile(f, n, x0, xn, k*dt, fo);
	//}
	//fo.close();


	return 0;
}


