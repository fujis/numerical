/*! 
  @file wave.cpp
	
  @brief 波動方程式のソルバー
 
  @author Makoto Fujisawa
  @date 2019-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

// 差分法用関数ポインタ
typedef int(*FUNC)(vector<double> &, vector<double> &, vector<double> &, double, double, int, double, double);
// 境界値設定用関数ポインタ
typedef void(*BC)(vector<double> &, int);


// 境界値
double alpha = 0.0, beta = 0.0;


//-----------------------------------------------------------------------------
//! 波動方程式のソルバー
//-----------------------------------------------------------------------------

/*!
 * 境界条件設定(ノイマン境界条件)
 * @param[inout] f 未知関数fの各グリッドでの値(初期値を入れておく)
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 */
void setbc_d(vector<double> &f, int n)
{
	f[0] = f[1];
	f[n] = f[n-1];
}

/*!
 * FTCS法(前進オイラー法+中心差分)による波動方程式のソルバー
 * @param[inout] f1 未知関数fの各グリッドでの値(ステップi+1での値)
 * @param[inout] f0 未知関数fの各グリッドでの値(ステップiでの値)
 * @param[in] c 波の速度を表す係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 * @return 問題なければ0を返す
 */
int wave1d_ftcs(vector<double> &f2, vector<double> &f1, vector<double> &f0, double c, double dt, int n, double x0, double xn)
{
	if(n < 0) return 1;
	double dx = (xn-x0)/n; // 空間刻み幅
	double eta = c*c*dt*dt/(dx*dx);
	for(int i = 1; i <= n-1; ++i){
		f2[i] = 2*f1[i]-f0[i]+eta*(f1[i+1]-2*f1[i]+f1[i-1]);
	}

	return 0;
}


/*!
 * 1ステップ進める
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[inout] f 未知関数fの各グリッドでの値
 * @param[in] c 波の速度を表す係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 */
void wave1d_step(FUNC sdfunc, BC bcfunc, vector<double> &f2, vector<double> &f1, double c, double dt, int n, double x0, double xn)
{
	vector<double> f0(f1.size());
	for(int i = 0; i < f2.size(); ++i){ f0[i] = f1[i]; f1[i] = f2[i]; }
	bcfunc(f0, n); bcfunc(f1, n);
	sdfunc(f2, f1, f0, c, dt, n, x0, xn);
	bcfunc(f2, n);
}

/*!
 * 時間ステップを進めていったときの結果をファイル出力するための関数
 * @param[in] tdfunc 時間差分を行う関数
 * @param[in] sdfunc 空間差分を行う関数
 * @param[in] bcfunc 境界条件設定用関数
 * @param[in] f 未知関数fの各グリッドでの値(元のデータに影響しないようにコピー渡しにしている)
 * @param[in] c 波の速度を表す係数
 * @param[in] dt 時間ステップ幅
 * @param[in] n 計算範囲内での分割数(h=(xn-x0)/n)
 * @param[in] x0,xn 計算範囲
 * @param[in] filename 出力ファイル名
 */
void wave1d(FUNC sdfunc, BC bcfunc, vector<double> f, double c, double dt, int n, double x0, double xn, const string &filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	fo << "#1d,-0.4,0.8" << endl;

	// 初期値の出力
	OutputValueToFile(f, n, x0, xn, 0, fo);

	vector<double> f0(f);

	// タイムステップを進める
	for(int k = 1; k < 4000; ++k){
		wave1d_step(sdfunc, bcfunc, f, f0, c, dt, n, x0, xn);

		// 結果の出力
		OutputValueToFile(f, n, x0, xn, k*dt, fo);
	}
	fo.close();

	cout << "wave1d : " << filename << endl;
}

//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double dt = 0.01;
	double x0 = 0.0, xn = 1.0;
	double c = 0.02;
	int n = 100;
	BC bcfunc = setbc_d;

	vector<double> f(n+1, 0.0);

	// 初期値設定
	for(int i = 0; i < n+1; ++i) f[i] = 0.0;
	f[n/2] = 0.5;

	bcfunc(f, n);

	// データ保存パス(ファイル名の最後の_ignはGitで除外ファイルにするためのもの)
	string path = "../../bin/data/";

	// 前進オイラー+中心差分
	wave1d(wave1d_ftcs, bcfunc, f, c, dt, n, x0, xn, path+"wave1d_ftcs_ign.txt");


	return 0;
}


