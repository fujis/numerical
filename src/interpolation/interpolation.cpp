/*! 
  @file interpolation.chpp
	
  @brief 補間法
 
  @author Makoto Fujisawa
  @date 2019-08
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"


//-----------------------------------------------------------------------------
// 補間法
//-----------------------------------------------------------------------------
/*!
 * 線形補間
 * @param[in] f 関数値(2つ)
 * @param[in] a,b 関数値fに対応する位置x
 * @param[in] x 補間した値が必要な位置x
 * @param[out] ans 解
 * @return 
 */
int linear_interpolation(const vector<double> &f, double a, double b, double x, double &ans)
{
	double l = b-a;
	ans = ((x-a)*f[1]+(b-x)*f[0])/l;
	return 0;
}


/*!
 * ラグランジュ補間
 *  - 線形補間は2点の時のラグランジュ補間
 * @param[in] fi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列
 * @param[in] n データ数
 * @param[in] x 補間した値が必要な位置x
 * @param[out] f 位置xでの関数値
 * @return
 */
int lagrangian_interpolation(const vector<double> &fi, const vector<double> &xi, int n, double x, double &f)
{
	double Ln = 0.0;
	for(int i = 0; i < n; ++i){
		// 補間係数(Π(x-xj)/(xi-xj) (j!=i) の計算)
		double l = 1;
		for(int j = 0; j < n; ++j){
			if(j == i) continue;
			l *= (x-xi[j])/(xi[i]-xi[j]);
		}
		// 補間値
		Ln += fi[i]*l;
	}
	f = Ln;
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
	xi.push_back(1.0);
	yi.push_back(func(xi.back()));

	cout << "sampling points : ";
	for(int i = 0; i < xi.size(); ++i){
		cout << "(" << xi[i] << ", " << yi[i] << ")" << (i == xi.size()-1 ? "" : ",  ");
	}
	cout << endl;

	// 線形補間
	linear_interpolation(yi, xi[0], xi[1], x, fx);
	cout << "f_linear(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	// 2点ラグランジュ補間(=線形補間)
	lagrangian_interpolation(yi, xi, 2, x, fx);
	cout << "f_lagrangian(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	// 補間用の点の追加
	xi.push_back(0.33);
	yi.push_back(func(xi.back()));
	xi.push_back(0.66);
	yi.push_back(func(xi.back()));

	cout << "sampling points : ";
	for(int i = 0; i < xi.size(); ++i){
		cout << "(" << xi[i] << ", " << yi[i] << ")" << (i == xi.size()-1 ? "" : ",  ");
	}
	cout << endl;

	// 4点ラグランジュ補間
	lagrangian_interpolation(yi, xi, 4, x, fx);
	cout << "f_lagrangian(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	cout << "ground truth = " << gt << endl;

		
	return 0;
}


