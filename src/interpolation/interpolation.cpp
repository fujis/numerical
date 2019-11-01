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
 * @return 補間値
 */
double linear_interpolation(const vector<double> &f, double a, double b, double x)
{
	double l = b-a;
	return ((x-a)*f[1]+(b-x)*f[0])/l;
}


/*!
 * ラグランジュ補間
 *  - 線形補間は2点の時のラグランジュ補間
 * @param[in] fi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列
 * @param[in] n データ数
 * @param[in] x 補間した値が必要な位置x
 * @return 位置xでの補間値
 */
double lagrangian_interpolation(const vector<double> &fi, const vector<double> &xi, int n, double x)
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
	return Ln;
}


/*!
 * サンプリング点(データ点)の生成(1次元)
 *  - チェビシェフ節点
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] dx サンプリング間隔
 * @param[in] func 関数値を与える関数ポインタ
 * @param[out] xi,yi サンプリングデータ
 * @return 生成されたデータ個数
 */
static int MakeChebyshevNodes(double x0, double x1, double dx, double func(double), vector<double> &xi, vector<double> &yi)
{
	xi.clear(); yi.clear();
	int cnt = 0;
	int n = (x1-x0)/dx+1;
	for(int i = n; i >= 1; --i){
		double x = cos((2.0*i-1.0)/(2.0*n)*RX_PI);
		x = x0+(x/2.0+0.5)*(x1-x0);
		xi.push_back(x);
		yi.push_back(func(x));
	}
	return n;
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

	// 補間用の点
	vector<double> xi, yi;
	MakeSamplingPoints(x0, x1, x1-x0, func, xi, yi);
	OutputSamplingPoints(xi, yi); // サンプリング点の画面出力

	// 線形補間
	fx = linear_interpolation(yi, xi[0], xi[1], x);
	cout << "f_linear(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;

	// 2点ラグランジュ補間(=線形補間)
	fx = lagrangian_interpolation(yi, xi, xi.size(), x);
	cout << "f_lagrangian(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	cout << endl;


	// 補間用の点を増やす(6点⇒5区間なので5次多項式で補間)
	MakeSamplingPoints(x0, x1, (x1-x0)/5.0, func, xi, yi);
	OutputSamplingPoints(xi, yi); // サンプリング点の画面出力

	// 6点ラグランジュ補間(5次多項式補間)
	fx = lagrangian_interpolation(yi, xi, xi.size(), x);
	cout << "f_lagrangian(" << x << ") = " << fx << ",  error = " << fabs(fx-gt) << endl;
	cout << "ground truth = " << gt << endl;



	// グラフ描画用にデータ出力
	int m = 6; // データ点数(次数はデータ点数-1)
	MakeSamplingPoints(x0, x1, (x1-x0)/(m-1.0), func, xi, yi);
	OutputSamplingPoints(xi, yi, "dat/lagrangian"+TOSTR(m)+"_data.txt");
	OutputFunction(x0, x1, (x1-x0)/200, std::bind(lagrangian_interpolation, yi, xi, xi.size(), std::placeholders::_1), "dat/lagrangian"+TOSTR(m)+".txt");
	// チェビシェフ節点を用いる場合
	MakeChebyshevNodes(x0, x1, (x1-x0)/(m-1.0), func, xi, yi);
	OutputSamplingPoints(xi, yi, "dat/lagrangian"+TOSTR(m)+"c_data.txt");
	OutputFunction(x0, x1, (x1-x0)/200, std::bind(lagrangian_interpolation, yi, xi, xi.size(), std::placeholders::_1), "dat/lagrangian"+TOSTR(m)+"c.txt");

	// 真値のグラフ作成用ファイル出力
	OutputFunction(x0, x1, (x1-x0)/200, func, "dat/lagrangian_ground_truth.txt");


	return 0;
}


