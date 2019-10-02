/*! 
  @file error.chpp
	
  @brief 数値の精度と誤差
 
  @author Makoto Fujisawa
  @date 2019-09
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
//#include <float.h>


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// 整数の範囲
	cout << "char : " << CHAR_MIN << " - " << CHAR_MAX << "  (" << sizeof(char)*8 << "bit)" << endl;
	cout << "short : " << SHRT_MIN << " - " << SHRT_MAX << "  (" << sizeof(short)*8 << "bit)" << endl; // 正確にはshort int
	cout << "int : " << INT_MIN << " - " << INT_MAX << "  (" << sizeof(int)*8 << "bit)" << endl;
	cout << "longlong : " << LLONG_MIN << " - " << LLONG_MAX << "  (" << sizeof(long long)*8 << "bit)" << endl;
	//cout << "float : " << FLT_MIN << " - " << FLT_MAX << "  (" << sizeof(float)*8 << "bit)" << endl;
	//cout << "double : " << DBL_MIN << " - " << DBL_MAX << "  (" << sizeof(double)*8 << "bit)" << endl;
	cout << endl;

	// float型の有効桁数
	float x = 123456789;
	float y = 123456700;
	float z = x-y;
	cout.precision(10); // 表示桁数を10桁にする
	cout << "x = " << x << endl;
	cout << "y = " << y << endl;
	cout << "x-y = " << z << endl;
	cout << endl;

	// 丸め誤差
	//cout << "\n[round-off error]" << endl;
	float a = 1.1;
	cout.precision(20); // 表示する桁数を20桁に設定
	cout << "a = " << a << endl; // 1.1ぴったりになるはずだけど...
	cout << endl;

	// 丸め誤差による計算への影響
	cout << "\n[effect by round-off error]" << endl;
	float b = a*3;
	if(b == 3.3) cout << "b is equal to 3.3" << endl;
	else cout << "b is not equal to 3.3" << endl;
	cout << endl;

	// 演算で発生する誤差
	cout << "\n[effect by round-off error 2]" << endl;
	float a1 = 1234567;
	float a2 = 0.00123;
	float a3 = 1234567;
	cout << "(a1+a2)-a3=" << (a1+a2)-a3 << endl;
	cout << "(a1-a3)+a2=" << (a1-a3)+a2 << endl;

	// 倍精度と単精度の違い
	cout << "\n[single/double precision]" << endl;
	float x1 = 1.1;
	double x2 = 1.1;
	cout << "x1 = " << x1 << "  (single precision - " << sizeof(x1)*8 << "bit)" << endl;
	cout << "x2 = " << x2 << "  (double precision - " << sizeof(x2)*8 << "bit)" << endl;

	// 打ち切り誤差
	// ライプニッツの公式でπを計算
	cout << "\n[trancation error]" << endl;
	int n = 1000000;
	double pi = 0;
	int sgn = 1;
	for(int i = 0; i <= n; ++i){
		pi += sgn/(2.0*i+1.0);
		sgn *= -1;
		if(i%(n/10) == 0) cout << "n=" << i << " : " << 4*pi << endl;
	}
	double pi0 = 3.141592653589793; // 真値
	cout << "error = " << 4*pi-pi0 << endl;

	// 真値が分からない場合の収束判定
	cout << "\n[convergence test]" << endl;
	double eps = 1e-5;
	pi = 0;
	sgn = 1;
	for(int i = 0; i <= 1000000; ++i){
		double pi0 = pi;	// 収束判定のために前の反復の値を確保しておく
		pi += sgn/(2.0*i+1.0)*4;
		sgn *= -1;
		if(fabs(pi-pi0) <= eps) break; // 収束判定
	}
	cout << "pi = " << pi << endl;


	return 0;
}


