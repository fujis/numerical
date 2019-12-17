/*! @file rx_utils.h
	
	@brief 数値計算テストの共通ヘッダ
 
	@author Makoto Fujisawa
	@date   2019
*/


#ifndef _RX_UTILS_H_
#define _RX_UTILS_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <functional>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#define RX_USE_MM
#ifdef WIN32
    #include <windows.h>
 
    #ifdef RX_USE_MM
    #include <mmsystem.h>	// 時間計測用
    #pragma comment (lib, "winmm.lib")
    #endif
#else
	#include <ctime>
#endif


using namespace std;

 
//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------

// 時間計測用
#ifdef WIN32
    #ifdef RX_USE_MM
        #define RXTIME DWORD
        #define RXGETTIME timeGetTime
        #define RXTIME2SEC 1.0e-3
        //#define RXTIME2SEC 1.0
    #else
        #define RXTIME DWORD
        #define RXGETTIME GetTickCount
        #define RXTIME2SEC 1.0e-3
        //#define RXTIME2SEC 1.0
    #endif
#else
    #define RXTIME clock_t
    #define RXGETTIME clock
    #define RXTIME2SEC (1.0/CLOCKS_PER_SEC)
#endif


//--------------------------------------------------------------------
// 定数
//--------------------------------------------------------------------
// 円周率
const double RX_PI = 3.14159265358979323846;

//! 許容誤差
const double RX_FEQ_EPS = 1.0e-10;


//--------------------------------------------------------------------
// マクロ(テンプレート関数)
//--------------------------------------------------------------------
//! ゼロ判定
template<class T> 
inline bool RX_IS_ZERO(const T &x){ return (fabs(x) < RX_FEQ_EPS); }

//! 許容誤差を含めた等値判定
template<class T> 
inline bool RX_FEQ(const T &a, const T &b){ return (fabs(a-b) < RX_FEQ_EPS); }

//! 最大値判定(2値)
template<class T> 
inline T RX_MAX(const T &a, const T &b){ return ((a > b) ? a : b); }

//! 最大値判定(3値)
template<class T> 
inline T RX_MAX3(const T &a, const T &b, const T &c){ return ( (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c)); }

//! 最小値判定(2値)
template<class T> 
inline T RX_MIN(const T &a, const T &b){ return ((a < b) ? a : b); }

//! 最小値判定(3値)
template<class T> 
inline T RX_MIN3(const T &a, const T &b, const T &c){ return ( (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c)); }

//! 値のクランプ(クランプした値を返す)
template<class T> 
inline T RX_CLAMP(const T &x, const T &a, const T &b){ return ((x < a) ? a : (x > b) ? b : x); }

//! 1次元線型補間
template<class T>
inline T RX_LERP(const T &a, const T &b, const T &t){ return a + t*(b-a); }

//! degree -> radian の変換係数(pi/180.0)
const double RX_DEGREES_TO_RADIANS = 0.0174532925199432957692369076848;

//! radian -> degree の変換係数(180.0/pi)
const double RX_RADIANS_TO_DEGREES = 57.295779513082320876798154814114;

//! degree -> radian の変換
template<class T>
inline T RX_TO_RADIANS(const T &x){ return static_cast<T>((x)*RX_DEGREES_TO_RADIANS); }

//! radian -> degree の変換
template<class T>
inline T RX_TO_DEGREES(const T &x){ return static_cast<T>((x)*RX_RADIANS_TO_DEGREES); }

//-----------------------------------------------------------------------------
// 関数定義
//-----------------------------------------------------------------------------
template<class T>
inline string RX_TO_STRING(T x)
{
	stringstream ss;
	ss << x;
	return ss.str();
}

template<class T>
inline string TOSTR(T x)
{
	stringstream ss;
	ss << x;
	return ss.str();
}

//! string型に<<オペレータを設定
template<typename T>
inline string &operator<<(string &cb, const T &x)
{
    cb += RX_TO_STRING(x);
    return cb;
}


//! vector型に<<オペレータを設定
template<typename T>
inline std::ostream &operator<<(std::ostream &out, const vector<T> &x)
{
	for(size_t i = 0; i < x.size(); ++i){
		out << x[i] << (i == x.size()-1 ? "" : ", ");
	}
	return out;
}



//-----------------------------------------------------------------------------
// テキストファイル処理
//-----------------------------------------------------------------------------

/*!
 * stringからpos以降で最初の区切り文字までを抽出
 *  - もし，"(ダブルクオーテーション)で囲まれていたらその範囲を抽出
 * @param[in] src 元の文字列
 * @param[out] sub 抽出文字列
 * @param[in] sep 区切り文字
 * @param[in] pos 探索開始位置
 * @return 次の抽出開始位置(sepの後にスペースがあればそのスペースの後)
 */
inline size_t GetNextString(const string &src, string &sub, string sep, size_t pos)
{
	bool extracted = false;
	if(src[pos] == '\"'){	// ダブルクオーテーションのチェック
		size_t j = src.find("\"", pos+1);
		if(j != string::npos){
			sub = src.substr(pos+1, j-(pos+1));
			pos = j+1;
			extracted = true;
		}
	}
 
	size_t i = src.find(sep, pos);
	if(i == string::npos){		
		if(!extracted) sub = src.substr(pos, string::npos);
		return string::npos;
	}
	else{
		int cnt = 1;
		while(src[i+cnt] == ' '){	// sepの後のスペースを消す
			cnt++;
		}
		if(!extracted) sub = src.substr(pos, i-pos);
		return (i+cnt >= src.size() ? string::npos : i+cnt);
	}
}
 
/*!
 * stringから最初の区切り文字までを抽出
 * @param[in] src 元の文字列
 * @param[out] sub 抽出文字列
 * @param[in] sep 区切り文字
 * @return 次の抽出開始位置(sepの後にスペースがあればそのスペースの後)
 */
inline size_t GetFirstString(const string &src, string &sub, string sep)
{
	return GetNextString(src, sub, sep, 0);
}

/*!
 * 文字列が実数値を表しているかを調べる
 * @param[inout] str 文字列
 * @return 実数値ならtrue
 */
inline bool IsNumeric(const string &str)
{
	if(str.find_first_not_of("-+0123456789. Ee\t") != string::npos){
		return false;
	}
	return true;
}

/*!
 * テキストファイルから行列要素を読み込む
 * @param[in] file_name ファイル名
 * @param[in] sep 区切り文字
 * @param[out] mat 行列要素
 * @return ファイルが開けなかったらfalseを返す
 */
static bool ReadMatrix(string file_name, string sep, vector< vector<double> > &mat)
{
	ifstream file;
 
	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "ReadMatrix : Invalid file specified" << endl;
		return false;
	}
 
	string buf, grid;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);
 
		// '#'以降はコメントとして無視
		if( (comment_start = buf.find('#')) != string::size_type(-1) )
			buf = buf.substr(0, comment_start);
 
		// 空行は無視
		if(buf.empty())
			continue;

		vector<double> mat_line;
		size_t pos = 0;
		do{
			string sub;
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				mat_line.push_back(atof(sub.c_str()));
			}
		}while(pos != string::npos);

		if(!mat_line.empty()){
			mat.push_back(mat_line);
		}
	}
 
	file.close();
 
	return true;
}
 
/*!
 * テキストファイルから代数方程式の係数を読み込む
 * @param[in] file_name ファイル名
 * @param[in] sep 区切り文字
 * @param[out] c 係数列
 * @return ファイルが開けなかったらfalseを返す
 */
static bool ReadAlgebra(string file_name, string sep, vector<double> &c)
{
	ifstream file;
 
	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "ReadAlgebra : Invalid file specified" << endl;
		return false;
	}
 
	string buf, grid;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);
 
		// '#'以降はコメントとして無視
		if( (comment_start = buf.find('#')) != string::size_type(-1) )
			buf = buf.substr(0, comment_start);
 
		// 空行は無視
		if(buf.empty())
			continue;

		c.clear();
		size_t pos = 0;
		do{
			string sub;
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				c.push_back(atof(sub.c_str()));
			}
		}while(pos != string::npos);

		if(!c.empty()) break;
	}
 
	file.close();
 
	return true;
}

//-----------------------------------------------------------------------------
// 行列処理
//-----------------------------------------------------------------------------
/*!
 * 行列のファイル出力
 * @param[in] fp ファイルポインタ
 * @param[in] matrix 行列を格納した配列
 * @param[in] nx,ny  行列の大きさ
 */
inline void OutputMatrix(vector< vector<double> > matrix, int nx, int ny)
{
	for(int i = 0; i < nx; ++i){
		for(int j = 0; j < ny; ++j){
			cout << matrix[i][j] << (j == ny-1 ? "" : " ");
		}
		cout << endl;
	}
}



/*!
 * 2次元vectorコンテナに格納された正方行列同士の掛け算
 * @param[in] a,b nxn行列
 * @param[out] y  結果の行列(nxn)
 * @param[in] n 行列の大きさ
 * @return 1:成功
 */
inline int MulMatrix(const vector< vector<double> > &a, const vector< vector<double> > &b, vector< vector<double> > &y, int n)
{
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			y[i][j] = 0.0;
			for(int k = 0; k < n; ++k){
				y[i][j] += a[i][k]*b[k][j];
			}
		}
	}

	return 1;
}
inline vector< vector<double> > mul_mm(const vector< vector<double> > &a, const vector< vector<double> > &b, int n)
{
	vector< vector<double> > y(a);
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			y[i][j] = 0.0;
			for(int k = 0; k < n; ++k){
				y[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return y;
}

/*!
 * 2次元vectorコンテナに格納された行列の転置
 * @param[in] a nxn行列
 * @param[in] n 行列の大きさ
 * @return 転置行列
 */
inline vector< vector<double> > transpose(const vector< vector<double> > &a, int n)
{
	vector< vector<double> > t(n, vector<double>(n, 0.0));
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			t[i][j] = a[j][i];
		}
	}
	return t;
}

/*!
 * 2次元vectorコンテナに格納された行列とベクトルの掛け算
 * @param[in] a n×n行列
 * @param[in] b n次元ベクトル
 * @param[out] y 結果のベクトル(n)
 * @param[in] n 行列の大きさ
 * @return 1:成功
 */
inline int MulMatrixVector(const vector< vector<double> > &a, const vector<double> &b, vector<double> &y, int n)
{
	for(int i = 0; i < n; ++i){
		y[i] = 0.0;
		for(int k = 0; k < n; ++k){
			y[i] += a[i][k]*b[k];
		}
	}
	return 1;
}
inline vector<double> mul_mv(const vector< vector<double> > &a, const vector<double> &b, int n)
{
	vector<double> y(n);
	for(int i = 0; i < n; ++i){
		y[i] = 0.0;
		for(int k = 0; k < n; ++k){
			y[i] += a[i][k]*b[k];
		}
	}
	return y;
}



/*!
 * 1次元vectorコンテナに格納されたベクトル同士の内積
 * @param[in] a,b n次元ベクトル
 * @param[in] n ベクトルの大きさ
 * @return 内積
 */
inline double dot(const vector<double> &a, const vector<double> &b, int n)
{
	double d = 0.0;
	for(int i = 0; i < n; ++i){
		d += a[i]*b[i];
	}
	return d;
}

/*!
 * 1次元vectorコンテナに格納されたベクトルとスカラー値の掛け算
 * @param[in] a スカラー値
 * @param[in] b n次元ベクトル
 * @param[out] y 結果のベクトル(n)
 * @param[in] n 行列の大きさ
 * @return 1:成功
 */
inline vector<double> mul_sv(const double &a, const vector<double> &b, int n)
{
	vector<double> y(n);
	for(int i = 0; i < n; ++i){
		y[i] = a*b[i];
	}
	return y;
}

/*!
 * 1次元vectorコンテナに格納されたベクトルの正規化
 * @param[inout] a n次元ベクトル
 * @param[in] n ベクトルの大きさ
 * @return ノルムの値
 */
inline double normalize(vector<double> &a, int n)
{
	double l = sqrt(dot(a, a, n)); // |a|の計算
	if(fabs(l) < 1e-10) return 0.0;
	for(int i = 0; i < n; ++i){
		a[i] /= l;
	}
	return l;
}

inline void unit(vector< vector<double> > &a, int n)
{
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			a[i][j] = (i == j ? 1.0 : 0.0);
		}
	}
}




//-----------------------------------------------------------------------------
// サンプリング点(データ点)
//-----------------------------------------------------------------------------
/*!
 * サンプリング点(データ点)の生成(1次元)
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] dx サンプリング間隔
 * @param[in] func 関数値を与える関数ポインタ
 * @param[out] xi,yi サンプリングデータ
 * @return 生成されたデータ個数
 */
inline int MakeSamplingPoints(double x0, double x1, double dx, double func(double), vector<double> &xi, vector<double> &yi)
{
	xi.clear(); yi.clear();
	int cnt = 0;
	double x = x0;
	while(x <= x1){
		xi.push_back(x);
		yi.push_back(func(x));
		x += dx;
		cnt++;
	}
	return cnt;
}

/*!
 * サンプリング点(データ点)のファイル出力
 * @param[in] xi,yi サンプリングデータ
 * @param[in] filename 出力ファイル名
 */
inline void OutputSamplingPoints(vector<double> &xi, vector<double> &yi, string filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	for(int i = 0; i < xi.size(); ++i){
		fo << xi[i] << ", " << yi[i] << endl;
	}
	fo.close();
}

/*!
 * サンプリング点(データ点)の画面出力
 * @param[in] xi,yi サンプリングデータ
 * @param[in] filename 出力ファイル名
 */
inline void OutputSamplingPoints(vector<double> &xi, vector<double> &yi)
{
	cout << "sampling points : ";
	for(int i = 0; i < xi.size(); ++i){
		cout << "(" << xi[i] << ", " << yi[i] << ")" << (i == xi.size()-1 ? "" : ",  ");
	}
	cout << endl;
}

/*!
 * 関数値のファイル出力
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] dx サンプリング間隔
 * @param[in] func 関数値を与える関数ポインタ(std::bindをつかうためにstd::functionにしている)
 * @param[in] filename 出力ファイル名
 * @return 生成されたデータ個数
 */
inline int OutputFunction(double x0, double x1, double dx, std::function<double(double)> func, string filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	int cnt = 0;
	double x = x0;
	while(x <= x1){
		fo << x << ", ";
		fo << func(x) << endl;
		x += dx;
		cnt++;
	}
	fo.close();
	return cnt;
}

/*!
 * 関数値のファイル出力
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] dx サンプリング間隔
 * @param[in] func 関数値を与える関数ポインタ(std::bindをつかうためにstd::functionにしている)
 * @param[in] filename 出力ファイル名
 * @return 生成されたデータ個数
 */
inline int OutputFunction(double x0, double x1, double dx, vector<double> f, string filename)
{
	ofstream fo;
	fo.open(filename.c_str(), ios::out);
	double x = x0;
	for(int i = 0; i < f.size(); ++i){
		fo << x << "," << f[i] << endl;
		x += dx;
	}
	fo.close();
	return 0;
}

/*!
 * サンプリング点(データ点)の生成(1次元,ホワイトノイズ付き)
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] dx サンプリング間隔
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] nwidth ノイズのサイズ([-0.5nsize, 0.5nsize]のノイズを追加する)
 * @param[out] xi,yi サンプリングデータ
 * @return 生成されたデータ個数
 */
inline int MakeSamplingPointsWithWhiteNoise(double x0, double x1, double dx, double func(double), double nwidth, vector<double> &xi, vector<double> &yi)
{
	srand((unsigned int)time(0)); // 乱数のシード値を時間によって変える
	xi.clear(); yi.clear();
	int cnt = 0;
	double x = x0;
	while(x <= x1){
		double noise = ((double)rand()/(double)RAND_MAX-0.5)*nwidth;
		xi.push_back(x);
		yi.push_back(func(x)+noise);
		x += dx;
		cnt++;
	}
	return cnt;
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
inline int MakeChebyshevNodes(double x0, double x1, double dx, double func(double), vector<double> &xi, vector<double> &yi)
{
	xi.clear(); yi.clear();
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
// 偏微分方程式のための初期値設定とファイル出力
//-----------------------------------------------------------------------------
/*!
 * 関数の初期値として矩形波を設定
 * @param[out] f 未知関数fの各グリッドでの値
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
inline int SetValueRectangle(vector<double> &f, int n, double a, double b)
{
	if(n <= 0) return 1;
	double h = (b-a)/n; // 空間刻み幅
	double l = b-a; // 計算範囲全体の長さ
	double x = 0.0; // 計算半以内の相対位置(原点をaとする)
	for(int i = 0; i <= n; ++i){
		if(x >= 0.05*l && x <= 0.25*l){
			f[i] = 1.0;
		}
		else{
			f[i] = 0.0;
		}
		x += h;
	}
	return 0;
}

/*!
 * 関数の初期値として正弦波を設定
 * @param[out] f 未知関数fの各グリッドでの値
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 * @param[in] a,b 計算範囲
 * @return 問題なければ0を返す
 */
inline int SetValueSin(vector<double> &f, int n, double a, double b)
{
	if(n <= 0) return 1;
	double h = (b-a)/n; // 空間刻み幅
	double l = b-a; // 計算範囲全体の長さ
	double x = 0.0; // 計算半以内の相対位置(原点をaとする)
	for(int i = 0; i <= n; ++i){
		if(x >= 0.05*l && x <= 0.45*l){
			f[i] = 0.5*sin(RX_PI*(x-0.05*l)/(0.4*l));
		} else{
			f[i] = 0.0;
		}
		x += h;
	}
	return 0;
}

/*!
 * 関数値をファイル出力(1D)
 * @param[out] f 未知関数fの各グリッドでの値
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 * @param[in] a,b 計算範囲
 * @param[in] t 現在の時間(時間軸方向の現在地)
 * @return 問題なければ0を返す
 */
inline int OutputValueToFile(const vector<double> &f, int n, double a, double b, double t, ofstream &fo)
{
	if(n <= 0) return 1;

	// 結果の出力
	double x = a;		// 位置
	double h = (b-a)/n;	// 空間方向の刻み幅
	fo << t << ",";
	for(int i = 0; i <= n; ++i){
		fo << x << "," << f[i] << (i == n ? "" : ",");
		x += h;
	}
	fo << endl;
	return 0;
}

/*!
 * 関数値をファイル出力(2D)
 * @param[out] f 未知関数fの各グリッドでの値
 * @param[in] n 計算範囲内での分割数(h=(b-a)/n)
 * @param[in] a,b 計算範囲
 * @param[in] t 現在の時間(時間軸方向の現在地)
 * @return 問題なければ0を返す
 */
inline int OutputValueToFile(const vector< vector<double> > &f, int n, double x0, double xn, double y0, double yn, double t, ofstream &fo)
{
	if(n <= 0) return 1;
	fo << t << "," << n << "," << n << ",";

	double dx = (xn-x0)/n;	// 空間方向の刻み幅
	double dy = (yn-y0)/n;	// 空間方向の刻み幅

	fo << x0 << "," << y0 << "," << dx << "," << dy << ",";

	// 結果の出力
	for(int j = 0; j <= n; ++j){
		for(int i = 0; i <= n; ++i){
			fo << f[i] << ((j == n && i == n) ? "" : ",");
		}
	}
	fo << endl;
	return 0;
}


//-----------------------------------------------------------------------------
// 2次元座標データ
//-----------------------------------------------------------------------------
struct rxPoint2
{
	double x, y;
};

/*!
 * 2次元ベクトルデータの読み込み
 * @param[in] filename
 * @param[out] data
 * @return 正常に読み込めたら0を返す
 */
static int Read2d(const string &filename, vector<rxPoint2> &data, string sep = ",")
{
	ifstream file;
	file.open(filename.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "Read2d : Invalid file specified" << endl;
		return 1;
	}

	string buf, grid;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if((comment_start = buf.find('#')) != string::size_type(-1))
			buf = buf.substr(0, comment_start);

		// 空行は無視
		if(buf.empty())
			continue;

		rxPoint2 p;
		size_t pos = 0;
		string sub;

		// 座標値を1つずつ読み込んでいく
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			p.x = atof(sub.c_str());
		} else{
			break;
		}
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			p.y = atof(sub.c_str());
		} else{
			break;
		}

		data.push_back(p);
	}
	file.close();

	if(data.empty()) return 1;

	return 0;
}
/*!
 * 2次元ベクトルデータの読み込み
 * @param[in] filename
 * @param[out] data
 * @return 正常に読み込めたら0を返す
 */
static int Write2d(const string &filename, const vector<rxPoint2> &data, string header = "#p2d", string sep = ",")
{
	if(data.empty()) return 1;

	ofstream file;
	file.open(filename.c_str(), ios::out);
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "Write2d : Invalid file specified" << endl;
		return 1;
	}

	// 1行目はデータの種類を表す
	file << header << endl;

	// 点データを1座標/行で書き出す
	for(rxPoint2 p : data){
		file << p.x << sep << p.y << endl;
	}
	file.close();

	return 0;
}

/*!
 * 座標値データの範囲を求める
 * @param[in] data 2次元座標値データ
 * @param[out] min,max 最小,最大座標
 */
inline int SearchRange(vector<rxPoint2> &data, double min[2], double max[2])
{
	if(data.empty()) return 1;
	min[0] = max[0] = data[0].x;
	min[1] = max[1] = data[0].y;
	for(rxPoint2 p : data){
		if(p.x < min[0]) min[0] = p.x;
		if(p.x > max[0]) max[0] = p.x;
		if(p.y < min[1]) min[1] = p.y;
		if(p.y > max[1]) max[1] = p.y;
	}
	return 0;
}

#endif // #ifndef _RX_COMMON_H_