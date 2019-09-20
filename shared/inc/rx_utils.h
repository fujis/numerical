/*! @file rx_common.h
	
	@brief 数値計算テストの共通ヘッダ
 
	@author Makoto Fujisawa
	@date   2012
*/


#ifndef _RX_COMMON_H_
#define _RX_COMMON_H_

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
inline int GetNextString(const string &src, string &sub, string sep, int pos)
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
		return (int)string::npos;
	}
	else{
		int cnt = 1;
		while(src[i+cnt] == ' '){	// sepの後のスペースを消す
			cnt++;
		}
		if(!extracted) sub = src.substr(pos, i-pos);
		return (int)(i+cnt >= src.size() ? (int)string::npos : i+cnt);
	}
}
 
/*!
 * stringから最初の区切り文字までを抽出
 * @param[in] src 元の文字列
 * @param[out] sub 抽出文字列
 * @param[in] sep 区切り文字
 * @return 次の抽出開始位置(sepの後にスペースがあればそのスペースの後)
 */
inline int GetFirstString(const string &src, string &sub, string sep)
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
	if(str.find_first_not_of("-0123456789. Ee\t") != string::npos){
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
bool ReadMatrix(string file_name, string sep, vector< vector<double> > &mat)
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
bool ReadAlgebra(string file_name, string sep, vector<double> &c)
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
static void OutputMatrix(vector< vector<double> > matrix, int nx, int ny)
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

/*!
 * 2次元vectorコンテナに格納された行列の転置
 * @param[in] a nxn行列
 * @param[in] n 行列の大きさ
 * @return 転置行列
 */
inline vector< vector<double> > Transpose(const vector< vector<double> > &a, int n)
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

/*!
 * ベクトルの内積
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


#endif // #ifndef _RX_COMMON_H_