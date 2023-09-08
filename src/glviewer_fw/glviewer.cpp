/*! 
@file glviewer.cpp

@brief OpenGLによるシミュレーションデータのビューワ(GLFW&ImGUI版)

@author Makoto Fujisawa
@date 2023-05
*/


#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "glfw3dll.lib")

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#define GL_SILENCE_DEPRECATION	// mac環境でgluを使っている場合の非推奨warningの抑制


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_filelist.h"
#include "rx_bitmap.h"

// OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

// ImGUI
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"


using namespace std;

struct rxData1D
{
	vector<double> x;
	vector<double> f;
	int n;
	double t;
	double min[3], max[3];
	int SearchRange(void)
	{
		if(x.empty()) return 1;
		min[0] = max[0] = x[0];
		min[1] = max[1] = f[0];
		if(x.size() != f.size()){
			cout << "size of x and f is not matched" << endl;
			cout << "  x : " << x.size() << ", f : " << f.size() << endl;
			return 2;
		}
		for(int i = 1; i < x.size(); ++i){
			if(x[i] < min[0]) min[0] = x[i];
			if(x[i] > max[0]) max[0] = x[i];
			if(f[i] < min[1]) min[1] = f[i];
			if(f[i] > max[1]) max[1] = f[i];
		}
		min[2] = max[2] = 0;
		return 0;
	}
};

struct rxData2D
{
	//vector<double> x, y;
	vector<double> f;
	int nx, ny;
	double dx, dy, x0, y0;
	double t;
	double min[3], max[3];
	int SearchRange(void)
	{
		if(f.empty()) return 1;
		min[2] = max[2] = f[0];
		for(int i = 0; i < f.size(); ++i){
			if(f[i] < min[2]) min[2] = f[i];
			if(f[i] > max[2]) max[2] = f[i];
		}
		min[0] = max[0] = 0;
		min[1] = max[1] = 0;

		return 0;
	}

};

template<typename T> class rxTimeData
{
public:
	vector<T> data;
	int current_index;
	double min[3], max[3];
	double fscale, fmin;
	int state;
	int dim;

	rxTimeData(int d)
	{
		current_index = 0;
		for(int k = 0; k < 3; ++k){
			min[k] = 0.0;
			max[k] = 1.0;
		}
		fmin = 0.0;
		fscale = -1.0;
		state = 0;
		dim = d;
	}

	void SearchRange(void)
	{
		// データの範囲を調べる
		for(int k = 0; k < 3; ++k){
			min[k] = 1e10;
			max[k] = -1e10;
		}
		vector<T>::iterator itr = data.begin();
		for(; itr != data.end(); ++itr){
			T &d = *itr;
			d.SearchRange();
			for(int k = 0; k < 3; ++k){
				if(d.min[k] < min[k]) min[k] = d.min[k];
				if(d.max[k] > max[k]) max[k] = d.max[k];
			}
		}
	}

	bool Increment(int d = 1)
	{
		current_index += d;
		if(current_index >= data.size()){
			current_index = data.size()-1;
			return false;
		}
		return true;
	}
	bool Decrement(int d = 1)
	{
		current_index -= d;
		if(current_index < 0){
			current_index = 0;
			return false;
		}
		return true;
	}
};



//-----------------------------------------------------------------------------
// 定数・グローバル変数
//-----------------------------------------------------------------------------
// 光源位置/色 (shadow mapなどで使う場合用にグローバル変数にしている)
const glm::vec4 RX_LIGHT0_POS(0.0f, 0.0f, 1.0f, 1.0f);
const glm::vec4 RX_LIGHT1_POS(-1.0f, -10.0f, -1.0f, 0.0f);
const glm::vec4 RX_LIGHT_AMBI(0.3f, 0.3f, 0.3f, 1.0f);
const glm::vec4 RX_LIGHT_DIFF(1.0f, 1.0f, 1.0f, 1.0f);
const glm::vec4 RX_LIGHT_SPEC(1.0f, 1.0f, 1.0f, 1.0f);

// ウィンドウ/アニメーション/描画関連変数
int g_winw = 1200;							//!< 描画ウィンドウの幅
int g_winh = 600;							//!< 描画ウィンドウの高さ
bool g_animation_on = false;				//!< アニメーションON/OFF
int g_currentstep = 0;						//!< 現在のステップ数

//! データ格納用
rxTimeData<rxData1D> g_data1d(1);
rxTimeData<rxData2D> g_data2d(2);
int g_current_type = 1;

// ファイルリスト
vector<string> g_dfiles;
int g_dfile_idx = -1;

// 描画内容
int g_draw = 0;
double g_lc = 1.0;

// 画面更新間隔
float g_dt = 0.033f;
int g_graph_step = -1;
int g_step_inc = 1;
bool g_fast_step = true;

//-----------------------------------------------------------------------------
// シミュレーションデータの読み込み
//-----------------------------------------------------------------------------
/*!
* 1次元スカラーデータの読み込み
* @param[in] filename
* @param[out] data
* @return 正常に読み込めたら0を返す
*/
int Read1ds(const string &filename, vector<rxData1D> &data, string sep = ",")
{
	ifstream file;
	file.open(filename.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "Read1ds : Invalid file specified" << endl;
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

		rxData1D d;
		size_t pos = 0;

		// 最初の1つは時刻t
		string sub;
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.t = atof(sub.c_str());
		}
		else{
			break;
		}

		// 位置xとそこでの値fの読み込み
		do{
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				d.x.push_back(atof(sub.c_str()));
			}
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				d.f.push_back(atof(sub.c_str()));
			}
		} while(pos != string::npos);

		if(!d.x.empty()){
			d.n = d.x.size();
			data.push_back(d);
		}
	}
	file.close();

	if(data.empty()) return 1;

	return 0;
}

/*!
* 2次元スカラーデータの読み込み
* @param[in] filename
* @param[out] data
* @return 正常に読み込めたら0を返す
*/
int Read2ds(const string &filename, vector<rxData2D> &data, string sep = ",")
{
	ifstream file;
	file.open(filename.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "Read2ds : Invalid file specified" << endl;
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

		rxData2D d;
		size_t pos = 0;
		string sub;

		// 最初の1つは時刻t
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.t = atof(sub.c_str());
		} else{
			break;
		}

		// 次に続く2つはグリッド分割数nx,ny
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.nx = atof(sub.c_str());
		} else{
			break;
		}
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.ny = atof(sub.c_str());
		} else{
			break;
		}

		// 次の2つは左下座標x0,y0
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.x0 = atof(sub.c_str());
		} else{
			break;
		}
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.y0 = atof(sub.c_str());
		} else{
			break;
		}

		// 次の2つは空間刻み幅dx,dy
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.dx = atof(sub.c_str());
		} else{
			break;
		}
		pos = GetNextString(buf, sub, sep, pos);
		if(IsNumeric(sub)){
			d.dy = atof(sub.c_str());
		} else{
			break;
		}

		// 関数値fの読み込み
		do{
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				d.f.push_back(atof(sub.c_str()));
			}
		} while(pos != string::npos);

		if(!d.f.empty()){
			data.push_back(d);
		}
	}
	file.close();

	if(data.empty()) return 1;

	return 0;
}
/*!
* データファイルの読み込み
* @param[in] filename
* @return 正常に読み込めたら0を返す
*/
int Read(const string &filename)
{
	ifstream file;
	file.open(filename.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "Read : Invalid file specified" << endl;
		return 1;
	}

	// 1行目だけ読み込んでデータの種類を判別
	string buf;
	getline(file, buf);
	file.close();

	size_t pos = 0;
	string sub;
	pos = GetNextString(buf, sub, ",", pos);
	if(buf.find("#1d") != string::npos){ // 1次元スカラー値データ
		if(pos != string::npos){
			pos = GetNextString(buf, sub, ",", pos); // y軸最小値
			if(IsNumeric(sub)){
				g_data1d.fmin = atof(sub.c_str());
			}
			pos = GetNextString(buf, sub, ",", pos); // y軸スケール
			if(IsNumeric(sub)){
				g_data1d.fscale = atof(sub.c_str());
			}
		}
		g_data1d.data.clear();
		if(Read1ds(filename, g_data1d.data, ",")) return 1;
		g_data1d.state = 1;
		g_data1d.current_index = 0;
		g_data1d.SearchRange();
		g_current_type = 1;
	}
	else if(buf.find("#2d") != string::npos){ // 2次元スカラー値データ
		if(pos != string::npos){
			pos = GetNextString(buf, sub, ",", pos); // 関数値fの描画上の最小値
			if(IsNumeric(sub)){
				g_data2d.fmin = atof(sub.c_str());
			}
			pos = GetNextString(buf, sub, ",", pos); // 関数値fの描画上のスケール
			if(IsNumeric(sub)){
				g_data2d.fscale = atof(sub.c_str());
			}
		}
		g_data2d.data.clear();
		if(Read2ds(filename, g_data2d.data, ",")) return 1;
		g_data2d.state = 1;
		g_data2d.current_index = 0;
		g_data2d.SearchRange();
		g_current_type = 2;
	}

	return 0;
}


inline void ReadNext(int inc = 1)
{
	g_dfile_idx += inc;

	if(g_dfile_idx < 0) g_dfile_idx = 0;
	else if(g_dfile_idx > g_dfiles.size()-1) g_dfile_idx = g_dfiles.size()-1;
	else{ g_data1d.state = 0; g_data2d.state = 0; }

	Read(g_dfiles[g_dfile_idx]);
}

/*!
* 青->緑->赤->白と変化するサーモグラフ用の色生成
* @param[out] col 生成された色
* @param[in] x 値
* @param[in] xmin 最小値
* @param[in] xmax 最大値
*/
inline void CalThermograph(double col[3], double x, const double xmin = 0.0, const double xmax = 1.0)
{
	double l = xmax-xmin;
	if(fabs(l) < 1e-10) return;

	const int ncolors = 5;
	double base[ncolors][3] = { {0.0, 0.0, 1.0},
							   {0.0, 1.0, 1.0},
							   {0.0, 1.0, 0.0},
							   {1.0, 1.0, 0.0},
							   {1.0, 0.0, 0.0} };
	x = RX_CLAMP(((x-xmin)/l), 0.0, 1.0)*(ncolors-1);
	int i = (int)x;
	double dx = x-floor(x);
	col[0] = RX_LERP(base[i][0], base[i+1][0], dx);
	col[1] = RX_LERP(base[i][1], base[i+1][1], dx);
	col[2] = RX_LERP(base[i][2], base[i+1][2], dx);
}



/*!
* 2Dグラフの外枠描画
*/
void DrawFrame(double xmin, double xmax, double ymin, double ymax)
{
	glBegin(GL_LINE_LOOP);
	glVertex2d(xmin, ymin);
	glVertex2d(xmax, ymin);
	glVertex2d(xmax, ymax);
	glVertex2d(xmin, ymax);
	glEnd();
}


//-----------------------------------------------------------------------------
// アプリケーション制御関数
//-----------------------------------------------------------------------------
/*!
* アニメーションON/OFF
* @param[in] on trueでON, falseでOFF
*/
bool switchanimation(int on)
{
	g_animation_on = (on == -1) ? !g_animation_on : (on ? true : false);
	return g_animation_on;
}
/*!
* 現在の画面描画を画像ファイルとして保存(連番)
* @param[in] stp 現在のステップ数(ファイル名として使用)
*/
void savedisplay(void)
{
	string fn = "data_"+RX_TO_STRING(g_dfile_idx)+"_step"+RX_TO_STRING(g_current_type ==2 ? g_data2d.current_index : g_data1d.current_index)+".bmp";
	int c = 3;
	int wstep = (((g_winw+1)*c)/4)*4;
	vector<unsigned char> imm_buf(wstep*g_winh);
	glReadPixels(0, 0, g_winw, g_winh, GL_RGB, GL_UNSIGNED_BYTE, &imm_buf[0]);
	WriteBitmapFile(fn, &imm_buf[0], g_winw, g_winh, c, RX_BMP_WINDOWS_V3, wstep, false, true);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
* グラフ用ファイルの再読み込み
*/
void reload(void)
{
	if(g_dfile_idx >= 0) Read(g_dfiles[g_dfile_idx]);
}

/*!
* 初期化関数
*  - プログラム起動時に一回だけだけ実行したい処理はここに書く
*/
void Init(void)
{
	// OpenGLのバージョンチェック
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;

	// GLEWの初期化
	GLenum err = glewInit();
	if(err != GLEW_OK) cout << "GLEW Error : " << glewGetErrorString(err) << endl;

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClearDepth(1.0f);

	// 描画系フラグ設定(アンチエイリアス,デプステスト,隠面除去,法線計算,点描画)
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_POINT_SMOOTH);

	// 光源&材質描画設定
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, glm::value_ptr(RX_LIGHT0_POS));
	glLightfv(GL_LIGHT0, GL_DIFFUSE, glm::value_ptr(RX_LIGHT_DIFF));
	glLightfv(GL_LIGHT0, GL_SPECULAR, glm::value_ptr(RX_LIGHT_SPEC));
	glLightfv(GL_LIGHT0, GL_AMBIENT, glm::value_ptr(RX_LIGHT_AMBI));

	glShadeModel(GL_SMOOTH);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// データファイルリストの取得
	g_dfiles = GetFileList("./data/", "txt");
	for(vector<string>::iterator i = g_dfiles.begin(); i != g_dfiles.end(); ++i){
		std::cout << *i << std::endl;
	}
	if(g_dfiles.empty()){
		cout << "there is no data file!" << endl;
		g_dfile_idx = -1;
	}
	else{
		g_dfile_idx = 0;
	}

	// データ読み込み
	if(g_dfile_idx >= 0){
		Read(g_dfiles[g_dfile_idx]);
	}
	//Read("diffuse1d_eular.txt");
	//Read("advect1d_fe+upwind.txt");
	//Read("advect1d_heun+cen.txt");
}


//-----------------------------------------------------------------------------
// OpenGL/GLFWコールバック関数
//-----------------------------------------------------------------------------
/*!
* 再描画イベントコールバック関数
*/
void Display(void)
{
	// ビューポート,透視変換行列,モデルビュー変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glm::mat4 mp = glm::ortho(-1.0f, 1.0f, -1.0f, 1.0f);
	glMultMatrixf(glm::value_ptr(mp));
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 描画バッファのクリア
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glPushMatrix();

	g_graph_step = -1;

	// 1Dデータ描画(グラフ+サーモバー)
	if(g_current_type == 1 && g_data1d.state && !g_data1d.data.empty()){
		float xmin, xmax, ymin, ymax;
		xmin = g_data1d.min[0];
		xmax = g_data1d.max[0];
		if(g_data1d.fscale > 0.0){	// データファイルにy軸のスケールが設定されていたらそれを用いる
			ymin = g_data1d.fmin;
			ymax = g_data1d.fmin+g_data1d.fscale;
		}
		else{	// 設定されていなかったらデータの最大，最小値を使う
			ymin = g_data1d.min[1];
			ymax = g_data1d.max[1];
		}

		// 関数値データ
		g_graph_step = g_data1d.current_index;
		rxData1D &d = g_data1d.data[g_graph_step];

		float xmargin = 0.05*(xmax-xmin);	// x方向余白
		float ymargin = 0.2*(ymax-ymin);	// y方向余白
		float yoffset = 0.5*(ymax-ymin);	// 下部分に別のものを描画するためのオフセット
		// 透視変換行列の設定
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glm::mat4 mp = glm::ortho(xmin-xmargin, xmax+xmargin, ymin-yoffset-ymargin, ymax+ymargin);
		glMultMatrixf(glm::value_ptr(mp));

		// モデルビュー変換行列の設定
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		// 軸と目盛り
		double lc = fabs(g_lc-0.5);
		glColor3d(lc, lc, lc);
		glLineWidth(2.0);
		//DrawFrame(xmin, xmax, ymin, ymax);
		glBegin(GL_LINES);
		glVertex2d(xmin, ymin);	glVertex2d(xmax, ymin);
		glVertex2d(xmin, ymin);	glVertex2d(xmin, ymax);
		glVertex2d(xmax, ymin);	glVertex2d(xmax, ymax);
		double ym = 10.0, xm = 0.0;
		while(ym < ymax){
			glVertex2d(xmin, ym); glVertex2d(xmin+0.5, ym);
			ym += 10.0;
		}
		while(xm < xmax){
			glVertex2d(xm, ymin); glVertex2d(xm, ymin+1.0);
			xm += 10.0;
		}
		glEnd();

		//// 目盛数値
		//vector<string> num(1);
		//num[0] = RX_TO_STRING(xmin);
		//glColor3d(1.0, 1.0, 1.0);
		//float sx0 = g_winw*xmargin/(xmax-xmin)-36, sy0 = g_winh*(ymargin+yoffset)/(ymax-ymin)-36; // 原点位置(ウィンドウピクセル座標)
		//DrawStrings(num, g_winw, g_winh, sx0, sy0);
		//ym = 10.0; xm = 0.0;
		//while(ym < ymax){
		//	num[0] = RX_TO_STRING(ym);
		//	DrawStrings(num, g_winw, g_winh, sx0, sy0-g_winh*ym/(ymax-ymin+2*ymargin+yoffset));
		//	ym += 10.0;
		//}
		//while(xm < xmax){
		//	num[0] = RX_TO_STRING(xm);
		//	DrawStrings(num, g_winw, g_winh, sx0+g_winw*xm/(xmax-xmin+2*xmargin)+20, sy0+15);
		//	xm += 10.0;
		//}


		// グラフ
		glColor3d(g_lc, g_lc, g_lc);
		glLineWidth(3.0);
		glBegin(GL_LINE_STRIP);
		for(int i = 0; i < d.x.size(); ++i){
			double x = d.x[i];
			double f = d.f[i];
			glVertex2d(x, f);
		}
		glEnd();

		// サーモバー
		double y0 = ymin-yoffset;
		double y1 = ymin-ymargin;
		glBegin(GL_QUADS);
		for(int i = 0; i < d.x.size()-1; ++i){
			double x0 = d.x[i];
			double f0 = d.f[i];
			double x1 = d.x[i+1];
			double f1 = d.f[i+1];

			double c[3];
			CalThermograph(c, f0, ymin, ymax);
			glColor3dv(c);
			glVertex2d(x0, y1);
			glVertex2d(x0, y0);
			CalThermograph(c, f1, ymin, ymax);
			glColor3dv(c);
			glVertex2d(x1, y0);
			glVertex2d(x1, y1);
		}
		glEnd();
	}
	else if(g_current_type == 2 && g_data2d.state && !g_data2d.data.empty()){
		float dx, dy;
		dx = g_data2d.data[0].dx;
		dy = g_data2d.data[0].dy;
		int nx, ny;
		nx = g_data2d.data[0].nx;
		ny = g_data2d.data[0].ny;
		float xmin, xmax, ymin, ymax, fmin, fmax;
		xmin = g_data2d.data[0].x0;
		xmax = xmin+nx*dx;
		ymin = g_data2d.data[0].y0;
		ymax = ymin+ny*dy;
		if(g_data2d.fscale > 0.0){	// データファイルにy軸のスケールが設定されていたらそれを用いる
			fmin = g_data2d.fmin;
			fmax = g_data2d.fmin+g_data2d.fscale;
		} else{	// 設定されていなかったらデータの最大，最小値を使う
			fmin = g_data2d.min[2];
			fmax = g_data2d.max[2];
		}

		// 関数値データ
		g_graph_step = g_data2d.current_index;
		rxData2D &d = g_data2d.data[g_graph_step];

		if(g_draw == 1)
		{
			// 透視変換行列の設定
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glm::mat4 mp = glm::perspective(45.0f, (float)g_winw/g_winh, 0.05f, 100.0f);
			glMultMatrixf(glm::value_ptr(mp));

			// モデルビュー変換行列の設定
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
		}
		else{
			float xmargin = 0.05*(xmax-xmin);	// x方向余白
			float ymargin = 0.05*(ymax-ymin);	// y方向余白

			// 透視変換行列の設定
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glm::mat4 mp = glm::ortho(xmin-xmargin, xmax+xmargin, ymin-ymargin, ymax+ymargin);
			glMultMatrixf(glm::value_ptr(mp));

			// モデルビュー変換行列の設定
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			// 2Dサーモ表示
			double y0 = ymin;
			double y1 = ymin-2*ymargin;
			glBegin(GL_QUADS);
			for(int j = 0; j < d.ny; ++j){
				for(int i = 0; i < d.nx; ++i){
					int k0 = i+j*(d.nx+1);
					int k1 = (i+1)+(j+1)*(d.nx+1);
					double x0 = xmin+i*dx;
					double y0 = xmin+j*dy;
					double x1 = x0+dx;
					double y1 = y0+dy;

					double f00 = d.f[k0];
					double f11 = d.f[k1];
					double f01 = d.f[(i)+(j+1)*(d.nx+1)];
					double f10 = d.f[(i+1)+(j)*(d.nx+1)];

					double c[3];
					CalThermograph(c, f00, fmin, fmax);
					glColor3dv(c); glVertex2d(x0, y0);
					CalThermograph(c, f10, fmin, fmax);
					glColor3dv(c); glVertex2d(x1, y0);
					CalThermograph(c, f11, fmin, fmax);
					glColor3dv(c); glVertex2d(x1, y1);
					CalThermograph(c, f01, fmin, fmax);
					glColor3dv(c); glVertex2d(x0, y1);

				}
			}
			glEnd();
		}
	}

	glPopMatrix();
}


/*!
 * タイマーイベント処理関数(ある時間間隔で実行)
 */
void Timer(void)
{
	if(g_animation_on){
		g_step_inc = g_fast_step ? 10 : 1; 
		if(g_current_type == 2) g_animation_on = g_data2d.Increment(g_step_inc);
		else g_animation_on = g_data1d.Increment(g_step_inc);

		g_currentstep++;
	}
}

/*!
* ヘルプテキストを表示
*/
void help(void)
{
	static const char* help = "<ESC> key - quit the program\n"
		"'j','k' key - increase/decrease 1 step\n"
		"Shift+'j','k' key - increase/decrease 10 steps\n"
		"'n','m' key - select next/previous file\n"
		"'r' key - read selected file\n"
		"'l' key - show file list\n"
		"Shift+'f' key - fullscreen on/off\n"
		"'h' key - show this help\n";
	std::cout << help << std::endl;
}

/*!
* キーボードイベント処理関数
* @param[in] window コールバック関数を呼んだウィンドウハンドル
* @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
* @param[in] scancode キーのスキャンコード(プラットフォーム依存)
* @param[in] action アクション(GLFW_PRESS:キーを押す, GLFW_RELESE:キーを離す，GLFW_REPEAT:キーリピート機能時)
* @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
*/
void Keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if(ImGui::GetIO().WantCaptureKeyboard) return;	// ImGUIウィンドウ上でのキーボードイベント時
	g_step_inc = (mods & GLFW_MOD_SHIFT) ? 10 : 1;
	if(action == GLFW_PRESS || action == GLFW_REPEAT){
		switch(key){
		case GLFW_KEY_ESCAPE:	// ESC,Qキーでアプリケーション終了
		case GLFW_KEY_Q:
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;

		case GLFW_KEY_S: // SキーでアニメーションON/OFF
			switchanimation(-1);
			break;
		case GLFW_KEY_SPACE: // スペースキーでアニメーションを1ステップだけ進める
			g_animation_on = true; Timer(); g_animation_on = false;
			break;

		// 背景色変更
		case GLFW_KEY_W: 
			break;
			g_lc = 1-g_lc;
			glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0f);
			break;
		case GLFW_KEY_B:
			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			g_lc = 1.0;
			break;

		// データのステップを進める/戻す
		case GLFW_KEY_J:
			if(g_current_type == 2) g_data2d.Increment(g_step_inc);
			else g_data1d.Increment(g_step_inc);
			break;
		case GLFW_KEY_K:
			if(g_current_type == 2) g_data2d.Decrement(g_step_inc);
			else g_data1d.Decrement(g_step_inc);
			break;

		// リスト内の次/前のファイルを選択
		case GLFW_KEY_N:
			ReadNext(1);
			break;
		case GLFW_KEY_M:
			ReadNext(-1);
			break;

		// データファイル再読み込み
		case GLFW_KEY_R:
			reload();
			break;

		// データファイルリスト表示
		case GLFW_KEY_L:
			for(vector<string>::iterator i = g_dfiles.begin(); i != g_dfiles.end(); ++i){
				std::cout << *i << std::endl;
			}
			break;

		//case GLFW_KEY_O: // 画像ファイル出力
		//	SaveFrameBuffer("data_"+RX_TO_STRING(g_dfile_idx)+"_step"+RX_TO_STRING(g_current_type ==2 ? g_data2d.current_index : g_data1d.current_index)+".bmp", g_winw, g_winh);
		//	break;


		default:
			break;
		}
	}
}


/*!
* リサイズイベント処理関数
* @param[in] window コールバック関数を呼んだウィンドウハンドル
* @param[in] w キャンバス幅(ピクセル数)
* @param[in] h キャンバス高さ(ピクセル数)
*/
void Resize(GLFWwindow* window, int w, int h)
{
	g_winw = w; g_winh = h;
	glViewport(0, 0, g_winw, g_winh);
}

/*!
* ImGUIのウィジット配置
*  - ImGUI/imgui_demo.cppを参考に ( https://github.com/ocornut/imgui#demo )
* @param[in] window コールバック関数を呼んだウィンドウハンドル
*/
void SetImGUI(GLFWwindow* window)
{
	// Start the ImGui frame
	ImGui_ImplOpenGL2_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// GUI
	ImGui::Begin("ImGui Window");

	// テキスト表示
	ImGui::Text("data:");
	if(g_dfile_idx >= 0){
		ImGui::Text(g_dfiles[g_dfile_idx].c_str());
		if(ImGui::Button("next")) ReadNext(1);
		ImGui::SameLine();
		if(ImGui::Button("prev")) ReadNext(-1);
	}
	else{
		ImGui::Text("no data");
	}
	ImGui::Separator();

	ImGui::Text("animation:");
	if(g_graph_step >= 0) ImGui::Text(("step "+RX_TO_STRING(g_graph_step)).c_str());
	if(ImGui::Button("start/stop")){ switchanimation(-1); } ImGui::SameLine();
	if(ImGui::Button("run a step")){ g_animation_on = true; Timer(); g_animation_on = false; }
	if(ImGui::Checkbox("fast", &g_fast_step)){ g_step_inc = g_fast_step ? 10 : 1; }
	if(ImGui::Button("reload")) reload();
	ImGui::Separator();
	if(ImGui::Button("save screenshot")){ savedisplay(); }
	if(ImGui::Button("quit")){ glfwSetWindowShouldClose(window, GL_TRUE); }

	ImGui::End();

}

void Clean()
{
}


/*!
* GLFWのエラーコールバック関数
* @param[in] error エラーコード
* @param[in] description エラー内容
*/
static inline void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}


/*!
* メインルーチン
* @param[in] argc コマンドライン引数の数
* @param[in] argv コマンドライン引数
*/
int main(int argc, char *argv[])
{
	if(!glfwInit()) return 1;
	glfwSetErrorCallback(glfw_error_callback);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
#ifdef __APPLE__
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

	// Create window
	GLFWwindow* window = glfwCreateWindow(g_winw, g_winh, "OpenGL Application", NULL, NULL);
	if(window == NULL) return 1;

	// Set glfw window as current OpenGL rendering context
	glfwMakeContextCurrent(window);
	glewExperimental = GL_TRUE;
	glfwSwapInterval(0); // Disable vsync

	// Initilization
	Init();

	help();

	// Setup callback functions
	glfwSetKeyCallback(window, Keyboard);
	glfwSetFramebufferSizeCallback(window, Resize);
	Resize(window, g_winw, g_winh);

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsClassic();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL2_Init();

	// Settings for timer to display FPS on ImGUI window
	float cur_time = 0.0f, last_time = 0.0f, elapsed_time = 0.0f;
	glfwSetTime(0.0);	// Initialize the glfw timer

	// Main loop
	while(!glfwWindowShouldClose(window))
	{
		// Poll and handle events (inputs, window resize, etc.)
		glfwPollEvents();

		// OpenGL Rendering & Animation function
		Display();

		// Timer
		cur_time = glfwGetTime();
		elapsed_time = cur_time-last_time;
		if(elapsed_time >= g_dt){
			Timer();
			last_time = glfwGetTime();
		}

		// GUI
		SetImGUI(window);

		// Rendering of the ImGUI frame in opengl canvas
		ImGui::Render();
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	// Cleanup
	Clean();
	ImGui_ImplOpenGL2_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();


	return 0;
}


