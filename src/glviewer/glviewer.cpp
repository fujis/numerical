/*! 
  @file glviewer.cpp
	
  @brief OpenGLによるシミュレーションデータのビューワ
 
  @author Makoto Fujisawa
  @date 2019-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_filelist.h"
#include "rx_bitmap.h"

// OpenGL
#include <GL/glut.h>

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

	void Increment(int d = 1)
	{
		current_index += d;
		if(current_index >= data.size()) current_index = data.size()-1;
	}
	void Decrement(int d = 1)
	{
		current_index -= d;
		if(current_index < 0) current_index = 0;
	}
};


//-----------------------------------------------------------------------------
// グローバル変数
//-----------------------------------------------------------------------------
//! データ格納用
rxTimeData<rxData1D> g_data1d(1);
rxTimeData<rxData2D> g_data2d(2);
int g_current_type = 1;

// ファイルリスト
vector<string> g_dfiles;
int g_dfile_idx = -1;

// 描画領域サイズ(ピクセル数)
int g_w = 1200, g_h = 600;

// 描画内容
int g_draw = 0;
double g_lc = 1.0;

//-----------------------------------------------------------------------------
// 関数プロトタイプ宣言
//-----------------------------------------------------------------------------
void SwitchIdle(int on = -1);
void SwitchFullScreen(void);
void CleanGL(void);
bool SaveFrameBuffer(const string &fn, int w, int h);

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
 * 文字列描画
 * @param[in] static_str 文字列
 * @param[in] w,h ウィンドウサイズ
 */
static void DrawStrings(vector<string> &static_str, int w, int h)
{
	glDisable(GL_LIGHTING);
	// 平行投影にする
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	float x0 = 5.0f;
	float y0 = h-20.0f;

	// 画面上部にテキスト描画
	for(int j = 0; j < (int)static_str.size(); ++j){
		glRasterPos2f(x0, y0);

		int size = (int)static_str[j].size();
		for(int i = 0; i < size; ++i){
			char ic = static_str[j][i];
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ic);
		}

		y0 -= 20;
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
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
// OpenGLイベントハンドラ
//-----------------------------------------------------------------------------
/*!
 * 再描画イベント処理関数
 */
void Display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	int cur_step = -1;

	// 1Dデータ描画(グラフ+サーモバー)
	if(g_current_type == 1 && g_data1d.state && !g_data1d.data.empty()){
		double xmin, xmax, ymin, ymax;
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

		double lc = fabs(g_lc-0.5);
		glColor3d(lc, lc, lc);
		glLineWidth(2.0);
		//DrawFrame(xmin, xmax, ymin, ymax);
		glBegin(GL_LINES);
		glVertex2d(xmin, ymin);
		glVertex2d(xmax, ymin);
		glVertex2d(xmin, ymin);
		glVertex2d(xmin, ymax);
		glVertex2d(xmax, ymin);
		glVertex2d(xmax, ymax);
		glEnd();


		// 関数値データ
		cur_step = g_data1d.current_index;
		rxData1D &d = g_data1d.data[cur_step];

		double xmargin = 0.05*(xmax-xmin);	// x方向余白
		double ymargin = 0.2*(ymax-ymin);	// y方向余白
		double yoffset = 0.5*(ymax-ymin);	// 下部分に別のものを描画するためのオフセット
		// 透視変換行列の設定
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(xmin-xmargin, xmax+xmargin, ymin-yoffset-ymargin, ymax+ymargin);

		// モデルビュー変換行列の設定
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

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
		double dx, dy;
		dx = g_data2d.data[0].dx;
		dy = g_data2d.data[0].dy;
		int nx, ny;
		nx = g_data2d.data[0].nx;
		ny = g_data2d.data[0].ny;
		double xmin, xmax, ymin, ymax, fmin, fmax;
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
		cur_step = g_data2d.current_index;
		rxData2D &d = g_data2d.data[cur_step];

		if(g_draw == 1)
		{
			// 透視変換行列の設定
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(45.0f, (float)g_w/(float)g_h, 0.2f, 1000.0f);

			// モデルビュー変換行列の設定
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
		}
		else{
			double xmargin = 0.05*(xmax-xmin);	// x方向余白
			double ymargin = 0.05*(ymax-ymin);	// y方向余白

			// 透視変換行列の設定
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluOrtho2D(xmin-xmargin, xmax+xmargin, ymin-ymargin, ymax+ymargin);

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

	// テキスト表示
	vector<string> strs;
	if(g_dfile_idx >= 0) strs.push_back(g_dfiles[g_dfile_idx]);
	if(cur_step >= 0) strs.push_back("step "+RX_TO_STRING(cur_step));
	glColor3d(1.0, 1.0, 1.0);
	DrawStrings(strs, g_w, g_h);

	glutSwapBuffers();
}


/*!
 * リサイズイベント処理関数
 * @param[in] w キャンバス幅(ピクセル数)
 * @param[in] h キャンバス高さ(ピクセル数)
 */
void Resize(int w, int h)
{
	g_w = w; g_h = h;
	glViewport(0, 0, w, h);

	// 透視変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(45.0f, (float)w/(float)h, 0.2f, 1000.0f);
	//glOrtho(-1, 1, -1, 1, -1, 1);
	gluOrtho2D(-1, 1, -1, 1);

	// モデルビュー変換行列の設定
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

/*!
 * マウスイベント処理関数
 * @param[in] button マウスボタン(GLUT_LEFT_BUTTON,GLUT_MIDDLE_BUTTON,GLUT_RIGHT_BUTTON)
 * @param[in] state マウスボタンの状態(GLUT_UP, GLUT_DOWN)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void Mouse(int button, int state, int x, int y)
{
	if(x < 0 || y < 0) return;

	int mod = glutGetModifiers();	// SHIFT,CTRL,ALTの状態取得

	if(button == GLUT_LEFT_BUTTON){
		if(state == GLUT_DOWN){
		} else if(state == GLUT_UP){
		}
	} else if(button == GLUT_MIDDLE_BUTTON){
	} else if(button == GLUT_RIGHT_BUTTON){
	}

	glutPostRedisplay();
}

/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void Motion(int x, int y)
{
	if(x < 0 || y < 0) return;
	glutPostRedisplay();
}


/*!
 * アイドルイベント処理関数
 */
void Idle(void)
{
	glutPostRedisplay();
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
 * @param[in] key キーの種類
 * @param[in] x,y キーが押されたときのマウス座標(スクリーン座標系)
 */
void Keyboard(unsigned char key, int x, int y)
{
	int mod = glutGetModifiers();	// SHIFT,CTRL,ALTの状態取得
	int inc = 1;

	switch(key){
	case '\033':  // '\033' は ESC の ASCII コード
		CleanGL(); exit(1);
		break;

	case 'h':	// ヘルプテキスト表示
		help();
		break;

	case ' ':	// アニメーション1ステップだけ実行
		Idle();
		break;

	case 's':	// アニメーションON/OFF
		SwitchIdle(-1);
		break;

	case 'F':	// フルスクリーン
		SwitchFullScreen();
		break;

	// 背景色変更
	case 'w':
		g_lc = 1-g_lc;
		glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0f);
		break;
	case 'b':
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		g_lc = 1.0;
		break;

	// データのステップを進める/戻す
	case 'J':
		inc = 10;
	case 'j':
		if(g_current_type == 2) g_data2d.Increment(inc);
		else g_data1d.Increment(inc);
		break;
	case 'K':
		inc = 10;
	case 'k':
		if(g_current_type == 2) g_data2d.Decrement(inc);
		else g_data1d.Decrement(inc);
		break;

	// リスト内の次/前のファイルを選択
	case 'n':
		g_dfile_idx++;
		if(g_dfile_idx > g_dfiles.size()-1) g_dfile_idx = g_dfiles.size()-1;
		else{ g_data1d.state = 0; g_data2d.state = 0; }
		break;
	case 'm':
		g_dfile_idx--;
		if(g_dfile_idx < 0) g_dfile_idx = 0;
		else{ g_data1d.state = 0; g_data2d.state = 0; }
		break;

	// データファイル読み込み
	case 'r':
		if(g_dfile_idx >= 0) Read(g_dfiles[g_dfile_idx]);
		break;

	// データファイルリスト表示
	case 'l':
		for(vector<string>::iterator i = g_dfiles.begin(); i != g_dfiles.end(); ++i){
			std::cout << *i << std::endl;
		}
		break;

	case 'o': // 画像ファイル出力
		SaveFrameBuffer("data_"+RX_TO_STRING(g_dfile_idx)+"_step"+RX_TO_STRING(g_current_type ==2 ? g_data2d.current_index : g_data1d.current_index)+".bmp", g_w, g_h);
		break;

	default:
	break;
	}

	glutPostRedisplay();
}

/*!
 * GLの初期化関数
 */
void InitGL(void)
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClearDepth(1.0f);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	// ポリゴン法線設定
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

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

/*!
 * GLの終了関数
 */
void CleanGL(void)
{
}


/*!
 * メインルーチン
 * @param[in] argc コマンドライン引数の数
 * @param[in] argv コマンドライン引数
 */
int main(int argc, char *argv[])
{
	help();

	glutInitWindowPosition(100, 100);
	glutInitWindowSize(g_w, g_h);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ACCUM);
	glutCreateWindow(argv[0]);

	glutDisplayFunc(Display);
	glutReshapeFunc(Resize);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);
	glutKeyboardFunc(Keyboard);

	SwitchIdle(0);

	InitGL();

	glutMainLoop();

	CleanGL();


	return 0;
}




//-----------------------------------------------------------------------------
// アプリケーション制御
//-----------------------------------------------------------------------------

/*!
 * アイドル関数のON/OFF
 * @param[in] on trueでON, falseでOFF
 */
void SwitchIdle(int on)
{
	static bool idle_state = false;
	idle_state = (on == -1) ? !idle_state : (on ? true : false);
	glutIdleFunc((idle_state ? Idle : 0));
	cout << "idle " << (idle_state ? "on" : "off") << endl;
}


/*!
 * フルスクリーン/ウィンドウ表示の切り替え
 */
void SwitchFullScreen(void)
{
	static int fullscreen = 0;		// フルスクリーン状態
	static int pos0[2] = { 0, 0 };
	static int win0[2] = { 500, 500 };
	if(fullscreen){
		glutPositionWindow(pos0[0], pos0[1]);
		glutReshapeWindow(win0[0], win0[1]);
	} else{
		pos0[0] = glutGet(GLUT_WINDOW_X);
		pos0[1] = glutGet(GLUT_WINDOW_Y);
		win0[0] = glutGet(GLUT_WINDOW_WIDTH);
		win0[1] = glutGet(GLUT_WINDOW_HEIGHT);

		glutFullScreen();
	}
	fullscreen ^= 1;
}



/*!
 * 描画の画像保存
 */
bool SaveFrameBuffer(const string &fn, int w, int h)
{
	int c = 3;
	vector<unsigned char> img_buf(w*h*c);

	int format = GL_RGBA;
	if(c == 3){
		format = GL_RGB;
	}
	//glFlush();
	glutSwapBuffers();

	glRasterPos2i(0, h);
	glReadPixels(0, 0, w, h, format, GL_UNSIGNED_BYTE, &img_buf[0]);

	WriteBitmapFile(fn, &img_buf[0], w, h, c, RX_BMP_WINDOWS_V3, -1, false, true);

	//for(int i = 0; i < w; ++i){
	//	for(int j = 0; j < h; ++j){
	//		int idx = 3*(w*(h-j-1)+i);
	//		//int idx = 3*(j*w+i);

	//		int r, g, b;
	//		r = img_buf[idx+0];
	//		g = img_buf[idx+1];
	//		b = img_buf[idx+2];

	//		SetPixel(cvimg, i, j, r, g, b);
	//	}
	//}

	return true;
}
