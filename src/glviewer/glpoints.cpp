/*! 
  @file glpoints.cpp
	
  @brief OpenGLによる2次元サンプル点作成/描画
 
  @author Makoto Fujisawa
  @date 2019-12
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_filelist.h"
#include "rx_bitmap.h"

// QR法
#include "qr.h"

// OpenGL
#include <GL/glut.h>

using namespace std;



//-----------------------------------------------------------------------------
// グローバル変数
//-----------------------------------------------------------------------------
//! データ格納用
vector<rxPoint2> g_data;
double g_min[2] = { 0, 0 }, g_max[2] = { 1, 1 };

// ファイルリスト
vector<string> g_dfiles;
int g_dfile_idx = -1;

// 描画領域サイズ(ピクセル数)
int g_w = 800, g_h = 800;

// 描画内容
int g_draw = 0; // 0で入力モード,1で表示モード
double g_lc = 1.0;
double g_margin = 0.05;

// 固有値
vector<double> g_lambda; // 固有値
vector< vector<double> > g_v; // 固有ベクトル
rxPoint2 g_c;

//-----------------------------------------------------------------------------
// 関数プロトタイプ宣言
//-----------------------------------------------------------------------------
void SwitchIdle(int on = -1);
void SwitchFullScreen(void);
void CleanGL(void);
bool SaveFrameBuffer(const string &fn, int w, int h);




//-----------------------------------------------------------------------------
// データの読み込みと固有値計算
//-----------------------------------------------------------------------------
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
	if(buf.find("#p2d") != string::npos){ // 2次元ベクトル値データ
		g_data.clear();
		if(Read2d(filename, g_data, ",")) return 1;

		if(pos != string::npos){
			pos = GetNextString(buf, sub, ",", pos); // データ範囲最小値
			if(IsNumeric(sub)){
				g_min[0] = g_min[1] = atof(sub.c_str());
			}
			pos = GetNextString(buf, sub, ",", pos); // データ範囲最大値
			if(IsNumeric(sub)){
				g_max[0] = g_max[1] = atof(sub.c_str());
			}
		}
		else{
			SearchRange(g_data, g_min, g_max);
		}

	}

	return 0;
}

/*!
 * データファイルの書き出し
 * @param[in] filename
 * @return 正常に読み込めたら0を返す
 */
int Write(const string &filename)
{
	return Write2d(filename, g_data, "#p2d,0,1", ",");
}

/*!
 * 点群から固有値・固有ベクトルを計算
 */
void CalEigen(void)
{
	if(g_data.empty()) return;

	g_lambda.clear();
	g_v.clear();

	// 共分散行列の計算
	int n = 2;
	vector< vector<double> > A, A0;
	covariance_matrix2(g_data, A, g_c);
	A0 = A;

	// 固有値
	int max_iter = 100;
	double eps = 1e-6;
	qr(A, n, max_iter, eps);
	g_lambda.resize(n, 0.0);
	for(int i = 0; i < n; ++i){
		g_lambda[i] = A[i][i];
	}

	// 固有値の表示
	cout << "e = (";
	for(int i = 0; i < n; ++i){
		cout << g_lambda[i] << (i == n-1 ? ")" : ", ");
	}
	cout << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	// 逆反復法で固有ベクトルを計算
	max_iter = 100;
	eps = 1e-6;
	g_v.clear();
	inverse_iteration(A0, g_lambda, g_v, n, max_iter, eps);

	// 固有ベクトルの表示
	for(int i = 0; i < n; ++i){
		cout << "v" << i+1 << " = (";
		for(int j = 0; j < n; ++j){
			cout << g_v[i][j] << (j == n-1 ? ")" : ", ");
		}
		cout << endl;
	}
	cout << "iter = " << max_iter << ", eps = " << eps << endl;

}

//-----------------------------------------------------------------------------
// 描画関数
//-----------------------------------------------------------------------------
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
 * 矢印の描画
 * @param[in] s 矢印の始点
 * @param[in] e 矢印の始点
 * @param[in] scale 傘の部分の大きさ(全体の長さに対する係数)
 */
void DrawArrow2D(const double *s, const double *d, double scale = 0.2)
{
	// 始点，方向，長さ
	double origin[2], dir[2];
	origin[0] = s[0]; origin[1] = s[1];
	dir[0] = d[0]; dir[1] = d[1];
	
	double length = sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
	if(length > 1e-6){ dir[0] /= length; dir[1] /= length; }

	// ベクトル(1,0)との間の角度
	double theta = 180.0/RX_PI*acos(dir[0]*(dir[1] > 0 ? 1 : -1));

	// 矢印の傘部分の設定
	double arrow_x = scale*length;                // 軸方向の長さ
	//double arrow_y = arrow_x*0.174532925;        // 軸に垂直な方向の開き量(傘の開き角度20deg)
	double arrow_y = arrow_x*0.363970234;        // 軸に垂直な方向の開き量(傘の開き角度40deg)

	glPushMatrix();

	glTranslatef(origin[0], origin[1], 0.0);    // 矢印原点に移動
	glRotatef(theta, 0.0, 0.0, 1.0);            // 矢印方向に回転(z軸中心)

	glBegin(GL_LINES);
	// 軸
	glVertex2d(0.0, 0.0);
	glVertex2d(length, 0.0);

	// 傘
	glVertex2d(length, 0.0);
	glVertex2d(length-arrow_x, arrow_y);
	glVertex2d(length, 0.0);
	glVertex2d(length-arrow_x, -arrow_y);
	glEnd();

	glPopMatrix();
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

	double xmargin = g_margin*(g_max[0]-g_min[0]);	// x方向余白
	double ymargin = g_margin*(g_max[1]-g_min[1]);	// y方向余白

	// 透視変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(g_min[0]-xmargin, g_max[0]+xmargin, g_min[1]-ymargin, g_max[1]+ymargin);

	// モデルビュー変換行列の設定
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 周りのフレーム
	double lc = fabs(g_lc-0.5);
	glColor3d(lc, lc, lc);
	glLineWidth(2.0);
	DrawFrame(g_min[0], g_max[0], g_min[1], g_max[1]);

	// データ描画
	if(!g_data.empty()){
		glColor3d(g_lc, g_lc, g_lc);
		glPointSize(10.0);
		glBegin(GL_POINTS);
		for(rxPoint2 p : g_data){
			glVertex2d(p.x, p.y);
		}
		glEnd();
	}

	// 固有値/固有ベクトル描画
	if(!g_lambda.empty()){
		double cen[2];
		cen[0] = g_c.x; cen[1] = g_c.y;

		// 重心
		glColor3d(1.0, 1.0, 0.0);
		glPointSize(10.0);
		glBegin(GL_POINTS);
		glVertex2dv(cen);
		glEnd();
		glColor3d(g_lc, g_lc, g_lc);
		glPointSize(20.0);
		glBegin(GL_POINTS);
		glVertex2dv(cen);
		glEnd();

		// 固有ベクトル
		double ld, dir[2];
		ld = 2*sqrt(fabs(g_lambda[0])); // 分散は2乗されているので平方根をとって長さの目安にする(2倍は見た目を合わせるためのマジックナンバー)
		dir[0] = ld*g_v[0][0]; dir[1] = ld*g_v[0][1];
		glColor3d(0.0, 0.0, 1.0);
		glLineWidth(3.0);
		DrawArrow2D(cen, dir, 0.1);
		ld = 2*sqrt(fabs(g_lambda[1]));
		dir[0] = ld*g_v[1][0]; dir[1] = ld*g_v[1][1];
		glColor3d(1.0, 0.0, 0.0);
		glLineWidth(3.0);
		DrawArrow2D(cen, dir, 0.1);
	}

	glPopMatrix();

	// テキスト表示
	vector<string> strs;
	if(g_dfile_idx >= 0) strs.push_back(g_dfiles[g_dfile_idx]);
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
			double xmargin = g_margin*(g_max[0]-g_min[0]);	// x方向余白
			double ymargin = g_margin*(g_max[1]-g_min[1]);	// y方向余白
			double lx = g_max[0]-g_min[0]+2*xmargin;
			double ly = g_max[1]-g_min[1]+2*ymargin;

			rxPoint2 p;
			p.x = (x/(double)g_w)*lx+(g_min[0]-xmargin);
			p.y = ((g_h-y)/(double)g_h)*ly+(g_min[1]-ymargin);
			g_data.push_back(p);
		}
		else if(state == GLUT_UP){
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
		"'j','k' key - select next/previous file\n"
		"'r' key - read selected file\n"
		"'l' key - show file list\n"
		"'b' key - switch background color\n"
		"'e' key - calculate eigen value/vector\n"
		"'c' key - clear points\n"
		"'w' key - output points as text file\n"
		"'o' key - output as image file\n"
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

	case 'b':	// 背景色切り替え
		g_lc = 1-g_lc;
		glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0f);
		break;

	// リスト内の次/前のファイルを選択
	case 'j':
		g_dfile_idx++;
		if(g_dfile_idx > g_dfiles.size()-1) g_dfile_idx = g_dfiles.size()-1;
		break;
	case 'k':
		g_dfile_idx--;
		if(g_dfile_idx < 0) g_dfile_idx = 0;
		break;

	case 'r':	// データファイル読み込み
		if(g_dfile_idx >= 0) Read(g_dfiles[g_dfile_idx]);
		g_lambda.clear();
		g_v.clear();
		break;

	case 'l':	// データファイルリスト表示
		for(vector<string>::iterator i = g_dfiles.begin(); i != g_dfiles.end(); ++i){
			std::cout << *i << std::endl;
		}
		break;

	case 'c':	// データクリア
		g_data.clear();
		g_lambda.clear();
		g_v.clear();
		break;

	case 'w': // データファイル出力
		Write("./points/points.txt");
		break;

	case 'e': // 固有値/固有ベクトル計算
		CalEigen();
		break;

	case 'o': // 画像ファイル出力
		SaveFrameBuffer("points.bmp", g_w, g_h);
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
	//glEnable(GL_MULTISAMPLE);

	// 点描画設定
	glEnable(GL_POINT_SMOOTH);

	// データファイルリストの取得
	g_dfiles = GetFileList("./points/", "txt");
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
	//if(g_dfile_idx >= 0){
	//	Read(g_dfiles[g_dfile_idx]);
	//}
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
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ACCUM | GLUT_MULTISAMPLE);
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
