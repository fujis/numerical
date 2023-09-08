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

// QR法
#include "qr.h"

using namespace std;



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
int g_winw = 1000;							//!< 描画ウィンドウの幅
int g_winh = 1000;							//!< 描画ウィンドウの高さ
bool g_animation_on = false;				//!< アニメーションON/OFF
int g_currentstep = 0;						//!< 現在のステップ数
float g_dt = 0.033f;

//! データ格納用
vector<rxPoint2> g_data;
double g_min[2] = { 0, 0 }, g_max[2] = { 1, 1 };

// ファイルリスト
vector<string> g_dfiles;
int g_dfile_idx = -1;

// 描画内容
int g_draw = 0;
double g_lc = 1.0;
double g_margin = 0.05;

// 固有値
vector<double> g_lambda; // 固有値
vector< vector<double> > g_v; // 固有ベクトル
rxPoint2 g_c;

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
	string fn = "data_"+RX_TO_STRING(g_dfile_idx)+".bmp";
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
	g_dfiles = GetFileList("./points", "txt");
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


//-----------------------------------------------------------------------------
// OpenGL/GLFWコールバック関数
//-----------------------------------------------------------------------------
/*!
* 再描画イベントコールバック関数
*/
void Display(void)
{
	float xmargin = g_margin*(g_max[0]-g_min[0]);	// x方向余白
	float ymargin = g_margin*(g_max[1]-g_min[1]);	// y方向余白

	// ビューポート,透視変換行列,モデルビュー変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glm::mat4 mp = glm::ortho<float>(g_min[0]-xmargin, g_max[0]+xmargin, g_min[1]-ymargin, g_max[1]+ymargin);
	glMultMatrixf(glm::value_ptr(mp));
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 描画バッファのクリア
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glPushMatrix();

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
}


/*!
 * タイマーイベント処理関数(ある時間間隔で実行)
 */
void Timer(void)
{
	if(g_animation_on){
		g_currentstep++;
	}
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
* @param[in] window コールバック関数を呼んだウィンドウハンドル
* @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
* @param[in] scancode キーのスキャンコード(プラットフォーム依存)
* @param[in] action アクション(GLFW_PRESS:キーを押す, GLFW_RELESE:キーを離す，GLFW_REPEAT:キーリピート機能時)
* @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
*/
void Keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if(ImGui::GetIO().WantCaptureKeyboard) return;	// ImGUIウィンドウ上でのキーボードイベント時
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
		case GLFW_KEY_B: 
			break;
			g_lc = 1-g_lc;
			glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0f);
			break;

		// データのステップを進める/戻す
		case GLFW_KEY_J:
			g_dfile_idx++;
			if(g_dfile_idx > g_dfiles.size()-1) g_dfile_idx = g_dfiles.size()-1;
			break;
		case GLFW_KEY_K:
			g_dfile_idx--;
			if(g_dfile_idx < 0) g_dfile_idx = 0;
			break;

		// データファイル再読み込み
		case GLFW_KEY_R:
			if(g_dfile_idx >= 0) Read(g_dfiles[g_dfile_idx]);
			g_lambda.clear();
			g_v.clear();
			break;

		// データファイルリスト表示
		case GLFW_KEY_L:
			for(vector<string>::iterator i = g_dfiles.begin(); i != g_dfiles.end(); ++i){
				std::cout << *i << std::endl;
			}
			break;

		case GLFW_KEY_C:	// データクリア
			g_data.clear();
			g_lambda.clear();
			g_v.clear();
			break;

		case GLFW_KEY_W: // データファイル出力
			Write("./points/points.txt");
			break;

		case GLFW_KEY_E: // 固有値/固有ベクトル計算
			CalEigen();
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
* マウスイベント処理関数
* @param[in] window コールバック関数を呼んだウィンドウハンドル
* @param[in] button マウスボタン(GLFW_MOUSE_BUTTON_LEFT,GLFW_MOUSE_BUTTON_MIDDLE,GLFW_MOUSE_BUTTON_RIGHT)
* @param[in] action マウスボタンの状態(GLFW_PRESS, GLFW_RELEASE)
* @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
*/
void Mouse(GLFWwindow* window, int button, int action, int mods)
{
	if(ImGui::GetIO().WantCaptureMouse) return;	// ImGUIウィンドウ上でのマウスイベント時
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	if(button == GLFW_MOUSE_BUTTON_LEFT){
		if(action == GLFW_PRESS){
			double xmargin = g_margin*(g_max[0]-g_min[0]);	// x方向余白
			double ymargin = g_margin*(g_max[1]-g_min[1]);	// y方向余白
			double lx = g_max[0]-g_min[0]+2*xmargin;
			double ly = g_max[1]-g_min[1]+2*ymargin;

			rxPoint2 p;
			p.x = (x/(double)g_winw)*lx+(g_min[0]-xmargin);
			p.y = ((g_winh-y)/(double)g_winh)*ly+(g_min[1]-ymargin);
			g_data.push_back(p);
		}
		else if(action == GLFW_RELEASE){
		}
	}
}
/*!
* モーションイベント処理関数(マウスボタンを押したままドラッグ)
* @param[in] window コールバック関数を呼んだウィンドウハンドル
* @param[in] x,y マウス座標(スクリーン座標系)
*/
void Motion(GLFWwindow* window, double x, double y)
{
	if(ImGui::GetIO().WantCaptureMouse) return;	// ImGUIウィンドウ上でのマウスイベント時
	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
	   glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE){
		return;
	}

	if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS){
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
		if(ImGui::Button("next")){
			g_dfile_idx++;
			if(g_dfile_idx > g_dfiles.size()-1) g_dfile_idx = g_dfiles.size()-1;
		}
		ImGui::SameLine();
		if(ImGui::Button("prev")){
			g_dfile_idx--;
			if(g_dfile_idx < 0) g_dfile_idx = 0;
		}
	}
	else{
		ImGui::Text("no data");
	}
	if(ImGui::Button("reload")) reload();
	ImGui::Separator();

	ImGui::Text("animation:");
	if(ImGui::Button("start/stop")){ switchanimation(-1); } ImGui::SameLine();
	if(ImGui::Button("run a step")){ g_animation_on = true; Timer(); g_animation_on = false; }
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
	glfwSetCursorPosCallback(window, Motion);
	glfwSetMouseButtonCallback(window, Mouse);
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


