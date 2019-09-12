/*!
  @file graphgl.h

  @brief OpenGLによるグラフ描画

  @author Makoto Fujisawa
  @date 2019-06
*/

#ifndef _GRAPH_GL_H_
#define _GRAPH_GL_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// STL
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

// OpenGL
#include <GL/glew.h>
#include <GL/glut.h>

// 視点回転用
#include "rx_trackball.h"


using namespace std;

//-----------------------------------------------------------------------------
// プロトタイプ宣言とグローバル/スタティック変数
//-----------------------------------------------------------------------------
class GraphGL;
static GraphGL* g_GraphGLInstance = 0;

//-----------------------------------------------------------------------------
// OpenGLによるグラフ描画クラス
//  - freeglutでのウィンドウ生成も含む
//-----------------------------------------------------------------------------

class GraphGL
{
	struct Point
	{
		double x, y, z;
		double col[3];
		double size;
	};
	struct Line
	{
		Point p1, p2;
		double width;
	};
	struct TMesh
	{
		vector<Point> v;
		vector<int> t;
	};

	//! グリッドインデックスの計算
	inline int IDX(int i, int j){ return i+j*(m_iNx); }
	inline int IDX(int i, int j, int n){ return i+j*n; }

	int m_dim;
	double m_range[6];			//!< xmin,xmax,ymin,ymax,zmin,zmax
	double m_spacing[3];		//!< 目盛り間隔(x,y,z)
	double m_margin[3];			//!< 目盛数値と軸の間のマージン
	string m_label[3];			//!< 軸ラベル
	double m_label_margin[3];	//!< 軸ラベルと軸の間のマージン
	vector<Point> m_p;			//!< 点描画
	vector<Line> m_l;			//!< 線分描画

	vector<double> m_f;			//!< 関数値
	int m_draw;					//!< グラフセグメント描画方法(1:線,2:点,3(=1+2):点と線)

	TMesh m_mesh;				//!< 3Dグラフメッシュ
	int m_iNx, m_iNy;			//!< グリッド分割数
	double m_fDx, m_fDy;		//!< グリッド幅

	// ウィンドウ情報
	int m_iWinW;				//!< 描画ウィンドウの幅
	int m_iWinH;				//!< 描画ウィンドウの高さ
	int m_iWinX;				//!< 描画ウィンドウの位置x
	int m_iWinY;				//!< 描画ウィンドウの位置y

	bool m_bFullScreen = false;	//!< フルスクリーン表示
	bool m_bIdle = false;		//!< アイドル状態
	int m_iMouseButton = -1;	//!< マウスボタンの状態
	int m_iKeyMod = 0;			//!< 修飾キーの状態

	//! トラックボールによる視点移動
	rxTrackball m_tbView;

public:
	GraphGL(int w = 720, int h = 720, int x = 100, int y = 100, int dim = 2);
	~GraphGL(){}

public:
	// アニメーション,フルスクリーン切り替え
	void SwitchIdle(int on);
	void SwitchFullScreen(void);

	// GLUTイベントハンドラ
	void Display(void);
	void Resize(int w, int h);
	void Mouse(int button, int state, int x, int y);
	void Motion(int x, int y);
	void PassiveMotion(int x, int y);
	void Idle(void);
	void Keyboard(unsigned char key, int x, int y);
	void SpecialKey(int key, int x, int y);
	void MainMenu(int id);

	// イベントハンドラ設定用static関数
	static void DisplayCallback(){ g_GraphGLInstance->Display(); }
	static void ResizeCallback(int w, int h){ g_GraphGLInstance->Resize(w, h); }
	static void MouseCallback(int button, int state, int x, int y){ g_GraphGLInstance->Mouse(button, state, x, y); }
	static void MotionCallback(int x, int y){ g_GraphGLInstance->Motion(x, y); }
	static void PassiveMotionCallback(int x, int y){ g_GraphGLInstance->PassiveMotion(x, y); }
	static void IdleCallback(void){ g_GraphGLInstance->Idle(); }
	static void KeyboardCallback(unsigned char key, int x, int y){ g_GraphGLInstance->Keyboard(key, x, y); }
	static void SpecialKeyCallback(int key, int x, int y){ g_GraphGLInstance->SpecialKey(key, x, y); }
	static void MainMenuCallback(int id){ g_GraphGLInstance->MainMenu(id); }

protected:
	// GLの初期化関数
	void initGL(void);

	// GLUTの右クリックメニュー作成
	void initMenu(void);

	void drawString3D(string str, double * x);

	// グラフの各要素描画関数
	void drawFrame(void);
	void drawScale(void);
	void drawLabel(void);
	void drawPoints(void);
	void drawLines(void);

	void drawFrame3D(void);
	void drawScale3D(void);
	void drawLabel3D(void);
	void drawPoints3D(void);
	void drawLines3D(void);

	void drawGraphSegments(void);
	void drawGraphMesh(void);

	void generateMesh(double x1, double y1, double z1, double x2, double y2, double z2, int nx, int ny);

public:
	// ウィンドウ表示
	int Show(void);

	// グラフパラメータ設定
	void SetXRange(double xmin, double xmax){ m_range[0] = xmin; m_range[1] = xmax; }
	void SetYRange(double ymin, double ymax){ m_range[2] = ymin; m_range[3] = ymax; }
	void SetZRange(double zmin, double zmax){ m_range[4] = zmin; m_range[5] = zmax; }
	void SetXSpacing(double spacing){ m_spacing[0] = spacing; }
	void SetYSpacing(double spacing){ m_spacing[1] = spacing; }
	void SetZSpacing(double spacing){ m_spacing[2] = spacing; }
	void SetSpacing(double sx, double sy){ m_spacing[0] = sx; m_spacing[1] = sy; }
	void SetSpacing(double sx, double sy, double sz){ m_spacing[0] = sx; m_spacing[1] = sy; m_spacing[2] = sz; }
	void SetXMargin(double margin){ m_margin[0] = margin; }
	void SetYMargin(double margin){ m_margin[1] = margin; }
	void SetZMargin(double margin){ m_margin[2] = margin; }
	void SetMargin(double sx, double sy){ m_margin[0] = sx; m_margin[1] = sy; }
	void SetMargin(double sx, double sy, double sz){ m_margin[0] = sx; m_margin[1] = sy; m_margin[2] = sz; }
	void SetXLabel(const string &label){ m_label[0] = label; }
	void SetYLabel(const string &label){ m_label[1] = label; }
	void SetZLabel(const string &label){ m_label[2] = label; }
	void SetLabel(const string &sx, const string &sy){ m_label[0] = sx; m_label[1] = sy; }
	void SetLabel(const string &sx, const string &sy, const string &sz){ m_label[0] = sx; m_label[1] = sy; m_label[2] = sz; }
	void SetXLabelMargin(double margin){ m_label_margin[0] = margin; }
	void SetYLabelMargin(double margin){ m_label_margin[1] = margin; }
	void SetZLabelMargin(double margin){ m_label_margin[2] = margin; }
	void SetLabelMargin(double sx, double sy){ m_label_margin[0] = sx; m_label_margin[1] = sy; }
	void SetLabelMargin(double sx, double sy, double sz){ m_label_margin[0] = sx; m_label_margin[1] = sy; m_label_margin[2] = sz; }

	void SetGraph1D(double func(double), int nx, int drw = 1);
	void SetGraph2D(double func(vector<double>), int nx, int ny);

	// 描画用ポイント追加
	void AddPoint(double x, double y, double r = 1.0, double g = 1.0, double b = 1.0, double size = 1.0);
	void AddPoint3D(double x, double y, double z, double r = 1.0, double g = 1.0, double b = 1.0, double size = 1.0);

	// 描画用線分追加
	void AddLine(double x1, double y1, double x2, double y2, double r = 1.0, double g = 1.0, double b = 1.0, double width = 1.0);
	void AddLine3D(double x1, double y1, double z1, double x2, double y2, double z2, double r = 1.0, double g = 1.0, double b = 1.0, double width = 1.0);

};








#endif // #ifdef _GRAPH_GL_H_
