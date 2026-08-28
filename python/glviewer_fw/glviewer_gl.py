# -*- coding: utf-8 -*-
"""!
@file glviewer_gl.py

@brief OpenGLによるシミュレーションデータのビューワ(PyOpenGL+GLFW+簡易GUI版)
       C++版の glviewer.cpp (OpenGL+GLFW+ImGUI) をほぼそのまま移植したもの
       ImGUIの代わりに simple_gui.py の簡易GUIを使っている

       使い方: python glviewer_gl.py [データフォルダ]
        - フォルダを省略すると ../pde/data を読みに行く

       必要なパッケージ: pip install PyOpenGL glfw pillow

@author Makoto Fujisawa
@date 2023-05 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

import glfw
from OpenGL.GL import *

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
from viewer_data import *
from simple_gui import SimpleGUI


# -----------------------------------------------------------------------------
# 定数・グローバル変数
# -----------------------------------------------------------------------------
g_winw = 1200                # 描画ウィンドウの幅
g_winh = 600                 # 描画ウィンドウの高さ
g_animation_on = False       # アニメーションON/OFF
g_currentstep = 0            # 現在のステップ数

g_data1d = rxTimeData(1)     # データ格納用
g_data2d = rxTimeData(2)
g_current_type = 1

g_dfiles = []                # ファイルリスト
g_dfile_idx = -1

g_lc = 1.0                   # 線の色(1.0で白)

g_dt = 0.033                 # 画面更新間隔
g_graph_step = -1
g_step_inc = 1
g_fast_step = True

g_gui = None                 # 簡易GUI


# -----------------------------------------------------------------------------
# データの読み込み
# -----------------------------------------------------------------------------
def ReadFile(idx):
    """!
    ファイルリストのidx番目のファイルを読み込む
    @param[in] idx ファイル番号
    """
    global g_dfile_idx, g_current_type
    if not g_dfiles:
        return
    if idx < 0:
        idx = 0
    if idx > len(g_dfiles)-1:
        idx = len(g_dfiles)-1
    g_dfile_idx = idx
    g_data1d.state = 0
    g_data2d.state = 0
    ret, t = Read(g_dfiles[g_dfile_idx], g_data1d, g_data2d)
    if ret == 0:
        g_current_type = t
        n = len(g_data2d.data) if t == 2 else len(g_data1d.data)
        print("read : " + g_dfiles[g_dfile_idx] + "  (" + str(n) + " steps)")


def ReadNext(inc=1):
    """! リスト内の次/前のファイルを読み込む """
    ReadFile(g_dfile_idx+inc)


def current_time_data():
    """! 現在表示中の種類の時系列データを返す """
    return g_data2d if g_current_type == 2 else g_data1d


# -----------------------------------------------------------------------------
# 描画関数
# -----------------------------------------------------------------------------
def DrawFrame(xmin, xmax, ymin, ymax):
    """! 2Dグラフの外枠描画 """
    glBegin(GL_LINE_LOOP)
    glVertex2d(xmin, ymin)
    glVertex2d(xmax, ymin)
    glVertex2d(xmax, ymax)
    glVertex2d(xmin, ymax)
    glEnd()


def ortho2d(left, right, bottom, top):
    """! 平行投影(正射影)の設定 """
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(left, right, bottom, top, -1.0, 1.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()


# -----------------------------------------------------------------------------
# アプリケーション制御関数
# -----------------------------------------------------------------------------
def switchanimation(on=-1):
    """!
    アニメーションON/OFF
    @param[in] on Trueに相当する1でON, 0でOFF, -1でトグル
    """
    global g_animation_on
    g_animation_on = (not g_animation_on) if on == -1 else bool(on)
    return g_animation_on


def savedisplay():
    """! 現在の画面描画を画像ファイルとして保存 """
    from PIL import Image
    d = current_time_data()
    fn = "data_" + str(g_dfile_idx) + "_step" + str(d.current_index) + ".png"
    w, h = glfw.get_framebuffer_size(glfw.get_current_context())
    glPixelStorei(GL_PACK_ALIGNMENT, 1)
    buf = glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE)
    img = Image.frombytes("RGB", (w, h), buf).transpose(Image.FLIP_TOP_BOTTOM)
    img.save(fn)
    print("saved the screen image to " + fn)


def reload_file():
    """! グラフ用ファイルの再読み込み """
    if g_dfile_idx >= 0:
        ReadFile(g_dfile_idx)


def Init(folder):
    """!
    初期化関数
     - プログラム起動時に一回だけ実行したい処理はここに書く
    @param[in] folder データフォルダ
    """
    global g_dfiles, g_dfile_idx
    # OpenGLのバージョンチェック
    print("OpenGL version: " + glGetString(GL_VERSION).decode())
    print("Vendor: " + glGetString(GL_VENDOR).decode())
    print("Renderer: " + glGetString(GL_RENDERER).decode())

    glClearColor(0.0, 0.0, 0.0, 1.0)
    glClearDepth(1.0)

    # 描画系フラグ設定(アンチエイリアス,デプステスト,隠面除去,点描画)
    glEnable(GL_MULTISAMPLE)
    glEnable(GL_DEPTH_TEST)
    glDisable(GL_CULL_FACE)
    glEnable(GL_POINT_SMOOTH)
    glShadeModel(GL_SMOOTH)

    # データファイルリストの取得
    g_dfiles = GetFileList(folder, "txt")
    for f in g_dfiles:
        print(f)
    if not g_dfiles:
        print("there is no data file!")
        g_dfile_idx = -1
    else:
        # データ読み込み
        ReadFile(0)


# -----------------------------------------------------------------------------
# OpenGL/GLFWコールバック関数
# -----------------------------------------------------------------------------
def Display():
    """!
    再描画イベントコールバック関数
    """
    global g_graph_step

    # ビューポート,透視変換行列,モデルビュー変換行列の設定
    ortho2d(-1.0, 1.0, -1.0, 1.0)

    # 描画バッファのクリア
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glDisable(GL_LIGHTING)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

    g_graph_step = -1

    # 1Dデータ描画(グラフ+サーモバー)
    if g_current_type == 1 and g_data1d.state and g_data1d.data:
        Display1D()
    elif g_current_type == 2 and g_data2d.state and g_data2d.data:
        Display2D()


def Display1D():
    """! 1Dデータ描画(グラフ+サーモバー) """
    global g_graph_step

    xmin = g_data1d.min[0]
    xmax = g_data1d.max[0]
    if g_data1d.fscale > 0.0:  # データファイルにy軸のスケールが設定されていたらそれを用いる
        ymin = g_data1d.fmin
        ymax = g_data1d.fmin+g_data1d.fscale
    else:                      # 設定されていなかったらデータの最大，最小値を使う
        ymin = g_data1d.min[1]
        ymax = g_data1d.max[1]

    # 関数値データ
    g_graph_step = g_data1d.current_index
    d = g_data1d.data[g_graph_step]

    xmargin = 0.05*(xmax-xmin)   # x方向余白
    ymargin = 0.2*(ymax-ymin)    # y方向余白
    yoffset = 0.5*(ymax-ymin)    # 下部分に別のものを描画するためのオフセット

    # 透視変換行列とモデルビュー変換行列の設定
    ortho2d(xmin-xmargin, xmax+xmargin, ymin-yoffset-ymargin, ymax+ymargin)

    # 軸と目盛り
    lc = abs(g_lc-0.5)
    glColor3d(lc, lc, lc)
    glLineWidth(2.0)
    glBegin(GL_LINES)
    glVertex2d(xmin, ymin)
    glVertex2d(xmax, ymin)
    glVertex2d(xmin, ymin)
    glVertex2d(xmin, ymax)
    glVertex2d(xmax, ymin)
    glVertex2d(xmax, ymax)
    ym, xm = 10.0, 0.0
    while ym < ymax:
        glVertex2d(xmin, ym)
        glVertex2d(xmin+0.5, ym)
        ym += 10.0
    while xm < xmax:
        glVertex2d(xm, ymin)
        glVertex2d(xm, ymin+1.0)
        xm += 10.0
    glEnd()

    # グラフ
    glColor3d(g_lc, g_lc, g_lc)
    glLineWidth(3.0)
    glBegin(GL_LINE_STRIP)
    for i in range(len(d.x)):
        glVertex2d(d.x[i], d.f[i])
    glEnd()

    # サーモバー
    y0 = ymin-yoffset
    y1 = ymin-ymargin
    glBegin(GL_QUADS)
    for i in range(len(d.x)-1):
        x0 = d.x[i]
        f0 = d.f[i]
        x1 = d.x[i+1]
        f1 = d.f[i+1]

        c = CalThermograph(f0, ymin, ymax)
        glColor3d(*c)
        glVertex2d(x0, y1)
        glVertex2d(x0, y0)
        c = CalThermograph(f1, ymin, ymax)
        glColor3d(*c)
        glVertex2d(x1, y0)
        glVertex2d(x1, y1)
    glEnd()


def Display2D():
    """! 2Dデータ描画(サーモ表示) """
    global g_graph_step

    d0 = g_data2d.data[0]
    dx, dy = d0.dx, d0.dy
    nx, ny = d0.nx, d0.ny
    xmin = d0.x0
    xmax = xmin+nx*dx
    ymin = d0.y0
    ymax = ymin+ny*dy
    if g_data2d.fscale > 0.0:  # データファイルにスケールが設定されていたらそれを用いる
        fmin = g_data2d.fmin
        fmax = g_data2d.fmin+g_data2d.fscale
    else:                      # 設定されていなかったらデータの最大，最小値を使う
        fmin = g_data2d.min[2]
        fmax = g_data2d.max[2]

    # 関数値データ
    g_graph_step = g_data2d.current_index
    d = g_data2d.data[g_graph_step]

    xmargin = 0.05*(xmax-xmin)  # x方向余白
    ymargin = 0.05*(ymax-ymin)  # y方向余白

    # 透視変換行列とモデルビュー変換行列の設定
    ortho2d(xmin-xmargin, xmax+xmargin, ymin-ymargin, ymax+ymargin)

    # 2Dサーモ表示
    glBegin(GL_QUADS)
    for j in range(d.ny):
        for i in range(d.nx):
            x0 = xmin+i*dx
            y0 = ymin+j*dy
            x1 = x0+dx
            y1 = y0+dy

            f00 = d.f[(i)+(j)*(d.nx+1)]
            f10 = d.f[(i+1)+(j)*(d.nx+1)]
            f11 = d.f[(i+1)+(j+1)*(d.nx+1)]
            f01 = d.f[(i)+(j+1)*(d.nx+1)]

            glColor3d(*CalThermograph(f00, fmin, fmax))
            glVertex2d(x0, y0)
            glColor3d(*CalThermograph(f10, fmin, fmax))
            glVertex2d(x1, y0)
            glColor3d(*CalThermograph(f11, fmin, fmax))
            glVertex2d(x1, y1)
            glColor3d(*CalThermograph(f01, fmin, fmax))
            glVertex2d(x0, y1)
    glEnd()


def Timer():
    """!
    タイマーイベント処理関数(ある時間間隔で実行)
    """
    global g_animation_on, g_step_inc, g_currentstep
    if g_animation_on:
        g_step_inc = 10 if g_fast_step else 1
        g_animation_on = current_time_data().Increment(g_step_inc)
        g_currentstep += 1


def help_text():
    """! ヘルプテキストを表示 """
    print("<ESC>,'q' key - quit the program\n"
          "'s' key - start/stop animation\n"
          "<SPACE> key - run a step\n"
          "'j','k' key - increase/decrease 1 step\n"
          "Shift+'j','k' key - increase/decrease 10 steps\n"
          "'n','m' key - select next/previous file\n"
          "'r' key - read selected file\n"
          "'l' key - show file list\n"
          "'o' key - save screenshot\n"
          "'h' key - show this help")


def Keyboard(window, key, scancode, action, mods):
    """!
    キーボードイベント処理関数
    @param[in] window コールバック関数を呼んだウィンドウハンドル
    @param[in] key キーの種類
    @param[in] scancode キーのスキャンコード(プラットフォーム依存)
    @param[in] action アクション(GLFW_PRESS/GLFW_RELEASE/GLFW_REPEAT)
    @param[in] mods 修飾キー(CTRL,SHIFT,ALT)
    """
    global g_step_inc, g_animation_on, g_lc
    inc = 10 if (mods & glfw.MOD_SHIFT) else 1
    if action in (glfw.PRESS, glfw.REPEAT):
        d = current_time_data()
        if key in (glfw.KEY_ESCAPE, glfw.KEY_Q):   # ESC,Qキーでアプリケーション終了
            glfw.set_window_should_close(window, True)
        elif key == glfw.KEY_S:                    # SキーでアニメーションON/OFF
            switchanimation(-1)
        elif key == glfw.KEY_SPACE:                # スペースキーで1ステップだけ進める
            g_animation_on = True
            Timer()
            g_animation_on = False
        elif key == glfw.KEY_B:                    # 背景色変更
            g_lc = 1-g_lc
            glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0)
        elif key == glfw.KEY_J:                    # データのステップを進める
            d.Increment(inc)
        elif key == glfw.KEY_K:                    # データのステップを戻す
            d.Decrement(inc)
        elif key == glfw.KEY_N:                    # リスト内の次のファイルを選択
            ReadNext(1)
        elif key == glfw.KEY_M:                    # リスト内の前のファイルを選択
            ReadNext(-1)
        elif key == glfw.KEY_R:                    # データファイル再読み込み
            reload_file()
        elif key == glfw.KEY_L:                    # データファイルリスト表示
            for f in g_dfiles:
                print(f)
        elif key == glfw.KEY_O:                    # 画像ファイル出力
            savedisplay()
        elif key == glfw.KEY_H:                    # ヘルプ表示
            help_text()


def Resize(window, w, h):
    """!
    リサイズイベント処理関数
    @param[in] window コールバック関数を呼んだウィンドウハンドル
    @param[in] w,h キャンバスの幅,高さ(ピクセル数)
    """
    global g_winw, g_winh
    g_winw, g_winh = w, h
    glViewport(0, 0, g_winw, g_winh)


def SetGUI(window):
    """!
    簡易GUIのウィジット配置
     - C++版のImGUIウィンドウに相当
    @param[in] window コールバック関数を呼んだウィンドウハンドル
    """
    global g_fast_step, g_step_inc, g_animation_on

    gui = g_gui
    gui.begin(window, x=10, y=10, width=200, title="Viewer")
    gui.draw_panel_background()

    # テキスト表示
    gui.text("data:")
    if g_dfile_idx >= 0:
        gui.text(os.path.basename(g_dfiles[g_dfile_idx]))
        if gui.button("next", w=90):
            ReadNext(1)
        gui.same_line()
        if gui.button("prev", w=90):
            ReadNext(-1)
    else:
        gui.text("no data")
    gui.separator()

    gui.text("animation:")
    d = current_time_data()
    if g_graph_step >= 0:
        gui.text("step " + str(g_graph_step) + " / " + str(len(d.data)-1))
    if gui.button("start/stop"):
        switchanimation(-1)
    if gui.button("run a step"):
        g_animation_on = True
        Timer()
        g_animation_on = False
    fast = gui.checkbox("fast", g_fast_step)
    if fast != g_fast_step:
        g_fast_step = fast
        g_step_inc = 10 if g_fast_step else 1
    if d.data:
        v = gui.slider("step", d.current_index, 0, len(d.data)-1)
        d.current_index = int(RX_CLAMP(int(v), 0, len(d.data)-1))
    if gui.button("reload"):
        reload_file()
    gui.separator()
    if gui.button("save screenshot"):
        savedisplay()
    if gui.button("quit"):
        glfw.set_window_should_close(window, True)

    gui.end()


# -----------------------------------------------------------------------------
# メインルーチン
# -----------------------------------------------------------------------------
def main():
    global g_gui

    # データフォルダはコマンドライン引数で指定できる(省略時は ../pde/data )
    folder = sys.argv[1] if len(sys.argv) > 1 else default_data_folder()
    print("data folder : " + folder)

    if not glfw.init():
        return 1

    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 2)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 1)

    # ウィンドウ生成
    window = glfw.create_window(g_winw, g_winh, "glviewer (PyOpenGL)", None, None)
    if not window:
        glfw.terminate()
        return 1

    glfw.make_context_current(window)
    glfw.swap_interval(1)

    # 初期化
    Init(folder)
    g_gui = SimpleGUI()

    help_text()

    # コールバック関数の設定
    glfw.set_key_callback(window, Keyboard)
    glfw.set_framebuffer_size_callback(window, Resize)
    Resize(window, g_winw, g_winh)

    # メインループ
    last_time = glfw.get_time()
    while not glfw.window_should_close(window):
        # イベント処理
        glfw.poll_events()

        # 描画
        Display()

        # タイマー
        cur_time = glfw.get_time()
        if cur_time-last_time >= g_dt:
            Timer()
            last_time = cur_time

        # GUI
        SetGUI(window)

        glfw.swap_buffers(window)

    glfw.destroy_window(window)
    glfw.terminate()

    return 0


if __name__ == "__main__":
    main()
