# -*- coding: utf-8 -*-
"""!
@file glpoints_gl.py

@brief 2次元点群データのビューワ(PyOpenGL+GLFW+簡易GUI版)
       C++版の glpoints.cpp (OpenGL+GLFW+ImGUI) をほぼそのまま移植したもの
       主成分分析(PCA)のデモ用．マウス左クリックで点を追加できる．

       使い方: python glpoints_gl.py [データフォルダ]
        - フォルダを省略すると ./points を読みに行く

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
sys.path.append(os.path.join(HERE, "..", "shared"))
# QR法・逆反復法・共分散行列の計算は eigen フォルダのものをそのまま使う
#  (C++版では同じコードが glviewer_fw/qr.cpp にコピーされている)
sys.path.append(os.path.join(HERE, "..", "eigen"))

from rx_utils import *
from viewer_data import GetFileList
from simple_gui import SimpleGUI
from qr import qr, inverse_iteration
from pca import covariance_matrix2


# -----------------------------------------------------------------------------
# 定数・グローバル変数
# -----------------------------------------------------------------------------
g_winw = 800                 # 描画ウィンドウの幅
g_winh = 800                 # 描画ウィンドウの高さ

g_data = []                  # データ格納用(rxPoint2のリスト)
g_min = [0.0, 0.0]
g_max = [1.0, 1.0]

g_dfiles = []                # ファイルリスト
g_dfile_idx = -1

g_lc = 1.0                   # 線の色(1.0で白)
g_margin = 0.05              # 描画範囲の余白

# 固有値・固有ベクトル
g_lambda = []
g_v = []
g_c = rxPoint2(0.0, 0.0)

g_gui = None                 # 簡易GUI


# -----------------------------------------------------------------------------
# データの読み込みと固有値計算
# -----------------------------------------------------------------------------
def ReadPoints(filename):
    """!
    データファイルの読み込み
    @param[in] filename ファイル名
    @return 正常に読み込めたら0を返す
    """
    global g_data, g_min, g_max
    try:
        f = open(filename, "r", encoding="utf-8", errors="ignore")
    except IOError:
        print("Read : Invalid file specified")
        return 1

    # 1行目だけ読み込んでデータの種類を判別
    buf = f.readline().rstrip("\r\n")
    f.close()

    if buf.find("#p2d") == -1:  # 2次元ベクトル値データ以外は読まない
        print("Read : unknown data type (the first line must be #p2d)")
        return 1

    d = Read2d(filename, ",")
    if d is None:
        return 1
    g_data = d

    # ヘッダ行にデータ範囲が書かれていればそれを使う
    hdr = [v.strip() for v in buf.split(",")]
    if len(hdr) > 2 and IsNumeric(hdr[1]) and IsNumeric(hdr[2]):
        g_min[0] = g_min[1] = float(hdr[1])  # データ範囲最小値
        g_max[0] = g_max[1] = float(hdr[2])  # データ範囲最大値
    else:
        g_min, g_max = SearchRange(g_data)

    print("read : " + filename + "  (" + str(len(g_data)) + " points)")
    return 0


def WritePoints(filename):
    """!
    データファイルの書き出し
    @param[in] filename ファイル名
    @return 正常に書き出せたら0を返す
    """
    return Write2d(filename, g_data, "#p2d,0,1", ",")


def CalEigen():
    """!
    点群から固有値・固有ベクトルを計算(主成分分析)
    """
    global g_lambda, g_v, g_c
    if not g_data:
        return

    g_lambda = []
    g_v = []

    # 共分散行列の計算
    n = 2
    A, g_c = covariance_matrix2(g_data)
    A0 = [list(row) for row in A]

    # 固有値
    max_iter = 100
    eps = 1e-6
    A, max_iter, eps = qr(A, n, max_iter, eps, verbose=False)
    g_lambda = [A[i][i] for i in range(n)]

    # 固有値の表示
    s = "e = ("
    for i in range(n):
        s += fmt(g_lambda[i]) + (")" if i == n-1 else ", ")
    print(s)
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))
    print("")

    # 逆反復法で固有ベクトルを計算
    max_iter = 100
    eps = 1e-6
    g_v, max_iter, eps = inverse_iteration(A0, g_lambda, n, max_iter, eps)

    # 固有ベクトルの表示
    for i in range(n):
        s = "v" + str(i+1) + " = ("
        for j in range(n):
            s += fmt(g_v[i][j]) + (")" if j == n-1 else ", ")
        print(s)
    print("iter = " + str(max_iter) + ", eps = " + fmt(eps))


# -----------------------------------------------------------------------------
# 描画関数
# -----------------------------------------------------------------------------
def DrawArrow2D(s, d, scale=0.2):
    """!
    矢印の描画
    @param[in] s 矢印の始点
    @param[in] d 矢印の方向(長さ込み)
    @param[in] scale 傘の部分の大きさ(全体の長さに対する係数)
    """
    origin = [s[0], s[1]]
    dir = [d[0], d[1]]

    length = math.sqrt(dir[0]*dir[0]+dir[1]*dir[1])
    if length > 1e-6:
        dir[0] /= length
        dir[1] /= length

    # ベクトル(1,0)との間の角度
    theta = 180.0/RX_PI*math.acos(RX_CLAMP(dir[0]*(1 if dir[1] > 0 else -1), -1.0, 1.0))

    # 矢印の傘部分の設定
    arrow_x = scale*length         # 軸方向の長さ
    arrow_y = arrow_x*0.363970234  # 軸に垂直な方向の開き量(傘の開き角度40deg)

    glPushMatrix()

    glTranslatef(origin[0], origin[1], 0.0)  # 矢印原点に移動
    glRotatef(theta, 0.0, 0.0, 1.0)          # 矢印方向に回転(z軸中心)

    glBegin(GL_LINES)
    # 軸
    glVertex2d(0.0, 0.0)
    glVertex2d(length, 0.0)

    # 傘
    glVertex2d(length, 0.0)
    glVertex2d(length-arrow_x, arrow_y)
    glVertex2d(length, 0.0)
    glVertex2d(length-arrow_x, -arrow_y)
    glEnd()

    glPopMatrix()


def DrawFrame(xmin, xmax, ymin, ymax):
    """! 2Dグラフの外枠描画 """
    glBegin(GL_LINE_LOOP)
    glVertex2d(xmin, ymin)
    glVertex2d(xmax, ymin)
    glVertex2d(xmax, ymax)
    glVertex2d(xmin, ymax)
    glEnd()


# -----------------------------------------------------------------------------
# アプリケーション制御関数
# -----------------------------------------------------------------------------
def SelectFile(idx):
    """!
    ファイルリストのidx番目のファイルを読み込む
    @param[in] idx ファイル番号
    """
    global g_dfile_idx, g_lambda, g_v
    if not g_dfiles:
        return
    if idx < 0:
        idx = 0
    if idx > len(g_dfiles)-1:
        idx = len(g_dfiles)-1
    g_dfile_idx = idx
    ReadPoints(g_dfiles[g_dfile_idx])
    g_lambda = []
    g_v = []


def reload_file():
    """! データファイル再読み込み """
    global g_lambda, g_v
    if g_dfile_idx >= 0:
        ReadPoints(g_dfiles[g_dfile_idx])
        g_lambda = []
        g_v = []


def clear_points():
    """! データクリア """
    global g_data, g_lambda, g_v
    g_data = []
    g_lambda = []
    g_v = []


def write_points():
    """! データファイル出力 """
    fn = os.path.join(HERE, "points", "points.txt")
    WritePoints(fn)
    print("wrote points to " + fn)


def savedisplay():
    """! 現在の画面描画を画像ファイルとして保存 """
    from PIL import Image
    fn = "data_" + str(g_dfile_idx) + ".png"
    w, h = glfw.get_framebuffer_size(glfw.get_current_context())
    glPixelStorei(GL_PACK_ALIGNMENT, 1)
    buf = glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE)
    img = Image.frombytes("RGB", (w, h), buf).transpose(Image.FLIP_TOP_BOTTOM)
    img.save(fn)
    print("saved the screen image to " + fn)


def Init(folder):
    """!
    初期化関数
    @param[in] folder データフォルダ
    """
    global g_dfiles, g_dfile_idx
    print("OpenGL version: " + glGetString(GL_VERSION).decode())
    print("Renderer: " + glGetString(GL_RENDERER).decode())

    glClearColor(0.0, 0.0, 0.0, 1.0)
    glClearDepth(1.0)
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
        SelectFile(0)


# -----------------------------------------------------------------------------
# OpenGL/GLFWコールバック関数
# -----------------------------------------------------------------------------
def Display():
    """!
    再描画イベントコールバック関数
    """
    xmargin = g_margin*(g_max[0]-g_min[0])  # x方向余白
    ymargin = g_margin*(g_max[1]-g_min[1])  # y方向余白

    # ビューポート,透視変換行列,モデルビュー変換行列の設定
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(g_min[0]-xmargin, g_max[0]+xmargin, g_min[1]-ymargin, g_max[1]+ymargin, -1.0, 1.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    # 描画バッファのクリア
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glDisable(GL_LIGHTING)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

    # 周りのフレーム
    lc = abs(g_lc-0.5)
    glColor3d(lc, lc, lc)
    glLineWidth(2.0)
    DrawFrame(g_min[0], g_max[0], g_min[1], g_max[1])

    # データ描画
    if g_data:
        glColor3d(g_lc, g_lc, g_lc)
        glPointSize(10.0)
        glBegin(GL_POINTS)
        for p in g_data:
            glVertex2d(p.x, p.y)
        glEnd()

    # 固有値/固有ベクトル描画
    if g_lambda:
        cen = [g_c.x, g_c.y]

        # 重心
        glColor3d(1.0, 1.0, 0.0)
        glPointSize(10.0)
        glBegin(GL_POINTS)
        glVertex2d(cen[0], cen[1])
        glEnd()
        glColor3d(g_lc, g_lc, g_lc)
        glPointSize(20.0)
        glBegin(GL_POINTS)
        glVertex2d(cen[0], cen[1])
        glEnd()

        # 固有ベクトル
        #  分散は2乗されているので平方根をとって長さの目安にする
        #  (2倍は見た目を合わせるためのマジックナンバー)
        ld = 2*math.sqrt(abs(g_lambda[0]))
        glColor3d(0.0, 0.0, 1.0)
        glLineWidth(3.0)
        DrawArrow2D(cen, [ld*g_v[0][0], ld*g_v[0][1]], 0.1)
        ld = 2*math.sqrt(abs(g_lambda[1]))
        glColor3d(1.0, 0.0, 0.0)
        glLineWidth(3.0)
        DrawArrow2D(cen, [ld*g_v[1][0], ld*g_v[1][1]], 0.1)


def help_text():
    """! ヘルプテキストを表示 """
    print("<ESC>,'q' key - quit the program\n"
          "left click - add a point\n"
          "'j','k' key - select next/previous file\n"
          "'r' key - read selected file\n"
          "'l' key - show file list\n"
          "'e' key - calculate eigen value/vector\n"
          "'c' key - clear points\n"
          "'w' key - output points as text file\n"
          "'o' key - output as image file\n"
          "'h' key - show this help")


def Keyboard(window, key, scancode, action, mods):
    """!
    キーボードイベント処理関数
    """
    global g_lc
    if action in (glfw.PRESS, glfw.REPEAT):
        if key in (glfw.KEY_ESCAPE, glfw.KEY_Q):
            glfw.set_window_should_close(window, True)
        elif key == glfw.KEY_B:      # 背景色変更
            g_lc = 1-g_lc
            glClearColor(1-g_lc, 1-g_lc, 1-g_lc, 1.0)
        elif key == glfw.KEY_J:      # 次のファイル
            SelectFile(g_dfile_idx+1)
        elif key == glfw.KEY_K:      # 前のファイル
            SelectFile(g_dfile_idx-1)
        elif key == glfw.KEY_R:      # データファイル再読み込み
            reload_file()
        elif key == glfw.KEY_L:      # データファイルリスト表示
            for f in g_dfiles:
                print(f)
        elif key == glfw.KEY_C:      # データクリア
            clear_points()
        elif key == glfw.KEY_W:      # データファイル出力
            write_points()
        elif key == glfw.KEY_E:      # 固有値/固有ベクトル計算
            CalEigen()
        elif key == glfw.KEY_O:      # 画像ファイル出力
            savedisplay()
        elif key == glfw.KEY_H:      # ヘルプ表示
            help_text()


def Mouse(window, button, action, mods):
    """!
    マウスイベント処理関数
     - 左クリックで点を追加
    """
    if g_gui is not None and g_gui.want_capture_mouse:
        return  # GUIパネル上でのマウスイベント時は何もしない
    x, y = glfw.get_cursor_pos(window)
    ww, wh = glfw.get_window_size(window)
    if button == glfw.MOUSE_BUTTON_LEFT and action == glfw.PRESS:
        xmargin = g_margin*(g_max[0]-g_min[0])  # x方向余白
        ymargin = g_margin*(g_max[1]-g_min[1])  # y方向余白
        lx = g_max[0]-g_min[0]+2*xmargin
        ly = g_max[1]-g_min[1]+2*ymargin

        p = rxPoint2()
        p.x = (x/float(ww))*lx+(g_min[0]-xmargin)
        p.y = ((wh-y)/float(wh))*ly+(g_min[1]-ymargin)
        g_data.append(p)


def Resize(window, w, h):
    """!
    リサイズイベント処理関数
    """
    global g_winw, g_winh
    g_winw, g_winh = w, h
    glViewport(0, 0, g_winw, g_winh)


def SetGUI(window):
    """!
    簡易GUIのウィジット配置
     - C++版のImGUIウィンドウに相当
    """
    gui = g_gui
    gui.begin(window, x=10, y=10, width=200, title="Points")
    gui.draw_panel_background()

    gui.text("data:")
    if g_dfile_idx >= 0:
        gui.text(os.path.basename(g_dfiles[g_dfile_idx]))
        if gui.button("next", w=90):
            SelectFile(g_dfile_idx+1)
        gui.same_line()
        if gui.button("prev", w=90):
            SelectFile(g_dfile_idx-1)
    else:
        gui.text("no data")
    if gui.button("reload"):
        reload_file()
    gui.separator()

    gui.text(str(len(g_data)) + " points")
    if gui.button("eigen (PCA)"):
        CalEigen()
    if gui.button("clear points"):
        clear_points()
    if gui.button("write points"):
        write_points()
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

    # データフォルダはコマンドライン引数で指定できる(省略時は ./points )
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "points")
    print("data folder : " + folder)

    if not glfw.init():
        return 1

    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 2)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 1)

    window = glfw.create_window(g_winw, g_winh, "glpoints (PyOpenGL)", None, None)
    if not window:
        glfw.terminate()
        return 1

    glfw.make_context_current(window)
    glfw.swap_interval(1)

    Init(folder)
    g_gui = SimpleGUI()

    help_text()

    # コールバック関数の設定
    glfw.set_key_callback(window, Keyboard)
    glfw.set_mouse_button_callback(window, Mouse)
    glfw.set_framebuffer_size_callback(window, Resize)
    Resize(window, g_winw, g_winh)

    # メインループ
    while not glfw.window_should_close(window):
        glfw.poll_events()
        Display()
        SetGUI(window)
        glfw.swap_buffers(window)

    glfw.destroy_window(window)
    glfw.terminate()

    return 0


if __name__ == "__main__":
    main()
