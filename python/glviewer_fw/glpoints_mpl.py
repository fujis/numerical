# -*- coding: utf-8 -*-
"""!
@file glpoints_mpl.py

@brief 2次元点群データのビューワ(matplotlib版)
       C++版の glpoints.cpp (OpenGL+GLFW+ImGUI) をmatplotlibで書き直したもの
       主成分分析(PCA)のデモ用．マウスクリックで点を追加できる．

       使い方: python glpoints_mpl.py [データフォルダ]
        - フォルダを省略すると ./points を読みに行く

@author Makoto Fujisawa
@date 2023-05 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

import matplotlib.pyplot as plt
from matplotlib.widgets import Button

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)
sys.path.append(os.path.join(HERE, "..", "shared"))
# QR法・逆反復法・共分散行列の計算は eigen フォルダのものをそのまま使う
#  (C++版では同じコードが glviewer_fw/qr.cpp にコピーされている)
sys.path.append(os.path.join(HERE, "..", "eigen"))

from rx_utils import *
from viewer_data import GetFileList
from qr import qr, inverse_iteration
from pca import covariance_matrix2


# -----------------------------------------------------------------------------
# 定数・グローバル変数
# -----------------------------------------------------------------------------
g_winw, g_winh = 8.0, 8.0    # 描画ウィンドウの大きさ(インチ)

g_data = []                  # データ格納用(rxPoint2のリスト)
g_min = [0.0, 0.0]
g_max = [1.0, 1.0]

g_dfiles = []                # データファイルのリスト
g_dfile_idx = -1             # 現在のファイル番号

g_margin = 0.05              # 描画範囲の余白

# 固有値・固有ベクトル
g_lambda = []                # 固有値
g_v = []                     # 固有ベクトル
g_c = rxPoint2(0.0, 0.0)     # データの重心


# -----------------------------------------------------------------------------
# データの読み込みと固有値計算
# -----------------------------------------------------------------------------
def Read(filename):
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


def Write(filename):
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
# ビューワ本体
# -----------------------------------------------------------------------------
class PointsViewer:
    """! matplotlibによる点群ビューワ """

    def __init__(self, folder):
        global g_dfiles
        # データファイルリストの取得
        g_dfiles = GetFileList(folder, "txt")
        for f in g_dfiles:
            print(f)
        if not g_dfiles:
            print("there is no data file in " + folder + " !")
        else:
            select_file(0)

        self.fig = plt.figure(figsize=(g_winw, g_winh))
        self.fig.canvas.manager.set_window_title("glpoints (matplotlib)")
        self.ax = self.fig.add_axes([0.08, 0.08, 0.68, 0.86])

        self._make_widgets()

        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)

        self.draw()

    def _make_widgets(self):
        """! ボタンの作成(C++版のImGUIウィンドウに相当) """
        f = self.fig
        self.b_prev = Button(f.add_axes([0.79, 0.88, 0.09, 0.05]), "prev")
        self.b_next = Button(f.add_axes([0.89, 0.88, 0.09, 0.05]), "next")
        self.b_prev.on_clicked(lambda e: self.select(-1))
        self.b_next.on_clicked(lambda e: self.select(+1))

        self.b_reload = Button(f.add_axes([0.79, 0.80, 0.19, 0.05]), "reload")
        self.b_eigen = Button(f.add_axes([0.79, 0.72, 0.19, 0.05]), "eigen (PCA)")
        self.b_clear = Button(f.add_axes([0.79, 0.64, 0.19, 0.05]), "clear points")
        self.b_write = Button(f.add_axes([0.79, 0.56, 0.19, 0.05]), "write points")
        self.b_save = Button(f.add_axes([0.79, 0.48, 0.19, 0.05]), "screenshot")
        self.b_reload.on_clicked(lambda e: self.reload())
        self.b_eigen.on_clicked(lambda e: self.calc_eigen())
        self.b_clear.on_clicked(lambda e: self.clear())
        self.b_write.on_clicked(lambda e: self.write())
        self.b_save.on_clicked(lambda e: self.savedisplay())

    # -------------------------------------------------------------------------
    # 操作
    # -------------------------------------------------------------------------
    def select(self, inc):
        """! リスト内の次/前のファイルを選択して読み込む """
        global g_lambda, g_v
        if not g_dfiles:
            return
        select_file(g_dfile_idx+inc)
        g_lambda = []
        g_v = []
        self.draw()

    def reload(self):
        """! データファイル再読み込み """
        global g_lambda, g_v
        if g_dfile_idx >= 0:
            Read(g_dfiles[g_dfile_idx])
            g_lambda = []
            g_v = []
            self.draw()

    def calc_eigen(self):
        """! 固有値/固有ベクトル計算 """
        CalEigen()
        self.draw()

    def clear(self):
        """! データクリア """
        global g_data, g_lambda, g_v
        g_data = []
        g_lambda = []
        g_v = []
        self.draw()

    def write(self):
        """! データファイル出力 """
        fn = os.path.join(HERE, "points", "points.txt")
        Write(fn)
        print("wrote points to " + fn)

    def savedisplay(self):
        """! 現在の画面描画を画像ファイルとして保存 """
        fn = "data_" + str(g_dfile_idx) + ".png"
        self.fig.savefig(fn, dpi=100)
        print("saved the screen image to " + fn)

    # -------------------------------------------------------------------------
    # イベント処理
    # -------------------------------------------------------------------------
    def on_click(self, event):
        """!
        マウスイベント処理関数
         - グラフ領域内を左クリックすると点を追加する
        """
        if event.inaxes is not self.ax:
            return
        if event.button != 1:
            return
        g_data.append(rxPoint2(event.xdata, event.ydata))
        self.draw()

    def on_key(self, event):
        """!
        キーボードイベント処理関数
         - C++版の Keyboard() に相当
        """
        k = event.key
        if k in ("escape", "q"):
            plt.close(self.fig)
        elif k == "j":     # 次のファイル
            self.select(+1)
        elif k == "k":     # 前のファイル
            self.select(-1)
        elif k == "r":     # データファイル再読み込み
            self.reload()
        elif k == "l":     # データファイルリスト表示
            for f in g_dfiles:
                print(f)
        elif k == "c":     # データクリア
            self.clear()
        elif k == "w":     # データファイル出力
            self.write()
        elif k == "e":     # 固有値/固有ベクトル計算
            self.calc_eigen()
        elif k == "o":     # 画像ファイル出力
            self.savedisplay()
        elif k == "h":     # ヘルプ表示
            help_text()

    # -------------------------------------------------------------------------
    # 描画処理(C++版の Display() に相当)
    # -------------------------------------------------------------------------
    def draw(self):
        """! 点群と固有ベクトルの描画 """
        xmargin = g_margin*(g_max[0]-g_min[0])  # x方向余白
        ymargin = g_margin*(g_max[1]-g_min[1])  # y方向余白

        self.ax.cla()
        self.ax.set_facecolor("k")
        self.ax.set_xlim(g_min[0]-xmargin, g_max[0]+xmargin)
        self.ax.set_ylim(g_min[1]-ymargin, g_max[1]+ymargin)
        self.ax.set_aspect("equal")

        # 周りのフレーム
        self.ax.plot([g_min[0], g_max[0], g_max[0], g_min[0], g_min[0]],
                     [g_min[1], g_min[1], g_max[1], g_max[1], g_min[1]],
                     "-", color="0.5", lw=2)

        # データ描画
        if g_data:
            self.ax.plot([p.x for p in g_data], [p.y for p in g_data],
                         "o", color="w", ms=6)

        # 固有値/固有ベクトル描画
        if g_lambda:
            cen = [g_c.x, g_c.y]

            # 重心
            self.ax.plot([cen[0]], [cen[1]], "o", color="w", ms=16)
            self.ax.plot([cen[0]], [cen[1]], "o", color="y", ms=8)

            # 固有ベクトル
            #  分散は2乗されているので平方根をとって長さの目安にする
            #  (2倍は見た目を合わせるためのマジックナンバー)
            colors = ["b", "r"]
            w = max(g_max[0]-g_min[0], g_max[1]-g_min[1])  # 矢印の傘のサイズの基準
            for i in range(2):
                ld = 2*math.sqrt(abs(g_lambda[i]))
                self.ax.arrow(cen[0], cen[1], ld*g_v[i][0], ld*g_v[i][1],
                              color=colors[i], lw=2.0,
                              head_width=0.03*w, head_length=0.06*w,
                              length_includes_head=True)

        name = os.path.basename(g_dfiles[g_dfile_idx]) if g_dfile_idx >= 0 else "no data"
        self.ax.set_title(name + "    " + str(len(g_data)) + " points")
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.fig.canvas.draw_idle()


def select_file(idx):
    """!
    ファイルリストのidx番目のファイルを読み込む
    @param[in] idx ファイル番号
    """
    global g_dfile_idx
    if not g_dfiles:
        return
    if idx < 0:
        idx = 0
    if idx > len(g_dfiles)-1:
        idx = len(g_dfiles)-1
    g_dfile_idx = idx
    Read(g_dfiles[g_dfile_idx])


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


# -----------------------------------------------------------------------------
# メインルーチン
# -----------------------------------------------------------------------------
def main():
    # データフォルダはコマンドライン引数で指定できる(省略時は ./points )
    folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "points")
    print("data folder : " + folder)

    v = PointsViewer(folder)
    help_text()
    plt.show()

    return 0


if __name__ == "__main__":
    main()
