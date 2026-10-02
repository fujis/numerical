# -*- coding: utf-8 -*-
"""!
@file glviewer_mpl.py

@brief シミュレーションデータのビューワ(matplotlib版)
       C++版の glviewer.cpp (OpenGL+GLFW+ImGUI) をmatplotlibで書き直したもの

       使い方: python glviewer_mpl.py [データフォルダ]
        - フォルダを省略すると ../pde/data を読みに行く
        - 対応データ: 1行目が #1d もしくは #2d で始まるテキストファイル

@author Makoto Fujisawa
@date 2023-05 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, CheckButtons
from matplotlib.colors import LinearSegmentedColormap

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from viewer_data import *


# -----------------------------------------------------------------------------
# 定数・グローバル変数
# -----------------------------------------------------------------------------
g_winw, g_winh = 12.0, 6.5   # 描画ウィンドウの大きさ(インチ)
g_animation_on = False       # アニメーションON/OFF

g_data1d = rxTimeData(1)     # 1次元データ格納用
g_data2d = rxTimeData(2)     # 2次元データ格納用
g_current_type = 1           # 現在のデータの種類(1:1D, 2:2D)

g_dfiles = []                # データファイルのリスト
g_dfile_idx = -1             # 現在のファイル番号

g_dt = 0.033                 # 画面更新間隔[s]
g_step_inc = 10              # 1回の更新で進めるステップ数
g_fast_step = True           # 高速再生(10ステップずつ進める)

# サーモグラフ用のカラーマップ(C++版の CalThermograph と同じ色変化)
CMAP_THERMO = LinearSegmentedColormap.from_list("thermo", THERMO_COLORS)


# -----------------------------------------------------------------------------
# データの読み込み
# -----------------------------------------------------------------------------
def read_file(idx):
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


def current_time_data():
    """! 現在表示中の種類の時系列データを返す """
    return g_data2d if g_current_type == 2 else g_data1d


# -----------------------------------------------------------------------------
# 描画
# -----------------------------------------------------------------------------
class Viewer:
    """! matplotlibによるビューワ本体 """

    def __init__(self, folder):
        global g_dfiles
        # データファイルリストの取得
        g_dfiles = GetFileList(folder, "txt")
        for f in g_dfiles:
            print(f)
        if not g_dfiles:
            print("there is no data file in " + folder + " !")
        else:
            read_file(0)

        # 図と座標軸の作成
        self.fig = plt.figure(figsize=(g_winw, g_winh))
        self.fig.canvas.manager.set_window_title("glviewer (matplotlib)")
        # 上:グラフ(1D)もしくはサーモ表示(2D)，下:サーモバー(1Dのみ)
        self.ax = self.fig.add_axes([0.06, 0.34, 0.72, 0.60])
        self.axbar = self.fig.add_axes([0.06, 0.20, 0.72, 0.10])
        self.quad = None   # 2D表示用のpcolormesh

        self._make_widgets()

        # イベントの設定
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)

        # アニメーション用タイマー
        self.timer = self.fig.canvas.new_timer(interval=int(g_dt*1000))
        self.timer.add_callback(self.on_timer)
        self.timer.start()

        self.draw()

    # -------------------------------------------------------------------------
    # ウィジット(C++版のImGUIウィンドウに相当)
    # -------------------------------------------------------------------------
    def _make_widgets(self):
        """! ボタン・スライダ等の作成 """
        f = self.fig
        # ファイル選択
        self.b_prev = Button(f.add_axes([0.81, 0.86, 0.08, 0.05]), "prev")
        self.b_next = Button(f.add_axes([0.90, 0.86, 0.08, 0.05]), "next")
        self.b_prev.on_clicked(lambda e: self.select_file(-1))
        self.b_next.on_clicked(lambda e: self.select_file(+1))

        # アニメーション
        self.b_play = Button(f.add_axes([0.81, 0.78, 0.17, 0.05]), "start/stop")
        self.b_step = Button(f.add_axes([0.81, 0.71, 0.17, 0.05]), "run a step")
        self.b_play.on_clicked(lambda e: self.switch_animation())
        self.b_step.on_clicked(lambda e: self.one_step())

        # 高速再生のON/OFF
        self.c_fast = CheckButtons(f.add_axes([0.81, 0.62, 0.17, 0.06]), ["fast"], [g_fast_step])
        self.c_fast.on_clicked(lambda label: self.switch_fast())

        # 再読み込みと画像保存
        self.b_reload = Button(f.add_axes([0.81, 0.54, 0.17, 0.05]), "reload")
        self.b_save = Button(f.add_axes([0.81, 0.47, 0.17, 0.05]), "screenshot")
        self.b_reload.on_clicked(lambda e: self.reload())
        self.b_save.on_clicked(lambda e: self.savedisplay())

        # ステップ選択スライダ
        n = max(len(current_time_data().data)-1, 1)
        self.s_step = Slider(f.add_axes([0.10, 0.08, 0.68, 0.03]), "step", 0, n, valinit=0, valstep=1)
        self.s_step.on_changed(self.on_slider)
        self._slider_updating = False

    def switch_animation(self, on=-1):
        """!
        アニメーションON/OFF
        @param[in] on 1でON, 0でOFF, -1でトグル
        """
        global g_animation_on
        g_animation_on = (not g_animation_on) if on == -1 else bool(on)
        return g_animation_on

    def switch_fast(self):
        """! 高速再生(10ステップずつ進める)のON/OFF """
        global g_fast_step, g_step_inc
        g_fast_step = not g_fast_step
        g_step_inc = 10 if g_fast_step else 1

    def one_step(self):
        """! アニメーションを1ステップだけ進める """
        global g_animation_on
        g_animation_on = True
        self.timer_step()
        g_animation_on = False
        self.draw()

    def select_file(self, inc):
        """!
        リスト内の次/前のファイルを選択して読み込む
        @param[in] inc +1で次,-1で前
        """
        if not g_dfiles:
            return
        read_file(g_dfile_idx+inc)
        self.update_slider_range()
        self.reset_axes()
        self.draw()

    def reload(self):
        """! グラフ用ファイルの再読み込み """
        if g_dfile_idx >= 0:
            read_file(g_dfile_idx)
            self.update_slider_range()
            self.reset_axes()
            self.draw()

    def savedisplay(self):
        """! 現在の画面描画を画像ファイルとして保存 """
        d = current_time_data()
        fn = "data_" + str(g_dfile_idx) + "_step" + str(d.current_index) + ".png"
        self.fig.savefig(fn, dpi=100)
        print("saved the screen image to " + fn)

    def update_slider_range(self):
        """! データ数が変わったときにスライダの範囲を更新 """
        n = max(len(current_time_data().data)-1, 1)
        self.s_step.valmax = n
        self.s_step.ax.set_xlim(0, n)
        self._set_slider(0)

    def _set_slider(self, v):
        """! スライダの値をイベントを発生させずに設定 """
        self._slider_updating = True
        self.s_step.set_val(v)
        self._slider_updating = False

    def reset_axes(self):
        """! データの種類が変わったときに描画をリセット """
        self.ax.cla()
        self.axbar.cla()
        self.quad = None

    # -------------------------------------------------------------------------
    # イベント処理
    # -------------------------------------------------------------------------
    def on_slider(self, val):
        """! スライダ操作時のコールバック """
        if self._slider_updating:
            return
        d = current_time_data()
        if d.data:
            d.current_index = int(RX_CLAMP(int(val), 0, len(d.data)-1))
            self.draw()

    def timer_step(self):
        """!
        タイマーイベント処理関数(ある時間間隔で実行)
         - C++版の Timer() に相当
        """
        global g_animation_on
        if g_animation_on:
            d = current_time_data()
            g_animation_on = d.Increment(g_step_inc)

    def on_timer(self):
        """! アニメーション用タイマーのコールバック """
        if g_animation_on:
            self.timer_step()
            self._set_slider(current_time_data().current_index)
            self.draw()

    def on_key(self, event):
        """!
        キーボードイベント処理関数
         - C++版の Keyboard() に相当
        """
        k = event.key
        inc = g_step_inc
        d = current_time_data()
        if k in ("escape", "q"):
            plt.close(self.fig)
        elif k == "s":                    # SキーでアニメーションON/OFF
            self.switch_animation()
        elif k == " ":                    # スペースキーで1ステップだけ進める
            self.one_step()
        elif k in ("j", "J"):             # データのステップを進める
            d.Increment(10 if k == "J" else 1)
            self._set_slider(d.current_index)
            self.draw()
        elif k in ("k", "K"):             # データのステップを戻す
            d.Decrement(10 if k == "K" else 1)
            self._set_slider(d.current_index)
            self.draw()
        elif k == "n":                    # 次のファイル
            self.select_file(+1)
        elif k == "m":                    # 前のファイル
            self.select_file(-1)
        elif k == "r":                    # 再読み込み
            self.reload()
        elif k == "l":                    # ファイルリスト表示
            for f in g_dfiles:
                print(f)
        elif k == "o":                    # 画像ファイル出力
            self.savedisplay()
        elif k == "h":                    # ヘルプ表示
            help_text()

    # -------------------------------------------------------------------------
    # 描画処理(C++版の Display() に相当)
    # -------------------------------------------------------------------------
    def draw(self):
        """! 現在のステップのデータを描画 """
        if g_current_type == 1 and g_data1d.state and g_data1d.data:
            self.draw1d()
        elif g_current_type == 2 and g_data2d.state and g_data2d.data:
            self.draw2d()
        self.fig.canvas.draw_idle()

    def draw1d(self):
        """! 1Dデータ描画(グラフ+サーモバー) """
        xmin = g_data1d.min[0]
        xmax = g_data1d.max[0]
        if g_data1d.fscale > 0.0:   # データファイルにy軸のスケールが設定されていたらそれを用いる
            ymin = g_data1d.fmin
            ymax = g_data1d.fmin+g_data1d.fscale
        else:                       # 設定されていなかったらデータの最大，最小値を使う
            ymin = g_data1d.min[1]
            ymax = g_data1d.max[1]

        # 関数値データ
        d = g_data1d.data[g_data1d.current_index]

        ymargin = 0.05*(ymax-ymin)

        # グラフ
        self.ax.cla()
        self.ax.plot(d.x, d.f, "-", color="w", lw=2.5)
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin-ymargin, ymax+ymargin)
        self.ax.set_facecolor("k")
        self.ax.grid(True, color="0.3", lw=0.5)
        self.ax.set_ylabel("f(x)")
        self.ax.set_title(_title_text(d.t, g_data1d))

        # サーモバー(関数値を色で表したもの)
        self.axbar.cla()
        self.axbar.set_axis_on()
        v = np.array(d.f)[np.newaxis, :]
        self.axbar.imshow(v, cmap=CMAP_THERMO, vmin=ymin, vmax=ymax, aspect="auto",
                          extent=[d.x[0], d.x[-1], 0.0, 1.0], origin="lower", interpolation="bilinear")
        self.axbar.set_xlim(xmin, xmax)
        self.axbar.set_yticks([])
        self.axbar.set_xlabel("x")

    def draw2d(self):
        """! 2Dデータ描画(サーモ表示) """
        d0 = g_data2d.data[0]
        nx, ny = d0.nx, d0.ny
        dx, dy = d0.dx, d0.dy
        xmin = d0.x0
        xmax = xmin+nx*dx
        ymin = d0.y0
        ymax = ymin+ny*dy
        if g_data2d.fscale > 0.0:   # データファイルにスケールが設定されていたらそれを用いる
            fmin = g_data2d.fmin
            fmax = g_data2d.fmin+g_data2d.fscale
        else:                       # 設定されていなかったらデータの最大，最小値を使う
            fmin = g_data2d.min[2]
            fmax = g_data2d.max[2]

        # 関数値データ(1次元配列に k = i+j*(nx+1) の順で入っている)
        d = g_data2d.data[g_data2d.current_index]
        v = np.array(d.f[0:(nx+1)*(ny+1)]).reshape(ny+1, nx+1)

        self.ax.cla()
        self.ax.imshow(v, cmap=CMAP_THERMO, vmin=fmin, vmax=fmax, aspect="auto",
                       extent=[xmin, xmax, ymin, ymax], origin="lower", interpolation="bilinear")
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.ax.set_aspect("equal")
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_title(_title_text(d.t, g_data2d))

        # 2D表示のときは下のサーモバーは使わない
        self.axbar.cla()
        self.axbar.set_axis_off()


def _title_text(t, td):
    """! グラフ上部に表示する情報の文字列を作る """
    name = os.path.basename(g_dfiles[g_dfile_idx]) if g_dfile_idx >= 0 else "no data"
    return "{0}    step {1}/{2}    t = {3}".format(name, td.current_index, len(td.data)-1, fmt(t))


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


# -----------------------------------------------------------------------------
# メインルーチン
# -----------------------------------------------------------------------------
def main():
    # データフォルダはコマンドライン引数で指定できる(省略時は ../pde/data )
    folder = sys.argv[1] if len(sys.argv) > 1 else default_data_folder()
    print("data folder : " + folder)

    v = Viewer(folder)
    help_text()
    plt.show()

    return 0


if __name__ == "__main__":
    main()
