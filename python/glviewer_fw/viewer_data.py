# -*- coding: utf-8 -*-
"""!
@file viewer_data.py

@brief シミュレーションデータの読み込みと共通処理
       (C++版 glviewer.cpp のデータ読み込み部分を切り出したもの)
       matplotlib版(glviewer_mpl.py)とOpenGL版(glviewer_gl.py)の両方から使う

@author Makoto Fujisawa
@date 2023-05 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import os
import sys
import math
import glob

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from rx_utils import *


# -----------------------------------------------------------------------------
# データ格納用クラス
# -----------------------------------------------------------------------------
class rxData1D:
    """! 1次元スカラーデータ(ある時刻のもの) """
    def __init__(self):
        self.x = []    # 各グリッドの位置
        self.f = []    # 各グリッドでの関数値
        self.n = 0     # グリッド数
        self.t = 0.0   # 時刻
        self.min = [0.0, 0.0, 0.0]
        self.max = [0.0, 0.0, 0.0]

    def SearchRange(self):
        """! データの範囲(最大・最小値)を調べる """
        if not self.x:
            return 1
        self.min[0] = self.max[0] = self.x[0]
        self.min[1] = self.max[1] = self.f[0]
        if len(self.x) != len(self.f):
            print("size of x and f is not matched")
            print("  x : " + str(len(self.x)) + ", f : " + str(len(self.f)))
            return 2
        for i in range(1, len(self.x)):
            if self.x[i] < self.min[0]:
                self.min[0] = self.x[i]
            if self.x[i] > self.max[0]:
                self.max[0] = self.x[i]
            if self.f[i] < self.min[1]:
                self.min[1] = self.f[i]
            if self.f[i] > self.max[1]:
                self.max[1] = self.f[i]
        self.min[2] = self.max[2] = 0
        return 0


class rxData2D:
    """! 2次元スカラーデータ(ある時刻のもの) """
    def __init__(self):
        self.f = []    # 各グリッドでの関数値(1次元配列, k = i+j*(nx+1))
        self.nx = 0
        self.ny = 0
        self.dx = 0.0
        self.dy = 0.0
        self.x0 = 0.0
        self.y0 = 0.0
        self.t = 0.0
        self.min = [0.0, 0.0, 0.0]
        self.max = [0.0, 0.0, 0.0]

    def SearchRange(self):
        """! データの範囲(最大・最小値)を調べる """
        if not self.f:
            return 1
        self.min[2] = self.max[2] = self.f[0]
        for v in self.f:
            if v < self.min[2]:
                self.min[2] = v
            if v > self.max[2]:
                self.max[2] = v
        self.min[0] = self.max[0] = 0
        self.min[1] = self.max[1] = 0
        return 0


class rxTimeData:
    """! 時系列データ(各時刻のデータをまとめたもの) """
    def __init__(self, d):
        self.data = []           # 各時刻のデータ
        self.current_index = 0   # 現在表示中のステップ
        self.min = [0.0, 0.0, 0.0]
        self.max = [1.0, 1.0, 1.0]
        self.fmin = 0.0
        self.fscale = -1.0
        self.state = 0
        self.dim = d

    def SearchRange(self):
        """! データの範囲を調べる """
        for k in range(3):
            self.min[k] = 1e10
            self.max[k] = -1e10
        for d in self.data:
            d.SearchRange()
            for k in range(3):
                if d.min[k] < self.min[k]:
                    self.min[k] = d.min[k]
                if d.max[k] > self.max[k]:
                    self.max[k] = d.max[k]

    def Increment(self, d=1):
        """! 表示ステップを進める """
        self.current_index += d
        if self.current_index >= len(self.data):
            self.current_index = len(self.data)-1
            return False
        return True

    def Decrement(self, d=1):
        """! 表示ステップを戻す """
        self.current_index -= d
        if self.current_index < 0:
            self.current_index = 0
            return False
        return True


# -----------------------------------------------------------------------------
# シミュレーションデータの読み込み
# -----------------------------------------------------------------------------
def Read1ds(filename, sep=","):
    """!
    1次元スカラーデータの読み込み
    @param[in] filename ファイル名
    @param[in] sep 区切り文字
    @return rxData1Dのリスト(読み込めなかったらNone)
    """
    lines = _read_lines_raw(filename, "Read1ds")
    if lines is None:
        return None

    data = []
    for buf in lines:
        d = rxData1D()
        vals = _split_line_fast(buf, sep)
        if not vals:
            break

        # 最初の1つは時刻t
        d.t = vals[0]

        # 位置xとそこでの値fの読み込み(x,fの順に交互に並んでいる)
        for i in range(1, len(vals)-1, 2):
            d.x.append(vals[i])
            d.f.append(vals[i+1])

        if d.x:
            d.n = len(d.x)
            data.append(d)

    if not data:
        return None
    return data


def Read2ds(filename, sep=","):
    """!
    2次元スカラーデータの読み込み
    @param[in] filename ファイル名
    @param[in] sep 区切り文字
    @return rxData2Dのリスト(読み込めなかったらNone)
    """
    lines = _read_lines_raw(filename, "Read2ds")
    if lines is None:
        return None

    data = []
    for buf in lines:
        vals = _split_line_fast(buf, sep)
        if len(vals) < 8:
            break

        d = rxData2D()
        # 最初の1つは時刻t
        d.t = vals[0]
        # 次に続く2つはグリッド分割数nx,ny
        d.nx = int(vals[1])
        d.ny = int(vals[2])
        # 次の2つは左下座標x0,y0
        d.x0 = vals[3]
        d.y0 = vals[4]
        # 次の2つは空間刻み幅dx,dy
        d.dx = vals[5]
        d.dy = vals[6]
        # 関数値fの読み込み
        d.f = vals[7:]

        if d.f:
            data.append(d)

    if not data:
        return None
    return data


def _split_line_fast(buf, sep):
    """!
    1行分の文字列をsepで分割して数値のリストにする
     - rx_utils.py の GetNextString を使った処理と同じことをしているが，
       ビューワでは非常に長い行を大量に読むのでPython標準のsplitを使って高速化している
    @param[in] buf 1行分の文字列
    @param[in] sep 区切り文字
    @return 数値のリスト
    """
    vals = []
    for t in buf.split(sep):
        t = t.strip()
        if t == "":
            continue
        try:
            vals.append(float(t))
        except ValueError:
            pass
    return vals


def _read_lines_raw(filename, funcname):
    """!
    テキストファイルを1行ずつ読み込む(先頭の'#'から始まるヘッダ行は飛ばす)
    @param[in] filename ファイル名
    @param[in] funcname エラー表示用の関数名
    @return 行のリスト(ファイルが開けなかったらNone)
    """
    try:
        f = open(filename, "r", encoding="utf-8", errors="ignore")
    except IOError:
        print(funcname + " : Invalid file specified")
        return None
    lines = []
    for buf in f:
        buf = buf.rstrip("\r\n")
        pos = buf.find('#')
        if pos != -1:
            buf = buf[0:pos]
        if buf.strip() == "":
            continue
        lines.append(buf)
    f.close()
    return lines


def Read(filename, data1d, data2d):
    """!
    データファイルの読み込み
     - 1行目のヘッダ("#1d"/"#2d")でデータの種類を判別する
    @param[in] filename ファイル名
    @param[inout] data1d 1次元データの格納先(rxTimeData)
    @param[inout] data2d 2次元データの格納先(rxTimeData)
    @return (戻り値, データの種類) 戻り値は正常に読み込めたら0
    """
    try:
        f = open(filename, "r", encoding="utf-8", errors="ignore")
    except IOError:
        print("Read : Invalid file specified")
        return 1, 0

    # 1行目だけ読み込んでデータの種類を判別
    buf = f.readline().rstrip("\r\n")
    f.close()

    # ヘッダ行の2つ目以降の値(描画上の最小値とスケール)を取り出す
    hdr = [v.strip() for v in buf.split(",")]

    if buf.find("#1d") != -1:  # 1次元スカラー値データ
        if len(hdr) > 1 and IsNumeric(hdr[1]):
            data1d.fmin = float(hdr[1])    # y軸最小値
        if len(hdr) > 2 and IsNumeric(hdr[2]):
            data1d.fscale = float(hdr[2])  # y軸スケール
        d = Read1ds(filename, ",")
        if d is None:
            return 1, 0
        data1d.data = d
        data1d.state = 1
        data1d.current_index = 0
        data1d.SearchRange()
        return 0, 1
    elif buf.find("#2d") != -1:  # 2次元スカラー値データ
        if len(hdr) > 1 and IsNumeric(hdr[1]):
            data2d.fmin = float(hdr[1])    # 関数値fの描画上の最小値
        if len(hdr) > 2 and IsNumeric(hdr[2]):
            data2d.fscale = float(hdr[2])  # 関数値fの描画上のスケール
        d = Read2ds(filename, ",")
        if d is None:
            return 1, 0
        data2d.data = d
        data2d.state = 1
        data2d.current_index = 0
        data2d.SearchRange()
        return 0, 2

    print("Read : unknown data type (the first line must be #1d or #2d)")
    return 1, 0


def GetFileList(folder, ext="txt"):
    """!
    フォルダ内のファイルリストの取得(C++版 rx_filelist.h の GetFileList に相当)
    @param[in] folder フォルダ名
    @param[in] ext 拡張子
    @return ファイルパスのリスト
    """
    files = glob.glob(os.path.join(folder, "*." + ext))
    files.sort()
    return files


# -----------------------------------------------------------------------------
# 描画用の色
# -----------------------------------------------------------------------------
# 青->水色->緑->黄->赤と変化するサーモグラフ用の基本色
THERMO_COLORS = [(0.0, 0.0, 1.0),
                 (0.0, 1.0, 1.0),
                 (0.0, 1.0, 0.0),
                 (1.0, 1.0, 0.0),
                 (1.0, 0.0, 0.0)]


def CalThermograph(x, xmin=0.0, xmax=1.0):
    """!
    青->緑->赤->白と変化するサーモグラフ用の色生成
    @param[in] x 値
    @param[in] xmin 最小値
    @param[in] xmax 最大値
    @return (r,g,b) 生成された色
    """
    l = xmax-xmin
    if abs(l) < 1e-10:
        return (0.0, 0.0, 0.0)

    ncolors = len(THERMO_COLORS)
    base = THERMO_COLORS
    x = RX_CLAMP((x-xmin)/l, 0.0, 1.0)*(ncolors-1)
    i = int(x)
    if i >= ncolors-1:
        i = ncolors-2
    dx = x-math.floor(x)
    return (RX_LERP(base[i][0], base[i+1][0], dx),
            RX_LERP(base[i][1], base[i+1][1], dx),
            RX_LERP(base[i][2], base[i+1][2], dx))


def default_data_folder():
    """!
    既定のデータフォルダを返す
     - pde フォルダの計算結果(python/pde/data)を既定の読み込み先とする
    @return フォルダのパス
    """
    here = os.path.dirname(os.path.abspath(__file__))
    cand = [os.path.join(here, "data"),
            os.path.join(here, "..", "pde", "data")]
    for c in cand:
        if os.path.isdir(c) and GetFileList(c, "txt"):
            return os.path.normpath(c)
    return os.path.normpath(cand[-1])
