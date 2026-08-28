# -*- coding: utf-8 -*-
"""!
@file rx_utils.py

@brief 数値計算テストの共通モジュール(C++版 rx_utils.h のPython移植)

@author Makoto Fujisawa
@date   2019 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import math
import random
import time


# -----------------------------------------------------------------------------
# 定数
# -----------------------------------------------------------------------------
# 円周率
RX_PI = 3.14159265358979323846

#! 許容誤差
RX_FEQ_EPS = 1.0e-10

#! degree -> radian の変換係数(pi/180.0)
RX_DEGREES_TO_RADIANS = 0.0174532925199432957692369076848

#! radian -> degree の変換係数(180.0/pi)
RX_RADIANS_TO_DEGREES = 57.295779513082320876798154814114


# -----------------------------------------------------------------------------
# 数値の表示
#  - C++のcoutは既定で「有効桁数6桁」で実数を表示する.
#    Pythonのprintは既定でもっと多くの桁を表示してしまうので,
#    C++版と同じ見た目にするために下記のfmt()を通して表示する.
#  - cout.precision(n) に相当するのが set_precision(n).
# -----------------------------------------------------------------------------
_g_precision = 6  # 表示する有効桁数(C++のcoutの既定値と同じ6)


def set_precision(n):
    """!
    表示桁数の設定(C++の cout.precision(n) に相当)
    @param[in] n 有効桁数
    """
    global _g_precision
    _g_precision = n


def get_precision():
    """! 現在の表示桁数を返す """
    return _g_precision


def fmt(x, n=None):
    """!
    C++のcoutと同じ書式で数値を文字列化する
     - 実数は有効桁数n(既定は6)の%g形式, 整数はそのまま
    @param[in] x 数値
    @param[in] n 有効桁数(省略時は set_precision で設定した値)
    @return 文字列
    """
    if n is None:
        n = _g_precision
    if isinstance(x, bool):
        return "1" if x else "0"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, complex):
        # C++の complex は (実部,虚部) の形で表示される
        return "(" + fmt(x.real, n) + "," + fmt(x.imag, n) + ")"
    if isinstance(x, (list, tuple)):
        # C++版では vector に対する << 演算子を定義してカンマ区切りで表示している
        return ", ".join([fmt(v, n) for v in x])
    return "{0:.{1}g}".format(float(x), n)


def TOSTR(x):
    """! 数値を文字列に変換(C++版の RX_TO_STRING/TOSTR に相当) """
    return fmt(x)


RX_TO_STRING = TOSTR


# -----------------------------------------------------------------------------
# マクロ(C++版のテンプレート関数に相当)
# -----------------------------------------------------------------------------
def RX_IS_ZERO(x):
    """! ゼロ判定 """
    return abs(x) < RX_FEQ_EPS


def RX_FEQ(a, b):
    """! 許容誤差を含めた等値判定 """
    return abs(a-b) < RX_FEQ_EPS


def RX_MAX(a, b):
    """! 最大値判定(2値) """
    return a if a > b else b


def RX_MAX3(a, b, c):
    """! 最大値判定(3値) """
    return (a if a > c else c) if a > b else (b if b > c else c)


def RX_MIN(a, b):
    """! 最小値判定(2値) """
    return a if a < b else b


def RX_MIN3(a, b, c):
    """! 最小値判定(3値) """
    return (a if a < c else c) if a < b else (b if b < c else c)


def RX_CLAMP(x, a, b):
    """! 値のクランプ(クランプした値を返す) """
    return a if x < a else (b if x > b else x)


def RX_LERP(a, b, t):
    """! 1次元線型補間 """
    return a + t*(b-a)


def RX_TO_RADIANS(x):
    """! degree -> radian の変換 """
    return x*RX_DEGREES_TO_RADIANS


def RX_TO_DEGREES(x):
    """! radian -> degree の変換 """
    return x*RX_RADIANS_TO_DEGREES


def div(a, b):
    """!
    実数の除算(C++と同じ挙動にするためのもの)
     - Pythonでは 0 で割ると例外(ZeroDivisionError)が発生してプログラムが止まってしまうが，
       C++の実数型(double)では inf や nan になって計算が続く．
       授業中のデモで途中で止まらないよう，C++と同じ挙動にしている．
    @param[in] a,b 被除数,除数
    @return a/b (b=0のときは inf もしくは nan)
    """
    if b == 0.0:
        if a == 0.0:
            return float("nan")
        return float("inf") if a > 0.0 else float("-inf")
    return a/b


def RXGETTIME():
    """! 現在時刻[s]を返す(処理時間計測用, C++版の RXGETTIME に相当) """
    return time.perf_counter()


RXTIME2SEC = 1.0


# -----------------------------------------------------------------------------
# 配列(C++版の vector<double>, vector<vector<double> > に相当)の生成
# -----------------------------------------------------------------------------
def zeros(n, val=0.0):
    """!
    1次元配列の確保 (C++の vector<double> a(n, val) に相当)
    @param[in] n 要素数
    @param[in] val 初期値
    """
    return [val]*n


def zeros2(n, m, val=0.0):
    """!
    2次元配列の確保 (C++の vector<vector<double> > a(n, vector<double>(m, val)) に相当)
     - [[val]*m]*n としてしまうと同じリストがn個並ぶだけになるので注意
    @param[in] n,m 行数,列数
    @param[in] val 初期値
    """
    return [[val]*m for _ in range(n)]


# -----------------------------------------------------------------------------
# テキストファイル処理
# -----------------------------------------------------------------------------
def GetNextString(src, sep, pos):
    """!
    文字列からpos以降で最初の区切り文字までを抽出
     - もし, "(ダブルクオーテーション)で囲まれていたらその範囲を抽出
    @param[in] src 元の文字列
    @param[in] sep 区切り文字
    @param[in] pos 探索開始位置
    @return (抽出文字列, 次の抽出開始位置) 次が無ければ位置として-1を返す
    """
    sub = ""
    extracted = False
    if pos < len(src) and src[pos] == '"':  # ダブルクオーテーションのチェック
        j = src.find('"', pos+1)
        if j != -1:
            sub = src[pos+1:j]
            pos = j+1
            extracted = True

    i = src.find(sep, pos)
    if i == -1:
        if not extracted:
            sub = src[pos:]
        return sub, -1
    else:
        cnt = 1
        while i+cnt < len(src) and src[i+cnt] == ' ':  # sepの後のスペースを消す
            cnt += 1
        if not extracted:
            sub = src[pos:i]
        return sub, (-1 if i+cnt >= len(src) else i+cnt)


def GetFirstString(src, sep):
    """!
    文字列から最初の区切り文字までを抽出
    @param[in] src 元の文字列
    @param[in] sep 区切り文字
    @return (抽出文字列, 次の抽出開始位置)
    """
    return GetNextString(src, sep, 0)


def IsNumeric(s):
    """!
    文字列が実数値を表しているかを調べる
    @param[in] s 文字列
    @return 実数値ならTrue
    """
    if s is None or s.strip() == "":
        return False
    for c in s:
        if c not in "-+0123456789. Ee\t":
            return False
    try:
        float(s)
    except ValueError:
        return False
    return True


def _split_line(buf, sep):
    """!
    1行分の文字列をsepで分割して数値のリストにする(内部用)
    @param[in] buf 1行分の文字列
    @param[in] sep 区切り文字
    @return 数値のリスト
    """
    vals = []
    pos = 0
    while True:
        sub, pos = GetNextString(buf, sep, pos)
        if IsNumeric(sub):
            vals.append(float(sub))
        if pos == -1:
            break
    return vals


def _read_lines(file_name, funcname):
    """!
    テキストファイルを読み込んで, コメント('#'以降)と空行を除いた行のリストを返す(内部用)
    @param[in] file_name ファイル名
    @param[in] funcname エラー表示用の関数名
    @return 行のリスト(ファイルが開けなかったらNone)
    """
    try:
        f = open(file_name, "r", encoding="utf-8", errors="ignore")
    except IOError:
        print(funcname + " : Invalid file specified")
        return None

    lines = []
    for buf in f:
        buf = buf.rstrip("\r\n")
        # '#'以降はコメントとして無視
        comment_start = buf.find('#')
        if comment_start != -1:
            buf = buf[0:comment_start]
        # 空行は無視
        if buf.strip() == "":
            continue
        lines.append(buf)
    f.close()
    return lines


def ReadMatrix(file_name, sep=","):
    """!
    テキストファイルから行列要素を読み込む
    @param[in] file_name ファイル名
    @param[in] sep 区切り文字
    @return 行列要素(2次元リスト). ファイルが開けなかったらNoneを返す
    """
    lines = _read_lines(file_name, "ReadMatrix")
    if lines is None:
        return None

    mat = []
    for buf in lines:
        mat_line = _split_line(buf, sep)
        if mat_line:
            mat.append(mat_line)
    return mat


def ReadAlgebra(file_name, sep=","):
    """!
    テキストファイルから代数方程式の係数を読み込む
    @param[in] file_name ファイル名
    @param[in] sep 区切り文字
    @return 係数列(リスト). ファイルが開けなかったらNoneを返す
    """
    lines = _read_lines(file_name, "ReadAlgebra")
    if lines is None:
        return None

    for buf in lines:
        c = _split_line(buf, sep)
        if c:
            return c  # 最初に見つかった1行だけを使う
    return []


# -----------------------------------------------------------------------------
# 行列処理
# -----------------------------------------------------------------------------
def OutputMatrix(matrix, nx=None, ny=None):
    """!
    行列の画面出力
    @param[in] matrix 行列を格納した2次元リスト
    @param[in] nx,ny  行列の大きさ(省略時はmatrixの大きさをそのまま使う)
    """
    if nx is None:
        nx = len(matrix)
    if ny is None:
        ny = len(matrix[0]) if nx > 0 else 0
    for i in range(nx):
        line = ""
        for j in range(ny):
            line += fmt(matrix[i][j]) + ("" if j == ny-1 else " ")
        print(line)


def MulMatrix(a, b, y, n):
    """!
    2次元リストに格納された正方行列同士の掛け算
    @param[in] a,b nxn行列
    @param[out] y  結果の行列(nxn)
    @param[in] n 行列の大きさ
    @return 1:成功
    """
    for i in range(n):
        for j in range(n):
            y[i][j] = 0.0
            for k in range(n):
                y[i][j] += a[i][k]*b[k][j]
    return 1


def mul_mm(a, b, n):
    """!
    正方行列同士の掛け算(結果を戻り値で返す版)
    @param[in] a,b nxn行列
    @param[in] n 行列の大きさ
    @return 結果の行列(nxn)
    """
    y = zeros2(n, n)
    for i in range(n):
        for j in range(n):
            y[i][j] = 0.0
            for k in range(n):
                y[i][j] += a[i][k]*b[k][j]
    return y


def transpose(a, n):
    """!
    2次元リストに格納された行列の転置
    @param[in] a nxn行列
    @param[in] n 行列の大きさ
    @return 転置行列
    """
    t = zeros2(n, n)
    for i in range(n):
        for j in range(n):
            t[i][j] = a[j][i]
    return t


def MulMatrixVector(a, b, y, n):
    """!
    行列とベクトルの掛け算
    @param[in] a n×n行列
    @param[in] b n次元ベクトル
    @param[out] y 結果のベクトル(n)
    @param[in] n 行列の大きさ
    @return 1:成功
    """
    for i in range(n):
        y[i] = 0.0
        for k in range(n):
            y[i] += a[i][k]*b[k]
    return 1


def mul_mv(a, b, n):
    """!
    行列とベクトルの掛け算(結果を戻り値で返す版)
    @param[in] a n×n行列
    @param[in] b n次元ベクトル
    @param[in] n 行列の大きさ
    @return 結果のベクトル(n)
    """
    y = zeros(n)
    for i in range(n):
        y[i] = 0.0
        for k in range(n):
            y[i] += a[i][k]*b[k]
    return y


def dot(a, b, n):
    """!
    ベクトル同士の内積
    @param[in] a,b n次元ベクトル
    @param[in] n ベクトルの大きさ
    @return 内積
    """
    d = 0.0
    for i in range(n):
        d += a[i]*b[i]
    return d


def mul_sv(a, b, n):
    """!
    ベクトルとスカラー値の掛け算
    @param[in] a スカラー値
    @param[in] b n次元ベクトル
    @param[in] n ベクトルの大きさ
    @return 結果のベクトル(n)
    """
    y = zeros(n)
    for i in range(n):
        y[i] = a*b[i]
    return y


def normalize(a, n):
    """!
    ベクトルの正規化(aは破壊的に書き換えられる)
    @param[inout] a n次元ベクトル
    @param[in] n ベクトルの大きさ
    @return ノルムの値
    """
    l = math.sqrt(dot(a, a, n))  # |a|の計算
    if abs(l) < 1e-10:
        return 0.0
    for i in range(n):
        a[i] /= l
    return l


def unit(a, n):
    """!
    単位行列の設定(aは破壊的に書き換えられる)
    @param[inout] a nxn行列
    @param[in] n 行列の大きさ
    """
    for i in range(n):
        for j in range(n):
            a[i][j] = 1.0 if i == j else 0.0


# -----------------------------------------------------------------------------
# サンプリング点(データ点)
# -----------------------------------------------------------------------------
def MakeSamplingPoints(x0, x1, dx, func):
    """!
    サンプリング点(データ点)の生成(1次元)
    @param[in] x0,x1 サンプリング範囲
    @param[in] dx サンプリング間隔
    @param[in] func 関数値を与える関数
    @return (xi, yi) サンプリングデータ
    """
    xi = []
    yi = []
    x = x0
    while x <= x1:
        xi.append(x)
        yi.append(func(x))
        x += dx
    return xi, yi


def MakeSamplingPointsWithWhiteNoise(x0, x1, dx, func, nwidth):
    """!
    サンプリング点(データ点)の生成(1次元,ホワイトノイズ付き)
    @param[in] x0,x1 サンプリング範囲
    @param[in] dx サンプリング間隔
    @param[in] func 関数値を与える関数
    @param[in] nwidth ノイズのサイズ([-0.5nwidth, 0.5nwidth]のノイズを追加する)
    @return (xi, yi) サンプリングデータ
    """
    random.seed()  # 乱数のシード値を時間によって変える
    xi = []
    yi = []
    x = x0
    while x <= x1:
        noise = (random.random()-0.5)*nwidth
        xi.append(x)
        yi.append(func(x)+noise)
        x += dx
    return xi, yi


def MakeChebyshevNodes(x0, x1, dx, func):
    """!
    サンプリング点(データ点)の生成(1次元)
     - チェビシェフ節点
    @param[in] x0,x1 サンプリング範囲
    @param[in] dx サンプリング間隔
    @param[in] func 関数値を与える関数
    @return (xi, yi) サンプリングデータ
    """
    xi = []
    yi = []
    n = int((x1-x0)/dx)+1
    for i in range(n, 0, -1):
        x = math.cos((2.0*i-1.0)/(2.0*n)*RX_PI)
        x = x0+(x/2.0+0.5)*(x1-x0)
        xi.append(x)
        yi.append(func(x))
    return xi, yi


def OutputSamplingPointsToFile(xi, yi, filename):
    """!
    サンプリング点(データ点)のファイル出力
    @param[in] xi,yi サンプリングデータ
    @param[in] filename 出力ファイル名
    """
    with open(filename, "w") as fo:
        for i in range(len(xi)):
            fo.write(fmt(xi[i]) + ", " + fmt(yi[i]) + "\n")


def OutputSamplingPoints(xi, yi, filename=None):
    """!
    サンプリング点(データ点)の画面出力(filename指定時はファイル出力)
    @param[in] xi,yi サンプリングデータ
    @param[in] filename 出力ファイル名(省略時は画面に出力)
    """
    if filename is not None:
        OutputSamplingPointsToFile(xi, yi, filename)
        return
    s = "sampling points : "
    for i in range(len(xi)):
        s += "(" + fmt(xi[i]) + ", " + fmt(yi[i]) + ")"
        s += "" if i == len(xi)-1 else ",  "
    print(s)


def OutputFunction(x0, x1, dx, func, filename):
    """!
    関数値のファイル出力
    @param[in] x0,x1 サンプリング範囲
    @param[in] dx サンプリング間隔
    @param[in] func 関数値を与える関数, もしくは, 関数値を格納したリスト
    @param[in] filename 出力ファイル名
    @return 生成されたデータ個数
    """
    cnt = 0
    with open(filename, "w") as fo:
        if isinstance(func, (list, tuple)):  # 関数値のリストが渡された場合
            x = x0
            for i in range(len(func)):
                fo.write(fmt(x) + "," + fmt(func[i]) + "\n")
                x += dx
        else:  # 関数が渡された場合
            x = x0
            while x <= x1:
                fo.write(fmt(x) + ", " + fmt(func(x)) + "\n")
                x += dx
                cnt += 1
    return cnt


# -----------------------------------------------------------------------------
# 偏微分方程式のための初期値設定とファイル出力
# -----------------------------------------------------------------------------
def SetValueRectangle(f, n, a, b):
    """!
    関数の初期値として矩形波を設定
    @param[out] f 未知関数fの各グリッドでの値
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @param[in] a,b 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (b-a)/n  # 空間刻み幅
    l = b-a      # 計算範囲全体の長さ
    x = 0.0      # 計算範囲内での相対位置(原点をaとする)
    for i in range(n+1):
        if x >= 0.05*l and x <= 0.25*l:
            f[i] = 1.0
        else:
            f[i] = 0.0
        x += h
    return 0


def SetValueSin(f, n, a, b):
    """!
    関数の初期値として正弦波を設定
    @param[out] f 未知関数fの各グリッドでの値
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @param[in] a,b 計算範囲
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    h = (b-a)/n  # 空間刻み幅
    l = b-a      # 計算範囲全体の長さ
    x = 0.0      # 計算範囲内での相対位置(原点をaとする)
    for i in range(n+1):
        if x >= 0.05*l and x <= 0.45*l:
            f[i] = 0.5*math.sin(RX_PI*(x-0.05*l)/(0.4*l))
        else:
            f[i] = 0.0
        x += h
    return 0


def OutputValueToFile(f, n, a, b, t, fo):
    """!
    関数値をファイル出力(1D)
    @param[in] f 未知関数fの各グリッドでの値
    @param[in] n 計算範囲内での分割数(h=(b-a)/n)
    @param[in] a,b 計算範囲
    @param[in] t 現在の時間(時間軸方向の現在地)
    @param[in] fo 出力先ファイルオブジェクト
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1

    # 結果の出力
    x = a            # 位置
    h = (b-a)/n      # 空間方向の刻み幅
    fo.write(fmt(t) + ",")
    for i in range(n+1):
        fo.write(fmt(x) + "," + fmt(f[i]) + ("" if i == n else ","))
        x += h
    fo.write("\n")
    return 0


def OutputValueToFile2D(f, n, x0, xn, y0, yn, t, fo):
    """!
    関数値をファイル出力(2D)
    @param[in] f 未知関数fの各グリッドでの値(2次元リスト)
    @param[in] n 計算範囲内での分割数
    @param[in] x0,xn,y0,yn 計算範囲
    @param[in] t 現在の時間(時間軸方向の現在地)
    @param[in] fo 出力先ファイルオブジェクト
    @return 問題なければ0を返す
    """
    if n <= 0:
        return 1
    fo.write(fmt(t) + "," + str(n) + "," + str(n) + ",")

    dx = (xn-x0)/n  # 空間方向の刻み幅
    dy = (yn-y0)/n  # 空間方向の刻み幅

    fo.write(fmt(x0) + "," + fmt(y0) + "," + fmt(dx) + "," + fmt(dy) + ",")

    # 結果の出力
    for j in range(n+1):
        for i in range(n+1):
            fo.write(fmt(f[i][j]) + ("" if (j == n and i == n) else ","))
    fo.write("\n")
    return 0


# -----------------------------------------------------------------------------
# 2次元座標データ
# -----------------------------------------------------------------------------
class rxPoint2:
    """! 2次元座標データ """
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y


def Read2d(filename, sep=","):
    """!
    2次元ベクトルデータの読み込み
    @param[in] filename ファイル名
    @param[in] sep 区切り文字
    @return 座標値データ(rxPoint2のリスト). 読み込めなかったらNoneを返す
    """
    lines = _read_lines(filename, "Read2d")
    if lines is None:
        return None

    data = []
    for buf in lines:
        pos = 0
        # 座標値を1つずつ読み込んでいく
        sub, pos = GetNextString(buf, sep, pos)
        if not IsNumeric(sub):
            break
        x = float(sub)
        sub, pos = GetNextString(buf, sep, pos)
        if not IsNumeric(sub):
            break
        y = float(sub)
        data.append(rxPoint2(x, y))

    if not data:
        return None
    return data


def Write2d(filename, data, header="#p2d", sep=","):
    """!
    2次元ベクトルデータの書き出し
    @param[in] filename ファイル名
    @param[in] data 座標値データ(rxPoint2のリスト)
    @param[in] header 1行目に書き出すデータの種類
    @param[in] sep 区切り文字
    @return 正常に書き出せたら0を返す
    """
    if not data:
        return 1

    try:
        f = open(filename, "w")
    except IOError:
        print("Write2d : Invalid file specified")
        return 1

    # 1行目はデータの種類を表す
    f.write(header + "\n")

    # 点データを1座標/行で書き出す
    for p in data:
        f.write(fmt(p.x) + sep + fmt(p.y) + "\n")
    f.close()
    return 0


def SearchRange(data):
    """!
    座標値データの範囲を求める
    @param[in] data 2次元座標値データ(rxPoint2のリスト)
    @return (min, max) 最小,最大座標(それぞれ[x,y]のリスト). データが空ならNone
    """
    if not data:
        return None
    minp = [data[0].x, data[0].y]
    maxp = [data[0].x, data[0].y]
    for p in data:
        if p.x < minp[0]:
            minp[0] = p.x
        if p.x > maxp[0]:
            maxp[0] = p.x
        if p.y < minp[1]:
            minp[1] = p.y
        if p.y > maxp[1]:
            maxp[1] = p.y
    return minp, maxp
