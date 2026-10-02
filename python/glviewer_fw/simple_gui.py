# -*- coding: utf-8 -*-
"""!
@file simple_gui.py

@brief OpenGL(PyOpenGL)+GLFW用の簡易GUI
       C++版で使っていた ImGUI の代わりになる最小限のウィジット群
        - ボタン，チェックボックス，スライダ，テキスト表示
        - ImGUIと同じ「即時モード(immediate mode)GUI」の形で使う:
            gui.begin(window)
            gui.text("hello")
            if gui.button("push me"): ...
            gui.end()
       文字はPillow(PIL)でテクスチャに描いてからOpenGLで貼り付けている．

@author Makoto Fujisawa
@date 2026 (Python版)
"""

# -----------------------------------------------------------------------------
# インポート
# -----------------------------------------------------------------------------
import glfw
from OpenGL.GL import *
from PIL import Image, ImageDraw, ImageFont


class SimpleGUI:
    """! 即時モードGUI(ImGUIの簡易版) """

    def __init__(self, font_size=14):
        self._tex = {}          # 文字列→テクスチャのキャッシュ
        self._font = None
        self._font_size = font_size
        self._prev_mouse = False
        self.want_capture_mouse = False  # GUIパネル上にマウスがあるか

        # 配色
        self.col_panel = (0.15, 0.15, 0.18, 0.85)
        self.col_widget = (0.30, 0.32, 0.40, 1.0)
        self.col_hover = (0.40, 0.45, 0.60, 1.0)
        self.col_active = (0.25, 0.55, 0.85, 1.0)
        self.col_text = (1.0, 1.0, 1.0)

    # -------------------------------------------------------------------------
    # 文字描画用テクスチャ
    # -------------------------------------------------------------------------
    def _get_font(self):
        """! 文字描画用フォントの取得 """
        if self._font is None:
            try:
                self._font = ImageFont.load_default(size=self._font_size)
            except TypeError:
                # 古いPillowではサイズ指定ができない
                self._font = ImageFont.load_default()
        return self._font

    def _get_text_texture(self, s):
        """!
        文字列をテクスチャにして返す(一度作ったものはキャッシュする)
        @param[in] s 文字列
        @return (テクスチャID, 幅, 高さ)
        """
        if s in self._tex:
            return self._tex[s]

        font = self._get_font()
        # 文字列のサイズを調べてから画像を作る
        img = Image.new("RGBA", (1, 1))
        d = ImageDraw.Draw(img)
        box = d.textbbox((0, 0), s, font=font)
        w = max(box[2]-box[0], 1)+2
        h = max(box[3]-box[1], 1)+2
        img = Image.new("RGBA", (w, h), (255, 255, 255, 0))
        d = ImageDraw.Draw(img)
        d.text((1-box[0], 1-box[1]), s, font=font, fill=(255, 255, 255, 255))

        tex = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.tobytes())
        glBindTexture(GL_TEXTURE_2D, 0)

        self._tex[s] = (tex, w, h)
        return self._tex[s]

    # -------------------------------------------------------------------------
    # 基本的な描画
    # -------------------------------------------------------------------------
    def _rect(self, x, y, w, h, col):
        """! 塗りつぶし矩形の描画(ピクセル座標) """
        if len(col) == 3:
            col = (col[0], col[1], col[2], 1.0)
        glColor4f(*col)
        glBegin(GL_QUADS)
        glVertex2f(x, y)
        glVertex2f(x+w, y)
        glVertex2f(x+w, y+h)
        glVertex2f(x, y+h)
        glEnd()

    def _draw_text(self, s, x, y, col=None):
        """! 文字列の描画(ピクセル座標,左上基準) """
        if col is None:
            col = self.col_text
        tex, w, h = self._get_text_texture(s)
        glEnable(GL_TEXTURE_2D)
        glBindTexture(GL_TEXTURE_2D, tex)
        glColor4f(col[0], col[1], col[2], 1.0)
        glBegin(GL_QUADS)
        glTexCoord2f(0, 0); glVertex2f(x, y)
        glTexCoord2f(1, 0); glVertex2f(x+w, y)
        glTexCoord2f(1, 1); glVertex2f(x+w, y+h)
        glTexCoord2f(0, 1); glVertex2f(x, y+h)
        glEnd()
        glBindTexture(GL_TEXTURE_2D, 0)
        glDisable(GL_TEXTURE_2D)
        return w, h

    # -------------------------------------------------------------------------
    # フレームの開始と終了
    # -------------------------------------------------------------------------
    def begin(self, window, x=10, y=10, width=180, title=None):
        """!
        GUIフレームの開始(描画の最後の方で呼ぶ)
        @param[in] window GLFWウィンドウ
        @param[in] x,y パネルの左上位置(ピクセル)
        @param[in] width パネルの幅(ピクセル)
        @param[in] title パネルのタイトル
        """
        self.window = window
        self.winw, self.winh = glfw.get_framebuffer_size(window)
        mx, my = glfw.get_cursor_pos(window)
        # HiDPI環境ではウィンドウ座標とフレームバッファ座標が異なるので補正
        ww, wh = glfw.get_window_size(window)
        self.mx = mx*(self.winw/ww if ww > 0 else 1)
        self.my = my*(self.winh/wh if wh > 0 else 1)
        down = glfw.get_mouse_button(window, glfw.MOUSE_BUTTON_LEFT) == glfw.PRESS
        self.mouse_down = down
        self.mouse_clicked = down and not self._prev_mouse  # 押した瞬間だけTrue
        self._prev_mouse = down

        self.px, self.py = x, y
        self.pw = width
        self.cx, self.cy = x+8, y+8    # 次のウィジットを置く位置
        self._last = None              # 直前に置いたウィジットの位置と大きさ

        self.title = title

        # パネル内にマウスがあるかを前フレームの高さで判定
        h = getattr(self, "_last_panel_h", 0)
        self.want_capture_mouse = (self.px <= self.mx <= self.px+self.pw and
                                   self.py <= self.my <= self.py+h)

        # ピクセル座標系での描画設定
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, self.winw, self.winh, 0, -1, 1)  # 左上原点
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)


    def end(self):
        """! GUIフレームの終了 """
        self._last_panel_h = (self.cy-self.py)+8
        glDisable(GL_BLEND)
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

    def draw_panel_background(self):
        """!
        パネルの背景を描画
         - ウィジットより先に描く必要があるので begin() の直後に呼ぶ
        """
        h = getattr(self, "_last_panel_h", 40)
        self._rect(self.px, self.py, self.pw, h, self.col_panel)
        if self.title:
            self._draw_text(self.title, self.px+8, self.py+6)
            self.cy += 20

    # -------------------------------------------------------------------------
    # ウィジット
    # -------------------------------------------------------------------------
    def _advance(self, w, h):
        """!
        ウィジットを1つ置いた後の描画位置の更新
         - 直前に置いたウィジットの位置と大きさを覚えておき，
           same_line() が呼ばれたときに横に並べられるようにする
        @param[in] w,h 置いたウィジットの大きさ
        """
        self._last = (self.cx, self.cy, w, h)
        self.cy += h+6            # 次の行へ
        self.cx = self.px+8       # 行の先頭へ戻す

    def _begin_widget(self):
        """! これから置くウィジットの左上位置を返す """
        return self.cx, self.cy

    def same_line(self):
        """! 次のウィジットを直前のウィジットと同じ行に置く """
        if self._last is not None:
            lx, ly, lw, lh = self._last
            self.cx = lx+lw+6
            self.cy = ly

    def text(self, s):
        """! テキスト表示 """
        x, y = self._begin_widget()
        w, h = self._draw_text(s, x, y)
        self._advance(w, h)

    def separator(self):
        """! 区切り線 """
        x, y = self._begin_widget()
        self._rect(self.px+6, y+3, self.pw-12, 1, (0.5, 0.5, 0.5, 0.8))
        self._advance(self.pw-12, 4)

    def button(self, label, w=None, h=22):
        """!
        ボタン
        @param[in] label ボタンのラベル
        @param[in] w,h ボタンの大きさ(wを省略するとパネル幅いっぱい)
        @return 押されたらTrue
        """
        if w is None:
            w = self.pw-16
        x, y = self._begin_widget()
        hover = (x <= self.mx <= x+w and y <= self.my <= y+h)
        col = self.col_active if (hover and self.mouse_down) else (self.col_hover if hover else self.col_widget)
        self._rect(x, y, w, h, col)
        tw, th = self._get_text_texture(label)[1:3]
        self._draw_text(label, x+(w-tw)*0.5, y+(h-th)*0.5)
        self._advance(w, h)
        return hover and self.mouse_clicked

    def checkbox(self, label, value, h=18):
        """!
        チェックボックス
        @param[in] label ラベル
        @param[in] value 現在の値
        @return 新しい値
        """
        x, y = self._begin_widget()
        hover = (x <= self.mx <= x+self.pw-16 and y <= self.my <= y+h)
        self._rect(x, y, h, h, self.col_hover if hover else self.col_widget)
        if value:
            self._rect(x+4, y+4, h-8, h-8, self.col_active)
        self._draw_text(label, x+h+6, y+1)
        self._advance(self.pw-16, h)
        if hover and self.mouse_clicked:
            return not value
        return value

    def slider(self, label, value, vmin, vmax, w=None, h=16):
        """!
        スライダ
        @param[in] label ラベル
        @param[in] value 現在の値
        @param[in] vmin,vmax 値の範囲
        @return 新しい値
        """
        if w is None:
            w = self.pw-16
        x, y = self._begin_widget()
        hover = (x <= self.mx <= x+w and y <= self.my <= y+h)
        self._rect(x, y, w, h, self.col_widget)
        t = 0.0 if vmax <= vmin else (value-vmin)/float(vmax-vmin)
        t = min(max(t, 0.0), 1.0)
        self._rect(x, y, w*t, h, self.col_active)
        self._draw_text(label, x+4, y+1)
        self._advance(w, h)
        if hover and self.mouse_down:
            t = (self.mx-x)/float(w)
            t = min(max(t, 0.0), 1.0)
            return vmin+t*(vmax-vmin)
        return value
