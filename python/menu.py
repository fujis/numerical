# -*- coding: utf-8 -*-
"""!
@file numerical_menu.py

@brief 授業デモ用のメニュープログラム
       各サンプルコードを番号を選ぶだけで実行できるようにしたもの

       使い方:
         python numerical_menu.py           … 全体のメニューを表示
         python numerical_menu.py pde       … pdeフォルダのメニューだけを表示

       Windowsでは run_numerical.bat (もしくは各フォルダの run.bat) から呼ばれる．

@author Makoto Fujisawa
@date 2026 (Python版)
"""

import os
import sys
import subprocess

ROOT = os.path.dirname(os.path.abspath(__file__))

# フォルダ名 : (見出し, [(スクリプト名, 説明), ...])
MENU = [
    ("error", "数値の精度と誤差", [
        ("error.py", "整数/実数の範囲, 丸め誤差, 打ち切り誤差, 収束判定"),
    ]),
    ("rootfinding", "求根問題", [
        ("bisection.py", "二分法"),
        ("newton-raphson.py", "ニュートン・ラフソン法 (多次元版のデモ)"),
        ("dka.py", "DKA法 (ワイヤストラス法+アバースの初期値)"),
    ]),
    ("linearsystem", "連立一次方程式", [
        ("gauss-elimination.py", "ガウスの消去法"),
        ("gauss-jordan.py", "ガウス・ジョルダン法 (逆行列)"),
        ("lu.py", "LU分解"),
        ("cholesky.py", "コレスキー分解"),
        ("jacobi.py", "ヤコビ反復法"),
        ("gauss-seidel.py", "ガウス・ザイデル反復法"),
        ("sor.py", "SOR法 (加速緩和係数を入力する)"),
        ("conjugate-gradient.py", "共役勾配法 (ICCG法)"),
    ]),
    ("eigen", "固有値・固有ベクトル", [
        ("power.py", "べき乗法"),
        ("jacobi.py", "ヤコビ法"),
        ("householder.py", "ハウスホルダー法"),
        ("qr.py", "QR法 + 逆反復法"),
        ("pca.py", "主成分分析(PCA)"),
    ]),
    ("interpolation", "補間・近似", [
        ("interpolation.py", "線形補間, ラグランジュ補間"),
        ("leastsquare.py", "最小二乗近似"),
        ("spline.py", "3次スプライン補間"),
    ]),
    ("integration", "数値積分", [
        ("integration.py", "区分求積法, 台形公式"),
        ("simpson.py", "シンプソン公式"),
        ("romberg.py", "ロンバーグ法"),
        ("gauss.py", "ガウス型積分公式"),
        ("monte-carlo.py", "モンテカルロ積分"),
    ]),
    ("differential", "常微分方程式", [
        ("eular.py", "オイラー法, ホイン法, 後退オイラー法"),
        ("rk.py", "ルンゲ・クッタ法 (RK3, RK4)"),
        ("adams.py", "アダムス・バッシュホース法"),
        ("predictor-corrector.py", "予測子修正子法"),
        ("odesystem.py", "連立常微分方程式 (単振り子)"),
    ]),
    ("optimization", "最適化", [
        ("goldensection.py", "黄金分割探索法"),
        ("gradientdecent.py", "最急降下法 (alphaを入力する)"),
        ("quasinewton.py", "準ニュートン法 (BFGS)"),
        ("simplex.py", "シンプレックス法"),
    ]),
    ("pde", "偏微分方程式", [
        ("poisson.py", "ポアソン方程式 (中心差分+CG法)"),
        ("advection.py", "移流方程式 (計算結果を data フォルダに出力)"),
        ("diffuse.py", "拡散方程式 (計算結果を data フォルダに出力)"),
        ("wave.py", "波動方程式 (計算結果を data フォルダに出力)"),
    ]),
    ("glviewer_fw", "可視化ビューワ", [
        ("glviewer_mpl.py", "1D/2Dデータのビューワ (matplotlib版)"),
        ("glviewer_gl.py", "1D/2Dデータのビューワ (PyOpenGL+GLFW版)"),
        ("glpoints_mpl.py", "点群+主成分分析のビューワ (matplotlib版)"),
        ("glpoints_gl.py", "点群+主成分分析のビューワ (PyOpenGL+GLFW版)"),
    ]),
]

LINE = "=" * 62


def clear():
    """! 画面のクリア """
    os.system("cls" if os.name == "nt" else "clear")


def ask(prompt):
    """!
    番号の入力
    @param[in] prompt 入力を促す文字列
    @return 入力された文字列(EOFやCtrl+Cの場合は"0")
    """
    try:
        return input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        print("")
        return "0"


def run_script(folder, script):
    """!
    サンプルコードの実行
    @param[in] folder フォルダ名
    @param[in] script スクリプト名
    """
    path = os.path.join(ROOT, folder, script)
    print("")
    print("-" * 62)
    print(" python " + folder + os.sep + script)
    print("-" * 62)
    print("")
    # 各スクリプトは自分のフォルダをカレントディレクトリとして動かす
    subprocess.run([sys.executable, path], cwd=os.path.join(ROOT, folder))
    print("")
    ask("Enterキーで戻る...")


def folder_menu(folder, title, items):
    """!
    1つのフォルダ内のメニュー
    @param[in] folder フォルダ名
    @param[in] title 見出し
    @param[in] items (スクリプト名, 説明)のリスト
    """
    while True:
        clear()
        print(LINE)
        print(" " + title + "  (" + folder + ")")
        print(LINE)
        for i, (script, desc) in enumerate(items):
            print("  {0:2d}) {1:<24} {2}".format(i+1, script, desc))
        print("   0) 戻る")
        print("")
        sel = ask("番号を入力してEnter: ")
        if sel == "0":
            return
        if sel.isdigit() and 1 <= int(sel) <= len(items):
            run_script(folder, items[int(sel)-1][0])


def top_menu():
    """! 全体のメニュー """
    while True:
        clear()
        print(LINE)
        print(" 数値計算 サンプルコード (Python版)")
        print(LINE)
        for i, (folder, title, items) in enumerate(MENU):
            print("  {0:2d}) {1:<16} {2}".format(i+1, folder, title))
        print("   0) 終了")
        print("")
        sel = ask("番号を入力してEnter: ")
        if sel == "0":
            return
        if sel.isdigit() and 1 <= int(sel) <= len(MENU):
            folder, title, items = MENU[int(sel)-1]
            folder_menu(folder, title, items)


def main():
    # 引数でフォルダ名が指定されていたらそのフォルダのメニューだけを表示
    if len(sys.argv) > 1:
        name = sys.argv[1]
        for folder, title, items in MENU:
            if folder == name:
                folder_menu(folder, title, items)
                return 0
        print("unknown folder : " + name)
        return 1

    top_menu()
    return 0


if __name__ == "__main__":
    main()
