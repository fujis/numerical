# 数値計算 サンプルコード (Python版)

C++で書かれた数値計算の授業用サンプルコード
([github.com/fujis/numerical](https://github.com/fujis/numerical)) を
Pythonに移植したもので，元のコードの構成とコメントをできるだけそのまま残している
(そのため，pythonとしては非効率なコードになっているところは結構ある)．

- 元のC++コード : `../src/`
- このPython版   : `./`(このフォルダ)

## フォルダ構成

C++版の `src` フォルダと同じ構成になっている．

| フォルダ | 内容 | 主なファイル |
| --- | --- | --- |
| `error/` | 数値の精度と誤差 | `error.py` |
| `linearsystem/` | 連立一次方程式 | `gauss-elimination.py`, `gauss-jordan.py`, `lu.py`, `cholesky.py`, `jacobi.py`, `gauss-seidel.py`, `sor.py`, `conjugate-gradient.py` |
| `rootfinding/` | 求根問題 | `bisection.py`, `newton-raphson.py`, `dka.py` |
| `optimization/` | 最適化 | `goldensection.py`, `gradientdecent.py`, `quasinewton.py`, `simplex.py` |
| `interpolation/` | 補間・近似 | `interpolation.py`, `leastsquare.py`, `spline.py` |
| `integration/` | 数値積分 | `integration.py`, `simpson.py`, `romberg.py`, `gauss.py`, `monte-carlo.py` |
| `differential/` | 常微分方程式 | `eular.py`, `rk.py`, `adams.py`, `predictor-corrector.py`, `odesystem.py` |
| `pde/` | 偏微分方程式 | `poisson.py`, `advection.py`, `diffuse.py`, `wave.py` |
| `eigen/` | 固有値・固有ベクトル | `power.py`, `jacobi.py`, `householder.py`, `qr.py`, `pca.py` |
||||
| `shared/` | 共通モジュール(C++版の`shared/inc`) | `rx_utils.py`, `rx_funcs.py` |
| `glviewer_fw/` | 計算結果の可視化ビューワ | `glviewer_mpl.py`, `glviewer_gl.py`, `glpoints_mpl.py`, `glpoints_gl.py` |

## 実行方法

### Windows (授業デモ用)

`run_menu.bat` をダブルクリックするとメニューが出るので，
番号を選ぶだけで各サンプルを実行できる．

```
run_menu.bat
```

各フォルダの `run.bat` は，そのフォルダのメニューだけを表示する．

```
pde\run.bat
```

可視化ビューワは以下のバッチファイルからも直接起動できる．

```
glviewer_fw\glviewer.bat      … 1D/2Dデータビューワ (matplotlib版)
glviewer_fw\glpoints.bat      … 点群+PCAビューワ    (matplotlib版)
glviewer_fw\glviewer_gl.bat   … 1D/2Dデータビューワ (PyOpenGL+GLFW版)
glviewer_fw\glpoints_gl.bat   … 点群+PCAビューワ    (PyOpenGL+GLFW版)
```

> バッチファイル自体はASCIIのみで書いてあり，日本語のメニュー表示は
> `menu.py` が行っている(cmd.exeがバッチファイル中の日本語を
> 正しく解釈できない場合があるための対策)．

### コマンドラインから直接実行

各スクリプトは，そのスクリプトが置かれているフォルダをカレントディレクトリにして実行．

```
cd linearsystem
python gauss-elimination.py
```

データファイル(`matrix.txt`など)はスクリプトからの相対パスで読むようにしてあるので，
相対パス関係が同じならばデータファイルと共に別のフォルダに移動させても実行できる(はず)．

## 必要なパッケージ

数値計算のコード本体(`error/` ～ `pde/`)は **標準ライブラリだけで動く**．
可視化ビューワだけ追加パッケージが必要．

```
pip install -r requirements.txt
```
もしくは，
```
pip install matplotlib numpy PyOpenGL glfw pillow
```
matplotlib版のみならmatplotlibとnumpyだけで良い(matplotlib版の方が環境依存は少ないかと)．

| パッケージ | 用途 |
| --- | --- |
| `matplotlib`, `numpy` | matplotlib版ビューワ (`*_mpl.py`) |
| `PyOpenGL`, `glfw`, `pillow` | PyOpenGL版ビューワ (`*_gl.py`) |

## 可視化ビューワの使い方

### glviewer (1D/2Dデータのビューワ)

`pde/` のプログラムが `pde/data/` に出力したテキストデータを読んでアニメーション表示．

1. まず計算結果を作る
   ```
   cd pde
   python advection.py
   python diffuse.py
   python wave.py
   ```
2. ビューワを起動する
   ```
   cd ..\glviewer_fw
   python glviewer_mpl.py
   ```
   引数でデータフォルダを指定することもできます．
   ```
   python glviewer_mpl.py ..\pde\data
   ```

キー操作(C++版とほぼ同じ):

| キー | 動作 |
| --- | --- |
| `s` | アニメーションの開始/停止 |
| `Space` | 1ステップだけ進める |
| `j` / `k` | 1ステップ進める/戻る |
| `Shift+j` / `Shift+k` | 10ステップ進める/戻る |
| `n` / `m` | 次/前のデータファイルを開く |
| `r` | 現在のファイルを再読み込み |
| `l` | ファイルリストを表示 |
| `o` | 画面を画像ファイルに保存 |
| `h` | ヘルプ表示 |
| `q`, `ESC` | 終了 |

画面右のボタン(matplotlib版)/左上のパネル(PyOpenGL版)からも同じ操作ができる．

### glpoints (点群 + 主成分分析のビューワ)

`glviewer_fw/points/` の点データを表示し，マウス左クリックで点を追加できる．
`e` キー(もしくは `eigen (PCA)` ボタン)で共分散行列の固有値・固有ベクトルを
QR法+逆反復法で計算し，第1主成分(青)と第2主成分(赤)の矢印を描画．

| キー | 動作 |
| --- | --- |
| 左クリック | 点を追加 |
| `e` | 固有値・固有ベクトルの計算(PCA) |
| `c` | 点をすべて消す |
| `w` | 点を `points/points.txt` に書き出す |
| `j` / `k` | 次/前のデータファイルを開く |
| `r` | 現在のファイルを再読み込み |
| `o` | 画面を画像ファイルに保存 |
| `q`, `ESC` | 終了 |

## C++版との違い(移植上のメモ)

授業でC++版と見比べるときのために，意図的に変えたところをまとめておく．

- **参照渡しの引数** : C++では `int &max_iter` のように参照で結果を返していた部分は，
  Pythonではタプルで返すようにしている(例: `ret, x, max_iter, eps = bisection(...)`)．
  リストや2次元リストで渡しているもの(行列など)は，Pythonでも中身が書き換わるので
  C++版と同じ書き方のまま．
- **表示桁数** : C++の `cout` は既定で有効桁数6桁だが，Pythonの `print` はもっと多く
  表示してしまう．C++版と同じ見た目にするため，`rx_utils.py` の `fmt()` を通して
  表示している．`cout.precision(n)` に相当するのが `set_precision(n)` ．
- **ゼロ除算** : C++の実数演算では `inf`/`nan` になって計算が続くが，Pythonでは
  例外が発生して止まってしまう．特異行列を扱う可能性のある箇所では
  `rx_utils.py` の `div()` を使ってC++と同じ挙動にしている(`eigen/qr.py`など)．
- **単精度(float)** : Pythonの `float` はC++の `double` 相当．単精度の挙動を見る
  `error/error.py` では `struct` で32bitに丸める `to_float32()` を使っている．
- **オーバーロード** : C++で同じ名前だった関数(1次元版/n次元版のニュートン法など)は，
  Pythonでは名前を変えている(`newton` / `newton_nd` など)．
- **`lambda`** : Pythonの予約語なので，`rx_funcs.py` の係数 λ は `lmbda` にしている．
  値を書き換えるときは `import rx_funcs` として `rx_funcs.lmbda = 25` のようにする．
- **2次元の拡散方程式** : クランク・ニコルソン法は係数行列が (n-1)^2 × (n-1)^2 の密行列に
  なるため，C++版と同じ n=50 だとPythonでは非常に時間がかかる．
  そのため，`pde/diffuse.py` の2次元計算はデフォルトを n=20 にしてある．
- **ImGUI** : PyOpenGL版のビューワでは，ImGUIの代わりに `glviewer_fw/simple_gui.py`
  (ボタン・チェックボックス・スライダだけの即時モードGUI)を使っている．
