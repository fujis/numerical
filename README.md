# numerical
筑波大学情報学群情報メディア創成学類の授業「情報数学C」の
受講生向け説明ページ(HTML)，サンプルコードのためのGitリポジトリです．

This is a git repository for lecture : Mathematics for Informatics C. 
It is including htmls for explanation and sample code written by C++. 

## 各フォルダについて
* docs : 説明ページのHTMLが入っています．Webページとしてみる場合は [https://fujis.github.io/numerical/](https://fujis.github.io/numerical/) にアクセスしてください．
* bin : glviewerの実行用フォルダ
* shared : サンプルプログラムを動かす際の共有ヘッダなどが入っています．VisualStudio以外でコンパイルする際はshared/incをインクルードフォルダとして追加するか，shared/inc内の*.hファイルをcppファイルのあるところにコピーするようにしてください．なお，shared/inc/GL,GLFW,glmとshared/libはglviewerを動かすときだけ必要です．
* src : サンプルソースです．それぞれの授業に対応するフォルダについては以下参照．
  * src/error : 第1回【数値計算の基礎】数の表現,数値誤差,桁落ち
  * src/linearsystem : 第2回【線形連立方程式の直接解法】ガウスの消去法,ピボット選択付きガウス消去法,LU分解,コレスキー分解 ＆ 第3回【線形連立方程式の反復解法】ヤコビ法,ガウス・ザイデル法,SOR法,共役勾配法
  * src/rootfinding : 第4回【非線形方程式の求根問題】 2分法,ニュートン法,初期値と収束性,DKA法
  * src/optimization : 第5回【最適化問題】 黄金分割探索,シンプレックス法,最急降下法,準ニュートン法
  * src/interpolation : 第6回【補間法と回帰分析】 ラグランジュ補間,スプライン補間,最小2乗法
  * src/integration : 第7回【数値積分】 区分求積法,台形公式,シンプソン公式,ガウス型求積公式
  * src/differential : 第8回【常微分方程式の数値解法】 前進/後退オイラー法,ルンゲ・クッタ法,予測子修正子法
  * src/pde : 第9回【偏微分方程式の数値解法】 放物型方程式,楕円型方程式,双曲型方程式
  * src/eigen : 第10回【行列の固有値計算】 べき乗法,ハウスホルダー変換,QR法

  * src/glviewer_fw : 偏微分方程式を数値計算出といた結果を見るためのOpenGL&glfw&imguiで書かれたビューワ．同梱のlib,dllはWindowsのVisual Studio 2022(x64)用なので，これを使いたい場合は自分の環境に合わせたlib,dllをとってくるか，自身でライブラリビルドしてください．結果の可視化のためだけのプログラムなので，Excelなど他のプログラムで結果を可視化するならば特に必要なしです．


## トラブルシューティング
 - Visual Studioで「error MSB8036: Windows SDK バージョン 10.0.xxxxxx.0 が見つかりませんでした。」というエラーが出てビルドできない．  
　⇒ 「プロジェクト」メニューから「プロジェクトの再ターゲット」で「Windows SDK バージョン:」 のところにその環境で対応するバージョンが出るので，問題なければそのまま「OK」をクリック
 - Visual Studioで実行しても特定のプロジェクトしか実行されない場合は，「プロジェクト」メニューから「スタートアップ プロジェクトの構成」で「現在の選択」にチェックしてください．

 - Mac環境下でエラーが出る場合は，以下のコマンドを試してみてください(error.cppや-o errorのところは適宜変えてください)．  
　`g++ error.cpp -O3 -std=c++11 -I../../shared/inc -o error`  
　ヘッダの中で一部C++11の機能を使っています(foreachとか). 上記はインクルードフォルダ(../../shared/inc)も設定しています．
 
 - 自身のPC環境(Windows)で動かしたい人向けにdocs/textフォルダに「Visual Studio Community インストール方法」というPDFを用意してあるので参考にしてください．
