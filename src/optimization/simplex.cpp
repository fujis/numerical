/*!
@file simplex.cpp

@brief シンプレックス法

@author Makoto Fujisawa
@date 2019-07
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"
#include <iomanip>


//-----------------------------------------------------------------------------
// シンプレックス法(線形計画問題の解法)
//-----------------------------------------------------------------------------
/*!
 * シンプレックス法(simplex method)
 	- 線形計画問題の解法
 	- 二段階シンプレックスではないので原点が解に含まれない場合は適応できない
 * @param[in] a 係数リスト(条件式の係数を格納した2次元配列,最後の1行は最適化式)
 * @param[in] eqn 格式の等号不等号(1:<=, -1:>=)
 * @param[in] n_cond 条件式の数
 * @param[in] m 変数の数
 * @param[out] ans 解
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @return 
 */
int simplex(const vector< vector<double> > &a, const vector<int> &eqn, int n_cond, int m, vector<double> &ans, int &max_iter)
{
	int n = n_cond+1; // 式の数(条件式の数+最適化式の数)
	int m_slack = n_cond; // スラック変数の数
	int m_all = m+m_slack; // スラック変数を含む全変数の数

	// 単体表の作成
	vector< vector<double> > s; // 単体表
	s.resize(n); // 単体表の行数は条件式の数+1=n，最後の行が最適化式
	vector<int> xi(n);	// それぞれの行の基底変数(x1=0,x2=1,x3=2,...とインデックス値を格納)
	for(int i = 0; i < n; ++i){
		s[i].resize(m_all+1); // 単体表の列数はスラック変数を含む全変数の数+1
		if(i != n-1){ // 条件式の行(0～n-2行)
			xi[i] = i+m; // 初期基底変数はスラック変数
			int sgn = (a[i][m] < 0 ? -1 : 1); // 条件式の右辺項の符号
			for(int j = 0; j < m; ++j) s[i][j] = sgn*a[i][j];						// 変数の係数を入れていく(右辺項の符号が正になるようにしている)
			for(int j = m; j < m_all; ++j) s[i][j] = (xi[i] == j ? sgn*eqn[i] : 0);	// スラック変数部分の初期係数は単位行列のような形になる(ただし条件式の符号が>=の場合は-1をセット)
			s[i][m_all] = sgn*a[i][m];												// 右辺項を初期基底可能解としてセット
		}
		else{	// 最終行は最適化式(ここだけ別で設定)
			xi[i] = -1; // 最適化式の変数インデックスには-1を格納しておく
			for(int j = 0; j < m; ++j) s[i][j] = -a[i][j];	// 最終行の最適化式は係数の符号を反転
			for(int j = m; j < m_all+1; ++j) s[i][j] = 0;	// 最後の行の最適化式では初期基底可能解に0をセット
		}
	}

	int k;
	for(k = 0; k < max_iter; ++k){
		// 非基底変数(初期状態ではスラック変数以外)から負で絶対値最大のものを選択
		int b = -1;
		double xmax = 0.0;
		for(int j = 0; j < m_all; ++j){
			if(s[n-1][j] < 0 && -s[n-1][j] > xmax){
				b = j; xmax = -s[n-1][j];
			}
		}
		if(b == -1) break; // 負の変数が見つからなければ収束したとしてループロ抜ける

		// ピボット要素の探索
		int p = 0;
		double ymin = s[0][m_all]/s[0][b];
		for(int i = 1; i < n-1; ++i){
			// 基底可能解をb列の値で割った値が最小のものを選択
			double y = s[i][m_all]/s[i][b];
			if(y < ymin){
				p = i; ymin = y;
			}
		}

		// ピボット要素行を選択された非基底変数(b)で置き換えて，ピボット要素でその行を割る
		xi[p] = b;
		double xp = s[p][b];	// ループ内で値が変わるので一時変数に値を確保しておく
		for(int j = 0; j < m_all+1; ++j){
			s[p][j] /= xp;
		}

		// ピボット行(p)以外の行を選択された非基底変数(b)の値が0になるようにピボット行の係数をx倍して引く
		//  -> s[p][b]は1になっているのでs[p][j]にs[i][b]を掛けて引けば良い
		for(int i = 0; i < n; ++i){
			if(i == p) continue; // ピボット行は飛ばす
			xp = s[i][b];	// ループ内で値が変わるので一時変数に値を確保しておく
			for(int j = 0; j < m_all+1; ++j){
				s[i][j] -= s[p][j]*xp;
			}
		}
	}

	// 解の格納
	for(int i = 0; i < n-1; ++i){
		if(xi[i] < m){
			ans[xi[i]] = s[i][m_all];
		}
	}
	ans[m] = s[n-1][m_all];

	// チェック用
	for(int i = 0; i < n; ++i){
		if(i == n-1) cout << " z: ";
		else cout << "x" << xi[i]+1 << ": ";
		for(int j = 0; j < m_all+1; ++j){
			cout << setw(8) << right << s[i][j] << " ";
		}
		cout << endl;
	}

	max_iter = k;

	return 0;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	//double a0[4][3] = { { 3, 1, 9 },{ 2.5, 2, 12.5 },{ 1, 2, 8 },{ 3, 2, 0 } };
	//double a0[4][3] = { { 1, 2, 800 },{ 3, 4, 1800 },{ 3, 1, 1500 },{ 20, 30, 0 } };
	double a0[3][3] = { { 2, 1, 8 },{ 1, 3, 9 },{ 1, 1, 0 } };
	//int eqn0[4] = { 1, 1, 1, 0 };
	int eqn0[3] = { 1, 1, 0 };

	int n_cond = 2; // 条件式の数
	int m = 2; // 変数の数
	vector<double> x(m+1, 0.0);
	vector< vector<double> > a(n_cond+1);
	vector<int> eqn(n_cond+1, 0);

	for(int i = 0; i < n_cond+1; ++i){
		eqn[i] = eqn0[i];
		a[i].resize(m+1);
		for(int j = 0; j < m+1; ++j){
			a[i][j] = a0[i][j];
		}
	}

	// シンプレックス法で線形計画問題を解く
	int max_iter = 100;
	simplex(a, eqn, n_cond, m, x, max_iter);

	// 結果の画面表示
	cout << "f(";
	for(int i = 0; i < m; ++i) cout << x[i] << (i == m-1 ? ") = " : ",");
	cout << x[m] << ", ";
	cout << "iter = " << max_iter << endl;
	cout << endl;

	return 0;
}


