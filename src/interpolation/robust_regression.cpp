/*!
  @file robust_regression.cpp

  @brief ロバスト推定(ノイズや外れ値を含むデータへの多項式フィッティング)
         - 最小二乗法(比較用)
         - 最小メディアン法(LMedS : Least Median of Squares)
         - RANSAC(RANdom SAmple Consensus)
         - M推定(Huber, Tukey biweight．反復重み付き最小二乗法(IRLS)で解く)

  @author Makoto Fujisawa
  @date 2026-11
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"
#include "rx_funcs.h"

#include <random>


//-----------------------------------------------------------------------------
// LU分解(最小自乗法での係数計算に用いる．leastsquare.cppと同じ)
//-----------------------------------------------------------------------------
/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
 *  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
 *  - LとUを一つの行列にまとめた形で結果を返す
 * @param[inout] A n×nの係数行列．LU分解した結果を格納する．
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int LUDecomp(vector< vector<double> > &A, int n)
{
	if(n <= 0) return 1;

	for(int i = 0; i < n; ++i){
		// l_ijの計算(i >= j)
		for(int j = 0; j <= i; ++j){
			double lu = A[i][j];
			for(int k = 0; k < j; ++k){
				lu -= A[i][k]*A[k][j];	// l_ik * u_kj
			}
			A[i][j] = lu;
		}

		// u_ijの計算(i < j)
		for(int j = i+1; j < n; ++j){
			double lu = A[i][j];
			for(int k = 0; k < i; ++k){
				lu -= A[i][k]*A[k][j];	// l_ik * u_kj
			}
			A[i][j] = lu/A[i][i];
		}
	}

	return 0;
}

/*!
 * LU分解した行列a(n×n)から前進代入・後退代入によりA・x=bを解く
 * @param[in] A LU分解された行列
 * @param[in] b 右辺ベクトル
 * @param[out] x 結果ベクトル
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int LUSolver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n)
{
	if(n <= 0) return 1;

	// 前進代入(forward substitution)
	//  LY=bからYを計算
	for(int i = 0; i < n; ++i){
		double bly = b[i];
		for(int j = 0; j < i; ++j){
			bly -= A[i][j]*x[j];
		}
		x[i] = bly/A[i][i];
	}

	// 後退代入(back substitution)
	//  UX=YからXを計算
	for(int i = n-1; i >= 0; --i){
		double yux = x[i];
		for(int j = i+1; j < n; ++j){
			yux -= A[i][j]*x[j];
		}
		x[i] = yux;
	}

	return 0;
}


//-----------------------------------------------------------------------------
// 共通の関数
//-----------------------------------------------------------------------------
/*!
 * 1次元n次多項式を計算
 * @param[in] x 値が必要な位置
 * @param[in] c 多項式の係数
 * @param[in] n 多項式の次数
 * @return 多項式を計算した結果
 */
double polynominal1d(double x, const vector<double> &c, int n)
{
	double fx = 0;
	double xb = 1;
	for(int i = 0; i <= n; ++i){
		fx += c[i]*xb; xb *= x;
	}
	return fx;
}

/*!
 * 重み付き最小二乗法による多項式フィッティング
 *  - 1次元データ(xi,fi)に対して，Σ w_k (f(x_k) - fi_k)^2 を最小にする1次元n次多項式の係数を求める
 *  - w_k = 1 なら通常の最小二乗法(leastsquare.cpp の leastsquare_interpolation と同じ)
 *  - w_k = 0 のデータは使われない(RANSAC，最小メディアン法で選択したものだけを使うときなど)
 * @param[in] fi 関数値を格納した配列
 * @param[in] xi 関数値に対応する位置を格納した配列
 * @param[in] w 各データの重み(空の場合は全て1)
 * @param[in] m データ数
 * @param[in] n フィッティングする多項式の次数
 * @param[out] c 係数ベクトル
 * @return 0:成功,1:失敗(重みが正のデータ数が係数の数より少ない)
 */
int weighted_leastsquare(const vector<double> &fi, const vector<double> &xi, const vector<double> &w, int m, int n, vector<double> &c)
{
	int dim_b = n+1; // 基底ベクトルの次元数(1次元の場合)

	// 重みが正のデータ数を確認
	int cnt = 0;
	for(int k = 0; k < m; ++k) if(w.empty() || w[k] > 0) cnt++;
	if(cnt < dim_b) return 1;

	// Ac=y, bは基底ベクトル
	c.assign(dim_b, 0.0);	// 結果の係数ベクトル
	vector<double> b(dim_b, 0.0), y(dim_b, 0.0); // 基底ベクトルと右辺項
	vector< vector<double> > A(dim_b, b);

	// 係数行列Aと右辺項yの計算
	for(int k = 0; k < m; ++k){ // データ分だけ反復
		double wk = (w.empty() ? 1.0 : w[k]);
		if(wk <= 0) continue;

		// 多項式の基底ベクトルの計算
		double xb = 1;
		for(int i = 0; i < dim_b; ++i){
			b[i] = xb; xb *= xi[k];
		}

		// 係数行列Aと右辺項bの計算(重みを掛ける)
		for(int i = 0; i < dim_b; ++i){
			for(int j = 0; j < dim_b; ++j){
				A[i][j] += wk*b[i]*b[j];
			}
			y[i] += wk*b[i]*fi[k];
		}
	}

	// Ac=yをLU分解で解く
	LUDecomp(A, dim_b);
	LUSolver(A, y, c, dim_b);

	return 0;
}

/*!
 * 各データの残差 r_k = fi_k - f(x_k) を計算
 * @param[in] fi,xi データ
 * @param[in] m データ数
 * @param[in] c 多項式の係数
 * @param[in] n 多項式の次数
 * @param[out] r 残差
 */
void residuals(const vector<double> &fi, const vector<double> &xi, int m, const vector<double> &c, int n, vector<double> &r)
{
	r.resize(m);
	for(int k = 0; k < m; ++k){
		r[k] = fi[k]-polynominal1d(xi[k], c, n);
	}
}

/*!
 * 中央値(メディアン)を計算
 * @param[in] v 値の配列(コピーして並べ替えるので元の配列は変わらない)
 * @return 中央値
 */
double median(vector<double> v)
{
	int m = (int)v.size();
	sort(v.begin(), v.end());
	return (m%2 == 1 ? v[m/2] : 0.5*(v[m/2-1]+v[m/2]));
}

/*!
 * m個のデータから重複なしでs個をランダムに選ぶ
 * @param[in] m データ数
 * @param[in] s 選ぶ個数
 * @param[in] rng 乱数生成器
 * @param[out] idx 選んだデータのインデックス
 */
void random_sample(int m, int s, std::mt19937 &rng, vector<int> &idx)
{
	idx.clear();
	std::uniform_int_distribution<int> dist(0, m-1);
	while((int)idx.size() < s){
		int k = dist(rng);
		if(find(idx.begin(), idx.end(), k) == idx.end()) idx.push_back(k);
	}
}


//-----------------------------------------------------------------------------
// ロバスト推定法
//-----------------------------------------------------------------------------
/*!
 * 最小メディアン法(LMedS : Least Median of Squares)
 *  - 残差の2乗の「和」ではなく「中央値」を最小にする係数を求める
 *  - 中央値はデータの半分までが外れ値でも影響を受けにくい(崩壊点50%)
 *  - 中央値の最小化は微分で解けないので，ランダムに選んだn+1点を通る多項式を候補として，中央値が最小になるものを探す
 *  - 最後に，求めた中央値から残差のばらつきσを推定して，|r| <= 2.5σ のインライアデータだけで最小二乗法をやり直す
 * @param[in] fi,xi データ
 * @param[in] m データ数
 * @param[in] n フィッティングする多項式の次数
 * @param[in] max_trial 候補を作る回数
 * @param[out] c 係数ベクトル
 * @param[out] inlier 各データがインライアなら1，外れ値なら0
 * @param[in] seed 乱数のシード
 * @return 最小となった残差2乗の中央値
 */
double lmeds(const vector<double> &fi, const vector<double> &xi, int m, int n, int max_trial, vector<double> &c, vector<int> &inlier, unsigned int seed = 1)
{
	std::mt19937 rng(seed);
	int s = n+1; // 多項式を決めるのに必要なデータ数
	double best_med = 1e30;
	vector<double> cs, r, r2(m), w(m);
	vector<int> idx;
	for(int t = 0; t < max_trial; ++t){
		// ランダムに選んだs点だけを使って多項式を決める(s点を通る多項式)
		random_sample(m, s, rng, idx);
		w.assign(m, 0.0);
		for(int k : idx) w[k] = 1.0;
		if(weighted_leastsquare(fi, xi, w, m, n, cs)) continue;

		// 全データの残差の2乗の中央値
		residuals(fi, xi, m, cs, n, r);
		for(int k = 0; k < m; ++k) r2[k] = r[k]*r[k];
		double med = median(r2);
		if(med < best_med){
			best_med = med; c = cs;
		}
	}

	// ロバストな標準偏差の推定(Rousseeuw & Leroy の式)と，それを使ったインライアの判定(インライア:inlier ⇔ アウトライア:outlier(外れ値))
	double sigma = 1.4826*(1.0+5.0/(m-s))*sqrt(best_med);
	residuals(fi, xi, m, c, n, r);
	inlier.assign(m, 0);
	for(int k = 0; k < m; ++k){
		w[k] = (fabs(r[k]) <= 2.5*sigma ? 1.0 : 0.0);
		inlier[k] = (int)w[k];
	}

	// インライアだけで最小二乗法
	weighted_leastsquare(fi, xi, w, m, n, c);

	return best_med;
}

/*!
 * RANSAC(RANdom SAmple Consensus)
 *  - ランダムに選んだn+1点を通る多項式を候補として，残差が閾値以下のデータ(インライア)の数が最大の候補を選ぶ
 *  - 最後にインライアだけで最小二乗法をやり直す
 *  - 閾値(どのくらいのずれまでを正しいデータとみなすか)はあらかじめ与える必要がある
 * @param[in] fi,xi データ
 * @param[in] m データ数
 * @param[in] n フィッティングする多項式の次数
 * @param[in] max_trial 候補を作る回数
 * @param[in] threshold インライアとみなす残差の大きさ
 * @param[out] c 係数ベクトル
 * @param[out] inlier 各データがインライアなら1，外れ値なら0
 * @param[in] seed 乱数のシード
 * @return インライアの数
 */
int ransac(const vector<double> &fi, const vector<double> &xi, int m, int n, int max_trial, double threshold, vector<double> &c, vector<int> &inlier, unsigned int seed = 1)
{
	std::mt19937 rng(seed);
	int s = n+1; // 多項式を決めるのに必要なデータ数
	int best_cnt = -1;
	double best_err = 1e30;
	vector<double> cs, r, w(m);
	vector<int> idx;
	for(int t = 0; t < max_trial; ++t){
		// ランダムに選んだs点だけを使って多項式を決める
		random_sample(m, s, rng, idx);
		w.assign(m, 0.0);
		for(int k : idx) w[k] = 1.0;
		if(weighted_leastsquare(fi, xi, w, m, n, cs)) continue;

		// インライアの数を数える(同じ数なら残差の2乗和が小さい方)
		residuals(fi, xi, m, cs, n, r);
		int cnt = 0;
		double err = 0.0;
		for(int k = 0; k < m; ++k){
			if(fabs(r[k]) <= threshold){ cnt++; err += r[k]*r[k]; }
		}
		if(cnt > best_cnt || (cnt == best_cnt && err < best_err)){
			best_cnt = cnt; best_err = err; c = cs;
		}
	}

	// インライアだけで最小二乗法
	residuals(fi, xi, m, c, n, r);
	inlier.assign(m, 0);
	for(int k = 0; k < m; ++k){
		w[k] = (fabs(r[k]) <= threshold ? 1.0 : 0.0);
		inlier[k] = (int)w[k];
	}
	weighted_leastsquare(fi, xi, w, m, n, c);

	return best_cnt;
}

/*!
 * M推定で使う重み関数(Huber)
 *  - |u| <= k では最小二乗法と同じ(重み1)，それより大きいと重みが 1/|u| に比例して小さくなる
 * @param[in] u 残差をスケールσで割った値
 * @param[in] k 調整定数(正規分布のとき効率95%になる 1.345 がよく使われる)
 */
double weight_huber(double u, double k = 1.345)
{
	return (fabs(u) <= k ? 1.0 : k/fabs(u));
}

/*!
 * M推定で使う重み関数(Tukey biweight)
 *  - |u| <= k では (1-(u/k)^2)^2，それより大きいと重み0(完全に無視する)
 * @param[in] u 残差をスケールσで割った値
 * @param[in] k 調整定数(正規分布のとき効率95%になる 4.685 がよく使われる)
 */
double weight_tukey(double u, double k = 4.685)
{
	if(fabs(u) > k) return 0.0;
	double a = 1.0-(u/k)*(u/k);
	return a*a;
}

/*!
 * M推定(反復重み付き最小二乗法 IRLS : Iteratively Reweighted Least Squares)
 *  - Σ ρ(r_k/σ) を最小にする係数を求める(ρ(u) = u^2/2 なら最小二乗法)
 *  - 残差が大きいデータほど重みを小さくした重み付き最小二乗法を，重みを更新しながら繰り返す
 *  - スケールσは残差の中央絶対偏差(MAD)から推定する : σ = median(|r - median(r)|)/0.6745
 *  - 初期値は通常の最小二乗法(外れ値が多いと初期値が悪くなり，うまくいかないことがある)
 * @param[in] fi,xi データ
 * @param[in] m データ数
 * @param[in] n フィッティングする多項式の次数
 * @param[in] wfunc 重み関数(weight_huber, weight_tukey)
 * @param[in] k 重み関数の調整定数(Huber : 1.345，Tukey : 4.685)
 * @param[out] c 係数ベクトル
 * @param[out] w 最終的な各データの重み
 * @param[inout] max_iter 最大反復数(反復回数を返す)
 * @param[inout] eps 許容誤差(係数の変化量を返す)
 * @return 0:成功，1:収束しなかった
 */
int m_estimation(const vector<double> &fi, const vector<double> &xi, int m, int n, double wfunc(double, double), double k,
				 vector<double> &c, vector<double> &w, int &max_iter, double &eps)
{
	vector<double> r, ar(m), c_old;

	// 初期値 : 通常の最小二乗法
	w.assign(m, 1.0);
	weighted_leastsquare(fi, xi, w, m, n, c);

	int ret = 1;
	int iter;
	double e = 0.0;
	for(iter = 0; iter < max_iter; ++iter){
		// 残差とスケールσ(MAD)
		residuals(fi, xi, m, c, n, r);
		double med = median(r);
		for(int i = 0; i < m; ++i) ar[i] = fabs(r[i]-med);
		double sigma = median(ar)/0.6745;
		if(sigma < 1e-12) { ret = 0; break; } // 残差がほぼ0(全データにぴったり合っている)

		// 重みを更新して重み付き最小二乗法
		for(int i = 0; i < m; ++i) w[i] = wfunc(r[i]/sigma, k);
		c_old = c;
		weighted_leastsquare(fi, xi, w, m, n, c);

		// 係数の変化量で収束判定
		e = 0.0;
		for(int i = 0; i <= n; ++i) e += (c[i]-c_old[i])*(c[i]-c_old[i]);
		e = sqrt(e);
		if(e < eps){ ret = 0; iter++; break; }
	}

	max_iter = iter;
	eps = e;
	return ret;
}


//-----------------------------------------------------------------------------
// データの作成と結果の表示
//-----------------------------------------------------------------------------
/*!
 * 外れ値を含むサンプリング点(データ点)の生成(1次元)
 *  - 等間隔の点で関数値に[-0.5nwidth, 0.5nwidth]の一様ノイズを加え，
 *    そのうち num_outlier 個のデータに [owidth/2, owidth] の大きさの外れ値(符号はランダム)を加える
 *  - 結果を再現できるように乱数のシードを指定する
 * @param[in] x0,x1 サンプリング範囲
 * @param[in] m データ数
 * @param[in] func 関数値を与える関数ポインタ
 * @param[in] nwidth ノイズの幅
 * @param[in] num_outlier 外れ値の数
 * @param[in] owidth 外れ値の大きさ
 * @param[out] xi,yi サンプリングデータ
 * @param[out] outlier 外れ値を加えたデータなら1
 * @param[in] seed 乱数のシード
 */
void MakeSamplingPointsWithOutliers(double x0, double x1, int m, double func(double), double nwidth, int num_outlier, double owidth,
									vector<double> &xi, vector<double> &yi, vector<int> &outlier, unsigned int seed = 12345)
{
	std::mt19937 rng(seed);
	std::uniform_real_distribution<double> uni(0.0, 1.0);
	xi.resize(m); yi.resize(m); outlier.assign(m, 0);
	for(int k = 0; k < m; ++k){
		xi[k] = x0+(x1-x0)*k/(m-1.0);
		yi[k] = func(xi[k])+(uni(rng)-0.5)*nwidth;
	}
	vector<int> idx;
	random_sample(m, num_outlier, rng, idx);
	for(int k : idx){
		double o = owidth*(0.5+0.5*uni(rng));
		yi[k] += (uni(rng) < 0.5 ? -o : o);
		outlier[k] = 1;
	}
}

/*!
 * フィッティング結果の表示
 *  - 係数，位置xでの値と真値との誤差，[x0,x1]での真の関数との最大誤差
 */
void PrintResult(string name, const vector<double> &c, int n, double func(double), double x, double x0, double x1)
{
	double fx = polynominal1d(x, c, n);
	double emax = 0.0;
	for(int i = 0; i <= 200; ++i){
		double xx = x0+(x1-x0)*i/200.0;
		emax = max(emax, fabs(polynominal1d(xx, c, n)-func(xx)));
	}
	cout << name << endl;
	cout << "  c = ";
	for(int i = 0; i <= n; ++i) cout << c[i] << (i == n ? "" : ", ");
	cout << endl;
	cout << "  f(" << x << ") = " << fx << ",  error = " << fabs(fx-func(x)) << ",  max error in [" << x0 << "," << x1 << "] = " << emax << endl;
}

/*!
 * インライア/外れ値の判定結果の表示(外れ値と判定されたデータの番号)
 */
void PrintInlier(const vector<int> &inlier, const vector<int> &outlier)
{
	int m = (int)inlier.size(), wrong = 0;
	cout << "  rejected : ";
	for(int k = 0; k < m; ++k){
		if(!inlier[k]) cout << k << " ";
		if(inlier[k] == outlier[k]) wrong++; // inlier=1かつoutlier=1，またはinlier=0かつoutlier=0なら判定ミス
	}
	cout << " (misclassified : " << wrong << ")" << endl;
}


//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	double(*func)(double) = FuncExp;
	double x0 = 0, x1 = 1;

	int d = 3;			// フィッティングする多項式の次数
	int m = 21;			// データ数
	int num_outlier = 5;// 外れ値の数
	double x = 0.5;		// 値を調べる位置
	cout.precision(6);

	// 外れ値を含むデータ(ノイズ幅0.05，外れ値の大きさ0.5～1.0)
	vector<double> xi, yi;
	vector<int> outlier;
	MakeSamplingPointsWithOutliers(x0, x1, m, func, 0.05, num_outlier, 1.0, xi, yi, outlier);
	OutputSamplingPoints(xi, yi); // サンプリング点の画面出力
	cout << "outliers : ";
	for(int k = 0; k < m; ++k) if(outlier[k]) cout << k << " ";
	cout << endl << "ground truth f(" << x << ") = " << func(x) << endl << endl;

	vector<double> c, w;
	vector<int> inlier;

	// 通常の最小二乗法(比較用)
	weighted_leastsquare(yi, xi, vector<double>(), m, d, c);
	PrintResult("[least squares]", c, d, func, x, x0, x1);

	// 最小メディアン法
	double med = lmeds(yi, xi, m, d, 500, c, inlier);
	PrintResult("[LMedS] (median of r^2 = "+TOSTR(med)+")", c, d, func, x, x0, x1);
	PrintInlier(inlier, outlier);

	// RANSAC(閾値はノイズの大きさから決める)
	int cnt = ransac(yi, xi, m, d, 500, 0.1, c, inlier);
	PrintResult("[RANSAC] (inliers = "+TOSTR(cnt)+")", c, d, func, x, x0, x1);
	PrintInlier(inlier, outlier);

	// M推定(Huber)
	int max_iter = 100;
	double eps = 1e-8;
	m_estimation(yi, xi, m, d, weight_huber, 1.345, c, w, max_iter, eps);
	PrintResult("[M-estimation (Huber)] (iter = "+TOSTR(max_iter)+")", c, d, func, x, x0, x1);

	// M推定(Tukey biweight)
	max_iter = 100; eps = 1e-8;
	m_estimation(yi, xi, m, d, weight_tukey, 4.685, c, w, max_iter, eps);
	PrintResult("[M-estimation (Tukey)] (iter = "+TOSTR(max_iter)+")", c, d, func, x, x0, x1);
	cout << "  weights : ";
	for(int k = 0; k < m; ++k) cout << w[k] << (k == m-1 ? "" : ", ");
	cout << endl;


	//// グラフ描画用にデータ出力
	//OutputSamplingPoints(xi, yi, "dat/robust_data.txt"); // サンプリング点のファイル出力
	//OutputFunction(x0, x1, (x1-x0)/200, func, "dat/robust_ground_truth.txt");
	//weighted_leastsquare(yi, xi, vector<double>(), m, d, c);
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, d), "dat/robust_ls.txt");
	//lmeds(yi, xi, m, d, 500, c, inlier);
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, d), "dat/robust_lmeds.txt");
	//ransac(yi, xi, m, d, 500, 0.1, c, inlier);
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, d), "dat/robust_ransac.txt");
	//max_iter = 100; eps = 1e-8;
	//m_estimation(yi, xi, m, d, weight_tukey, 4.685, c, w, max_iter, eps);
	//OutputFunction(x0, x1, (x1-x0)/200, std::bind(polynominal1d, std::placeholders::_1, c, d), "dat/robust_mest.txt");

	return 0;
}


