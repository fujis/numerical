/*! 
  @file householder.chpp
	
  @brief 固有値・固有ベクトル
         ハウスホルダー(Householder)法
 
  @author Makoto Fujisawa
  @date 2012-06
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"

/*!
 * ハウスホルダー変換で実対称行列を三重対角行列に変換
 * @param[inout] a 元の行列(n×n)．変換された三重対角行列を格納
 * @param[in] n 行列のサイズ
 */
void HouseholderTransformation(vector< vector<double> > &a, int n)
{
	vector<double> u(n, 0.0), v(n, 0.0), q(n, 0.0);
 
    for(int k = 0; k < n-2; ++k){
        // sの計算
        double s = 0;
        for(int i = k+1; i < n; ++i){
            s += a[i][k]*a[i][k];
        }
        s = -((a[k+1][k] >= 0) ? 1 : -1)*sqrt(s);
 
        // |x-y|の計算
        double alpha = sqrt(2.0*s*(s-a[k+1][k]));
        if(fabs(alpha) < 1e-8) continue;
 
        // uの計算
        u[k+1] = (a[k+1][k]-s)/alpha;
        for(int i = k+2; i < n; ++i){
            u[i] = a[i][k]/alpha;
        }
 
        // Auの計算
        q[k] = alpha/2.0;
        for(int i = k+1; i < n; ++i){
            q[i] = 0.0;
            for(int j = k+1; j < n; ++j){
                q[i] += a[i][j]*u[j];
            }
        }
 
        // v=2(Au-uu^T(Au))の計算
        alpha = 0.0;
        for(int i = k+1; i < n; ++i){
            alpha += u[i]*q[i];
        }
        v[k] = 2.0*q[k];
        for(int i = k+1; i < n; ++i){
            v[i] = 2.0*(q[i]-alpha*u[i]);
        }
 
        // A = PAP = A-uv^T-vu^Tの計算
        a[k][k+1] = a[k+1][k] = s;
        for(int i = k+2; i < n; ++i){
            a[k][i] = 0.0;
            a[i][k] = 0.0;
        }
        for(int i = k+1; i < n; ++i){
            a[i][i] = a[i][i]-2*u[i]*v[i];
            for(int j = i+1; j < n; ++j){
                a[i][j] = a[i][j]-u[i]*v[j]-v[i]*u[j];
                a[j][i] = a[i][j];
            }
        }
    }
}


/*!
 * スツルム関数列の符号変化回数をカウント
 *  - λが大きいときのf_k(λ)の計算でオーバーフローが発生するのを防ぐために，
 *    f_k でなく g_k = f_k/f_(k-1) を使う．
 *  - g_k < 0 となる回数を数える
 * @param[in] lambda 固有値候補
 * @param[inout] a 三重対角行列(n×n)
 * @param[in] n 行列のサイズ
 * @return 符号変化回数
 */
int NCS(double lambda, const vector< vector<double> > &A, int n)
{
	double g, g0 = lambda-A[0][0];
	int ncs = (g0 < 0 ? 1 : 0);
	for(int k = 1; k < n; ++k){
		if(g0 == 0) g0 = RX_FEQ_EPS;

		double a = A[k][k];
		double b = A[k-1][k];
		g = lambda-a-(b*b)/g0;
		if(g < 0) ncs++;

		g0 = g;
	}

	return ncs;
}

/*!
 * 三重対角行列の特性方程式の漸化式(スツルム関数列)を計算
 *  - ニュートン法で使うため，勾配値も同時に計算
 * @param[in] lambda 固有値候補(変数)
 * @param[out] f スツルム関数列f_nの値
 * @param[out] f スツルム関数列f_nの勾配値
 * @param[inout] a 三重対角行列(n×n)
 * @param[in] n 行列のサイズ
 */
void CFunc(double lambda, double &f, double &df, const vector< vector<double> > &A, int n)
{
	if(n == 0){
		f = 1.0;
		df = 0.0;
		return;
	}
	if(n == 1){
		f = lambda-A[0][0];
		df = 1.0;
		return;
	}

	double f0  = 1.0;	// f_0
	double f1  = lambda-A[0][0];	// f_1
	double df0 = 0.0;
	double df1 = 1.0;

	for(int k = 1; k < n; ++k){
		double a = A[k][k];
		double b = A[k-1][k];
		f  = (lambda-a)*f1-b*b*f0;
		df = (lambda-a)*df1-b*b*df0+f1;

		f0 = f1;
		f1 = f;
		df0 = df1;
		df1 = df;
	}
}

/*!
 * ニュートン法でf_n = 0となるλを探索
 * @param[in] x 探索開始位置
 * @param[in] eps 収束判定用許容誤差
 * @param[inout] a 三重対角行列(n×n)
 * @param[in] n 行列のサイズ
 * @return 根
 */
double Newton(double x, double eps, const vector< vector<double> > &A, int n)
{
	double x0;
	int iter = 0;
	do{
		x0 = x;

		double f, df;
		CFunc(x0, f, df, A, n);

		x = x0-f/df;
		iter++;
	}while(fabs(x-x0) > eps && iter < 100);

	return x;
}



/*!
 * ハウスホルダー法による固有値の算出
 *  - 実対称行列を三重対角行列に変換し，
 * @param[inout] H 実対称行列．計算後，対角要素に固有値が入る
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] eps 許容誤差
 * @return 反復回数
 */
int EigenHouseholder(vector< vector<double> > &H, int n, double &eps)
{
	//
	// ハウスホルダー変換で三重対角行列に変換
	//
	HouseholderTransformation(H, n);

	cout << "H(" << n << " x " << n << ") = " << endl;
	OutputMatrix(H, n, n);
	cout << endl;


	//
	// 三重対角行列の固有値をスツルム(Sturm)の定理と2分法,ニュートン法で算出
	//

	// 三重対角行列の対角成分と非対角成分
	vector<double> a(n), b(n-1);
	for(int i = 0; i < n-1; ++i){
		a[i] = H[i][i];
		b[i] = H[i+1][i];
	}
	a[n-1] = H[n-1][n-1];

	// 探索範囲[alpha, beta]をゲルシュゴリン(Gerschgorin)の定理より算出
	double alpha = a[0]-fabs(b[0]);
	double beta  = a[0]+fabs(b[0]);
	for(int i = 1; i < n-1; ++i){
		alpha = RX_MIN(alpha, a[i]-fabs(b[i])-fabs(b[i-1]));
		beta  = RX_MAX(beta,  a[i]+fabs(b[i])+fabs(b[i-1]));
	}
	alpha = RX_MIN(alpha, a[n-1]-fabs(b[n-2]));
	beta  = RX_MAX(beta,  a[n-1]+fabs(b[n-2]));


//	// ファイル出力
//	ofstream fo;
//	fo.open("_debug.txt");
//
//	double x1 = alpha;
//	while(x1 < beta){
//		double f, df;
//		CFunc(x1, f, df, H, n);
//
//		fo << x1 << ", " << f << ", " << NCS(x1, H, n) << endl;
//
//		x1 += 0.01;
//	}
//	fo.close();
	
	double ermax = fabs(beta-alpha)/n;

	int m = 0;
	int p = 0;
	vector<double> xsl(n), xsr(n);
	vector<int> nsl(n), nsr(n);

	// 探索範囲境界でのスツルム関数列の符号変化回数を計算
	xsl[p] = alpha;
	nsl[p] = NCS(alpha, H, n);
	xsr[p] = beta;
	nsr[p] = NCS(beta, H, n);

	// 二分法で符号変化回数の差が1となる領域を探索
	double x, xl = xsl[p], xr = xsr[p], xm;
	int nl, nr, nm;
	vector<double> eigen(n);
	while(p >= 0){
		if( (nsl[p]-nsr[p] == 1) && (xsr[p]-xsl[p] < ermax) ){
			// 符号変化回数の差が1となったらニュートン法で固有値を探索
			x = 0.5*(xsl[p]+xsr[p]);
			eigen[m] = Newton(x, eps, H, n);
			p--;
			m++;

			if(p < 0) break;
		}

		xl = xsl[p];
		nl = nsl[p];
		xr = xsr[p];
		nr = nsr[p];
		xm = 0.5*(xl+xr);
		nm = NCS(xm, H, n);
		p--;

		// 符号変化回数に変化があれば，中点を新たな探索境界にする
		if(nl > nm){
			p++;
			xsl[p] = xl;
			nsl[p] = nl;
			xsr[p] = xm;
			nsr[p] = nm;
		}
		if(nm > nr){
			p++;
			xsl[p] = xm;
			nsl[p] = nm;
			xsr[p] = xr;
			nsr[p] = nr;
		}
	}

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			H[i][j] = (i == j ? eigen[i] : 0.0);
		}
	}

	return 1;
}



//-----------------------------------------------------------------------------
//! メイン関数
//-----------------------------------------------------------------------------
int main(void)
{
	// ファイルから行列要素を読み込んで解く
	vector< vector<double> > A;

	ReadMatrix("matrix.txt", ",", A);

	int n = (int)A.size();	// n×n行列
	cout << "A(" << n << " x " << n << ") = " << endl;
	OutputMatrix(A, n, n);
	cout << endl;

	// ハウスホルダー法で固有値を算出
	double eps = 1e-6;
	EigenHouseholder(A, n, eps);

	// 固有値の表示
	cout << "e = (";
	for(int i = 0; i < n; ++i){
		cout << A[i][i] << (i == n-1 ? ")" : ", ");
	}
	cout << endl;

	return 0;
}


