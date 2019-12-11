/*! 
  @file jacobi.cpp
	
  @brief 固有値・固有ベクトル
         ヤコビ(Jacobi)法
 
  @author Makoto Fujisawa
  @date 2012-06
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"


/*!
 * Jacobi法による固有値の算出
 * @param[inout] a 実対称行列．計算後，対角要素に固有値が入る
 * @param[out] v 固有ベクトル(aと同じサイズ)
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 反復回数
 */
int eigen_jacobi(vector< vector<double> > &a, vector< vector<double> > &v, int n, int &max_iter, double &eps)
{
    vector<double> bim(n), bjm(n);
    double bii, bij, bjj, bji;
  
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            v[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
 
	double e;
    int k;
    for(k = 0; k < max_iter; ++k){
        int i, j;
 
		// 非対角要素で絶対値が最大の要素をa_ijとして0にする
        double x = 0.0;
        for(int ia = 0; ia < n; ++ia){
            for(int ja = 0; ja < n; ++ja){
                if(ia != ja && fabs(a[ia][ja]) > x){
                    i = ia;
                    j = ja;
                    x = fabs(a[ia][ja]);
                }
            }
        }
 
		// sinθ,cosθを算出
        double aii = a[i][i];
        double ajj = a[j][j];
        double aij = a[i][j];
 
        double alpha, beta;
        alpha = (aii-ajj)/2.0;
        beta  = sqrt(alpha*alpha+aij*aij);
 
        double st, ct;
        ct = sqrt((1.0+fabs(alpha)/beta)/2.0);    // sinθ
        st = (((aii-ajj) >= 0.0) ? 1.0 : -1.0)*aij/(2.0*beta*ct);    // cosθ
 
        // A = PAPの計算
        for(int m = 0; m < n; ++m){
            if(m == i || m == j) continue;
 
            double aim = a[i][m];
            double ajm = a[j][m];
 
            bim[m] =  aim*ct+ajm*st;
            bjm[m] = -aim*st+ajm*ct;
        }
 
        bii = aii*ct*ct+2.0*aij*ct*st+ajj*st*st;
        bij = 0.0;
 
        bjj = aii*st*st-2.0*aij*ct*st+ajj*ct*ct;
        bji = 0.0;
 
        for(int m = 0; m < n; ++m){
            a[i][m] = a[m][i] = bim[m];
            a[j][m] = a[m][j] = bjm[m];
        }
        a[i][i] = bii;
        a[i][j] = bij;
        a[j][j] = bjj;
        a[j][i] = bji;
 
        // V = PVの計算
        for(int m = 0; m < n; ++m){
            double vmi = v[m][i];
            double vmj = v[m][j];
 
            bim[m] =  vmi*ct+vmj*st;
            bjm[m] = -vmi*st+vmj*ct;
        }
        for(int m = 0; m < n; ++m){
            v[m][i] = bim[m];
            v[m][j] = bjm[m];
        }
 
		// 非対角要素の絶対値の和で収束を判定
        e = 0.0;
        for(int ja = 0; ja < n; ++ja){
            for(int ia = 0; ia < n; ++ia){
                if(ia != ja){
                    e += fabs(a[ja][ia]);
                }
            }
        }
        if(e < eps){
			k++;
			break;
		}
    }

	max_iter = k;
	eps = e;

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

	vector< vector<double> > V(n, vector<double>(n, 0.0));

	int max_iter = 100;
	double eps = 1e-6;
	eigen_jacobi(A, V, n, max_iter, eps);

	// 固有値の表示
	cout << "e = (";
	for(int i = 0; i < n; ++i){
		cout << A[i][i] << (i == n-1 ? ")" : ", ");
	}
	cout << endl;

	// 固有ベクトルの表示
	for(int j = 0; j < n; ++j){
		cout << "v" << j << " = (";
		for(int i = 0; i < n; ++i){
			cout << V[i][j] << (i == n-1 ? ")" : ", ");
		}
		cout << endl;
	}
	cout << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


