/*! 
  @file qr.cpp
	
  @brief 固有値・固有ベクトル
         QR法
 
  @author Makoto Fujisawa
  @date 2012-06
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"

/*!
 * ハウスホルダー変換で実対称行列をヘッセンベルグ行列に変換
 * @param[inout] a 元の行列(n×n)．変換された三重対角行列を格納
 * @param[in] n 行列のサイズ
 */
void HouseholderTransformationForQR(vector< vector<double> > &a, int n)
{
	vector< vector<double> > b(n, vector<double>(n, 0.0));
	vector< vector<double> > p(n, vector<double>(n, 0.0));
	vector< vector<double> > q(n, vector<double>(n, 0.0));
	vector<double> u(n, 0.0);
 
    for(int k = 0; k < n-2; ++k){
        // sの計算
        float s = 0;
        for(int i = k+1; i < n; ++i){
            s += a[i][k]*a[i][k];
        }
        s = -((a[k+1][k] >= 0) ? 1 : -1)*sqrt(s);
 
        // |x-y|の計算
        float alpha = sqrt(2.0*s*(s-a[k+1][k]));
        if(fabs(alpha) < 1e-8) continue;
 
        // uの計算
        u[k+1] = (a[k+1][k]-s)/alpha;
        for(int i = k+2; i < n; ++i){
            u[i] = a[i][k]/alpha;
        }
 
        // Pの計算
        for(int i = k+1; i < n; ++i){
            for(int j = i; j < n; ++j){
                if(j == i){
                    p[i][j] = 1.0-2.0*u[i]*u[i];
                }
                else{
                    p[i][j] = -2.0*u[i]*u[j];
                    p[j][i] = p[i][j];
                }
            }
        }
 
        // PAの計算
        for(int i = k+1; i < n; ++i){
            for(int j = k+1; j < n; ++j){
                q[i][j] = 0.0;
                for(int m = k+1; m < n; ++m){
                    q[i][j] += p[i][m]*a[m][j];
                }
            }
        }
 
        // A = PAP^Tの計算
        for(int i = 0; i <= k; ++i){
            b[i][k] = a[i][k];
        }
        b[k+1][k] = s;
        for(int i = k+2; i < n; ++i){
            b[i][k] = 0.0;
        }
        for(int j = k+1; j < n; ++j){
            for(int i = 0; i <= k; ++i){
                b[i][j] = 0.0;
                for(int m = k+1; m < n; ++m){
                    b[i][j] += a[i][m]*p[j][m];
                }
            }
            for(int i = k+1; i < n; ++i){
                b[i][j] = 0.0;
                for(int m = k+1; m < n; ++m){
                    b[i][j] += q[i][m]*p[j][m];
                }
            }
        }
 
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				a[i][j] = b[i][j];
			}
        }
    }
}


/*!
 * QR法による固有値の算出
 * @param[inout] h 元の行列(ヘッセンベルグ行列)．計算後，対角要素に固有値が入る
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 反復回数
 */
int EigenQR(vector< vector<double> > &h, int n, int &max_iter, double &eps)
{
	vector< vector<double> > r(n, vector<double>(n, 0.0));
	vector< vector<double> > q(n, vector<double>(n, 0.0));
	vector< vector<double> > t(n, vector<double>(n, 0.0));
	vector<double> u(n, 0.0), v(n, 0.0);
 
    // R = H : ヘッセンベルグ行列
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            r[i][j] = h[i][j];
        }
    }
 
	double e;
    int l;
    for(l = 0; l < max_iter; ++l){
        // Q=I (単位行列)
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                q[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
 
        for(int k = 0; k < n-1; ++k){
            // sinθ,cosθの計算
            double alpha = sqrt(r[k][k]*r[k][k]+r[k+1][k]*r[k+1][k]);
            if(fabs(alpha) < 1e-8) continue;
 
            double c = r[k][k]/alpha;
            double s = -r[k+1][k]/alpha;
 
            // Rの計算
            for(int j = k+1; j < n; ++j){
                u[j] = c*r[k][j]-s*r[k+1][j];
                v[j] = s*r[k][j]+c*r[k+1][j];
            }
            r[k][k] = alpha;
            r[k+1][k] = 0.0;
            for(int j = k+1; j < n; ++j){
                r[k][j]   = u[j];
                r[k+1][j] = v[j];
            }
 
            // Qの計算
            for(int j = 0; j <= k; ++j){
                u[j] = c*q[k][j];
                v[j] = s*q[k][j];
            }
            q[k][k+1] = -s;
            q[k+1][k+1] = c;
            for(int j = 0; j <= k; ++j){
                q[k][j]   = u[j];
                q[k+1][j] = v[j];
            }
        }
 
        // RQの計算
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                double rq = 0.0;
                for(int m = 0; m < n; ++m){
                    rq += r[i][m]*q[j][m];
                }
 
                t[i][j] = rq;
            }
        }
 
        // 収束判定
        e = 0.0;
        for(int i = 0; i < n; ++i){
            e += fabs(t[i][i]-h[i][i]);
        }
        if(e < eps){
			l++;
			break;
		}
 
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                r[i][j] = t[i][j];
                h[i][j] = t[i][j];
            }
        }
    }

	max_iter = l;
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

	HouseholderTransformationForQR(A, n);

	int max_iter = 100;
	double eps = 1e-6;
	EigenQR(A, n, max_iter, eps);

	// 固有値の表示
	cout << "e = (";
	for(int i = 0; i < n; ++i){
		cout << A[i][i] << (i == n-1 ? ")" : ", ");
	}
	cout << endl;

	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	return 0;
}


