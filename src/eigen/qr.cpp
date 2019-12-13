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

/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分解する
 *  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
 *  - LとUを一つの行列にまとめた形で結果を返す
 * @param[inout] A n×nの係数行列．LU分解した結果を格納する．
 * @param[in] n 行列の大きさ
 * @return 0:成功,1:失敗
 */
int lu_decomp(vector< vector<double> > &A, int n)
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
int lu_solver(const vector< vector<double> > &A, const vector<double> &b, vector<double> &x, int n)
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


/*!
 * ハウスホルダー変換でQR分解 (A=QR)
 * @param[in] A 元の行列(n×n)
 * @param[in] Q 分解された直行行列(n×n)
 * @param[in] R 分解された上三角行列(n×n)
 * @param[in] n 行列のサイズ
 */
void qr_decomposition(const vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R, int n)
{
	vector<double> u(n, 0.0);
	vector< vector<double> > H(n, u); // 行列H
	R = A; // A^(0)
	Q.resize(n, u); unit(Q, n); // Qを単位行列で初期化

	for(int k = 0; k < n-1; ++k){
		// a_kk^(k+1)の計算 (ベクトルxの大きさ)
		double akk = 0.0;
		for(int i = k; i < n; ++i) akk += R[i][k]*R[i][k];
		akk = sqrt(akk);

		// u=(x1-y1)/|x1-y1|の計算
		double l = 0.0;
		for(int i = k; i < n; ++i){
			u[i] = R[i][k]-(i == k ? akk : 0.0);
			l += u[i]*u[i];
		}
		if(fabs(l) < 1e-10) break;
		l = sqrt(l);
		for(int i = k; i < n; ++i) u[i] /= l;
		
		// ハウスホルダー行列H^(k)の計算(H=I-2uu^T)
		unit(H, n); // Hを単位行列で初期化
		for(int i = k; i < n; ++i){
			for(int j = k; j < n; ++j){
				H[i][j] -= 2*u[i]*u[j];
			}
		}

		// A^(k+1) = H^(k) A^(k)
		R = mul_mm(H, R, n);

		// Q^(k+1) = Q^(k) (H^(k))^T
		Q = mul_mm(Q, transpose(H, n), n);
	}
}

/*!
 * 逆反復法による固有ベクトルの算出
 * @param[inout] A 元の行列
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 反復回数
 */
int inverse_iteration(const vector< vector<double> > &A, const vector<double> lambda, vector< vector<double> > &v, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;

	vector< vector<double> > LU(A); // LU分解した行列を格納するための配列
	vector<double> y(n, 1.0), y1(n); // y^(k), y^(k+1)格納用
	double sum_e = 0.0;
	int sum_k = 0;
	v.resize(n);

	for(int l = 0; l < n; ++l){
		// (A-λI)を計算
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				LU[i][j] = A[i][j]-(i == j ? lambda[l] : 0.0);
			}
		}

		// (A-λI)をLU分解
		lu_decomp(LU, n);

		for(int i = 0; i < n; ++i) y[i] = 1.0;
		double l2 = dot(y, y, n);
		double mu, mu0 = 0.0;
		double e;
		int k;
		for(k = 0; k < max_iter; ++k){
			// |y^(k)|=1となるように正規化
			double ly = sqrt(l2); // |y^(k)|の計算
			y = mul_sv(1.0/ly, y, n); // |y^(k)|で割る(1/|y|を掛ける)

			// y^(k+1) = B y^(k)の計算(By^(k+1)=y^(k)をLU分解で解く)
			lu_solver(LU, y, y1, n);

			// 固有値1/|λ-λi|の計算(|y^(k)|=1で正規化されていることが前提)
			mu = dot(y, y1, n);

			// 収束判定
			if((e = fabs((mu-mu0)/mu)) < eps) break;

			y = y1;
			l2 = dot(y, y, n);
			mu0 = mu;
		}
		v[l] = mul_sv(1.0/sqrt(dot(y, y, n)), y, n); // 最後に正規化しておく
		sum_e += e;
		sum_k += k;
	}


	max_iter = sum_k/n;
	eps = sum_e/n;


	return 0;
}

/*!
 * QR法による固有値の算出
 * @param[inout] A 元の行列．計算後，対角要素に固有値が入る
 * @param[in] n 行列のサイズ(n×n)
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す)
 * @return 反復回数
 */
int qr(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
	if(n <= 0) return 1;
	vector< vector<double> > R(A), Q(A); // R,Q用配列を確保(Aと同じ大きさに初期化)
	vector<double> lambda(n, 0.0); // 収束判定用
	double e = eps*n;
	int k;
	for(k = 0; k < max_iter; ++k){
		// A^(k) ⇐ Q^(k) R^(k)
		qr_decomposition(A, Q, R, n);

		// A^(k+1) ⇐ R^(k) Q^(k)
		A = mul_mm(R, Q, n);

		// 確認用表示
		cout << k << " : lambda = ";
		for(int i = 0; i < n; ++i) cout << A[i][i] << (i == n-1 ? "" : ", ");
		cout << endl;

		// 収束判定
		e = 0;
		for(int i = 0; i < n; ++i) e += fabs(lambda[i]-A[i][i]);
		if(e/n < eps) break;

		for(int i = 0; i < n; ++i) lambda[i] = A[i][i];
	}
	eps = e/n;
	max_iter = k;

	return 0;
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
	vector< vector<double> > A0(A);
	cout << endl;

	// QR分解
	vector< vector<double> > Q(A), R(A);
	qr_decomposition(A, Q, R, n);
	cout << "Q = " << endl;
	OutputMatrix(Q, n, n);
	cout << "R = " << endl;
	OutputMatrix(R, n, n);

	vector< vector<double> > QR = mul_mm(Q, R, n);
	cout << "QR = " << endl;
	OutputMatrix(QR, n, n);

	vector< vector<double> > QQ = mul_mm(Q, transpose(Q, n), n);
	cout << "QQ = " << endl;
	OutputMatrix(QQ, n, n);
	cout << endl;


	// QR法で固有値を求める
	int max_iter = 100;
	double eps = 1e-6;
	qr(A, n, max_iter, eps);

	// 固有値の表示
	cout << "e = (";
	vector<double> lambda(n);
	for(int i = 0; i < n; ++i){
		cout << A[i][i] << (i == n-1 ? ")" : ", ");
		lambda[i] = A[i][i];
	}
	cout << endl;
	cout << "iter = " << max_iter << ", eps = " << eps << endl;
	cout << endl;

	// 逆反復法で固有ベクトルを計算
	max_iter = 100;
	eps = 1e-6;
	vector< vector<double> > v;
	inverse_iteration(A0, lambda, v, n, max_iter, eps);

	// 固有ベクトルの表示
	for(int i = 0; i < n; ++i){
		cout << "v" << i+1 << " = (";
		for(int j = 0; j < n; ++j){
			cout << v[i][j] << (j == n-1 ? ")" : ", ");
		}
		cout << endl;
	}
	cout << "iter = " << max_iter << ", eps = " << eps << endl;


	cout << endl;



	//HouseholderTransformationForQR(A, n);

	//int max_iter = 100;
	//double eps = 1e-6;
	//EigenQR(A, n, max_iter, eps);

	//// 固有値の表示
	//cout << "e = (";
	//for(int i = 0; i < n; ++i){
	//	cout << A[i][i] << (i == n-1 ? ")" : ", ");
	//}
	//cout << endl;

	//cout << "iter = " << max_iter << ", eps = " << eps << endl;
	//cout << endl;

	return 0;
}


