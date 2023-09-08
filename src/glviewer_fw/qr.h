/*! 
  @file qr.h

  @brief 固有値・固有ベクトル
		 QR法

  @author Makoto Fujisawa
  @date 2019-12
*/

#ifndef _QR_H_
#define _QR_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"

//-----------------------------------------------------------------------------
// 関数の宣言
//-----------------------------------------------------------------------------
void qr_decomposition(const vector< vector<double> > &A, vector< vector<double> > &Q, vector< vector<double> > &R, int n);
int inverse_iteration(const vector< vector<double> > &A, const vector<double> lambda, vector< vector<double> > &v, int n, int &max_iter, double &eps);
int qr(vector< vector<double> > &A, int n, int &max_iter, double &eps);
void covariance_matrix2(const vector<rxPoint2> &data, vector< vector<double> > &A, rxPoint2 &c);



#endif // #ifndef _QR_H_