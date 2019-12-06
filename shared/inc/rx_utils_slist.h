/*! @file rx_utils_slist.h

	@brief 数値計算テストの共通ヘッダ

	@author Makoto Fujisawa
	@date   2012
*/


#ifndef _RX_UTILS_SLIST_H_
#define _RX_UTILS_SLIST_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_utils.h"

//! 学籍番号とレポート用数値1,数値2
struct rxSList
{
	int student_number;
	double x1;
	double x2;
};


/*!
 * テキストファイルから学籍番号と数値を読み込む
 * @param[in] file_name ファイル名
 * @param[in] sep 区切り文字
 * @param[out] mat 行列要素
 * @return ファイルが開けなかったらfalseを返す
 */
bool ReadSList(string file_name, string sep, vector<rxSList> &slist)
{
	ifstream file;

	file.open(file_name.c_str());
	if(!file || !file.is_open() || file.bad() || file.fail()){
		cout << "ReadSList : Invalid file specified" << endl;
		return false;
	}

	string buf, grid;
	string::size_type comment_start = 0;
	while(!file.eof()){
		getline(file, buf);

		// '#'以降はコメントとして無視
		if((comment_start = buf.find('#')) != string::size_type(-1))
			buf = buf.substr(0, comment_start);

		// 空行は無視
		if(buf.empty())
			continue;

		vector<string> nums;
		size_t pos = 0;
		do{
			string sub;
			pos = GetNextString(buf, sub, sep, pos);
			if(IsNumeric(sub)){
				nums.push_back(sub);
			}
		} while(pos != string::npos);

		if(nums.size() == 3){
			rxSList s;
			s.student_number = atoi(nums[0].c_str());
			s.x1 = atof(nums[1].c_str());
			s.x2 = atof(nums[2].c_str());
			slist.push_back(s);
		}
	}

	file.close();

	return true;
}

#endif