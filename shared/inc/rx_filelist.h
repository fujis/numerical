/*! 
  @file rx_filelist.h
	
  @brief フォルダ内のファイルリストを取得する
 
  @author Makoto Fujisawa
  @date 2014-11
*/


#ifndef _RX_FILELIST_H_
#define _RX_FILELIST_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#ifdef WIN32
#include <cstdlib>
#include <windows.h>
#include <tchar.h>
 #include <shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#else
#include <sys/types.h>
#include <dirent.h>
#endif

#include <vector>
#include <string>


using namespace std;


#ifdef WIN32
//-----------------------------------------------------------------------------
// ファイルリスト取得関数(Win32api使用)
//-----------------------------------------------------------------------------

/*!
 * ワイド文字列(wstring)からマルチバイト文字列(string)へ変換
 *  - 用いる前に setlocale(LC_CTYPE, ""); としておくこと
 * @param[in] src ワイド文字列(wstring)
 * @return マルチバイト文字列(string)
 */
inline static std::string RX_W2S(const std::wstring &src)
{
    char *mbs = new char[src.length()*MB_CUR_MAX+1];
	size_t num;	// 変換された文字数
    wcstombs_s(&num, mbs, src.length()*MB_CUR_MAX+1, src.c_str(), _TRUNCATE);
    std::string dst = mbs;
    delete [] mbs;
    return dst;
}
 
/*!
 * マルチバイト文字列(string)からワイド文字列(wstring)へ変換
 *  - 用いる前に setlocale(LC_CTYPE, ""); としておくこと
 * @param[in] src マルチバイト文字列(string)
 * @return ワイド文字列(wstring)
 */
inline static std::wstring RX_S2W(const std::string &src)
{
    wchar_t *wcs = new wchar_t[src.length()+1];
	size_t num;	// 変換された文字数
    mbstowcs_s(&num, wcs, src.length()+1, src.c_str(), _TRUNCATE);
    std::wstring dst = wcs;
    delete [] wcs;
    return dst;
}
/*!
 * 再帰的に全ファイル(拡張子指定)を取り出す
 * @param[in] dpath フォルダパス
 * @param[out] paths 見つかったファイル一覧
 * @param[inout] d 現在の階層数
 * @param[in] n 最大階層数
 * @param[in] exts 拡張子指定
 */
static void SearchFiles(const std::string &dpath, std::vector<std::string> &paths, int d, const int n, 
                        const std::vector<std::string> &exts)
{
    HANDLE handle;
    WIN32_FIND_DATA fd;
 
    // search first file with the wildcard "*" to find the all type of file
    handle = FindFirstFile(RX_S2W(dpath+"\\*").c_str(), &fd);
 
    // if fail to find the file
    if(handle == INVALID_HANDLE_VALUE){
        return;
    }
 
    // search next files
    do{
        // file name
        std::string name = RX_W2S(static_cast<LPCTSTR>(fd.cFileName));
        std::string fpath = dpath+"\\"+name;
 
		if(name == "." || name == "..") continue;
 
		if((fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && (n == -1 || d < n)){
			// if the path is directory, recursively search in the directory
			SearchFiles(fpath, paths, d+1, n, exts);
		}
		else{
			vector<std::string>::const_iterator i;
			for(i = exts.begin(); i != exts.end(); ++i){
				if(fpath.find(*i, 0) != std::string::npos) break;
			}
			
			// store the file path if the extension was matched
			if(i != exts.end()){
				paths.push_back(fpath);
			}
		}
	}while(FindNextFile(handle, &fd));
 
	// terminate the search
	FindClose(handle);
}

static void SearchFiles(const std::string &dir, std::vector<std::string> &paths, 
						std::vector<std::string> exts, const int n = 0)
{
	if(PathIsDirectory(RX_S2W(dir).c_str())){	// dirがディレクトリであるかどうかのチェック
		int d = 0;
		SearchFiles(dir, paths, d, n, exts);
	}
}

//! 画像リストの取得
inline static std::vector<std::string> GetImageFileList(std::string dir)
{
	std::vector<std::string> exts;
	exts.push_back("bmp");
	exts.push_back("BMP");
	exts.push_back("jpg");
	exts.push_back("JPG");
	exts.push_back("png");
	exts.push_back("PNG");

	std::vector<std::string> files;

	SearchFiles(dir, files, exts, 0);

	files.push_back("");
	return files;
}

//! ファイルリストの取得
inline static std::vector<std::string> GetFileList(std::string dir, std::string ext)
{
	std::vector<std::string> exts;
	exts.push_back(ext);

	std::vector<std::string> files;
	SearchFiles(dir, files, exts, 0);

	return files;
}

#else

using namespace std;

inline static void SearchFiles(string path, vector<string> &filelist)
{
	DIR* dp = opendir(path.c_str()); // ディレクトリを開く
	if(dp != NULL){
		// ディレクトリ内のファイルを調べて、一覧に加えていく(ファイルorフォルダの判定などはしていない）
		dirent* dent;
		do{
			dent = readdir(dp);
			if(dent != NULL){
				string name = dent->d_name;
				if(name == "." || name == "..") continue;
				filelist.push_back(path+"/"+name);
			}	
		}while(dent != NULL);
		closedir(dp);
	}

  return;
}


//! マッチング画像リストの取得
inline static std::vector<std::string> GetImageFileList(std::string dir)
{
	std::vector<std::string> files;
	SearchFiles(dir, files);
	//files.push_back(dir+"/box.png");
	//files.push_back(dir+"/object.png");
	//files.push_back("");
	return files;
}

//! ファイルリストの取得
inline static std::vector<std::string> GetFileList(std::string dir)
{
	std::vector<std::string> files;
	SearchFiles(dir, files);
	return files;
}


#endif




#endif // #ifndef _RX_FILELIST_H_