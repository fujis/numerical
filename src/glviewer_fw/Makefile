# コンパイラ
COMPILER = g++
CXXFLAGS = -O3 -std=c++11

# ライブラリ関係
LDFLAGS = -lglfw -lGLEW -framework OpenGL
LIBS    = 
LIBDIR  = -L../../shared/lib

# インクルードフォルダ
INCLUDE = -I../../shared/inc -I./imgui

# 生成バイナリファイル名＆ディレクトリ
TARGETS = glviewer
TARGETDIR = ./bin

# ソースファイルの場所(サブフォルダも走査して*.cppファイルをリスト化)
#  - findのオプション : パス(-path)が"./inc"でない(-prune)場合(-o)にディレクトリを表示(-type d -print)
SRCROOT   = .
SRCDIRS  := $(shell find $(SRCROOT) -path "./inc" -prune -o -type d -print)
SOURCES   = $(foreach dir, $(SRCDIRS), $(wildcard $(dir)/*.cpp))

# 中間ファイル(*.o)を置く場所&ファイル名(cppファイルから決定)
OBJROOT   = .
OBJECTS   = $(addprefix $(OBJROOT)/, $(SOURCES:.cpp=.o)) 
OBJDIRS   = $(addprefix $(OBJROOT)/, $(SRCDIRS)) 

# 実行ファイルの作成
$(TARGETS): $(OBJECTS) $(LIBS)
	$(COMPILER) -o $(TARGETDIR)/$@ $^ $(LIBDIR) $(LDFLAGS)

# 中間バイナリのディレクトリ内を走査して.cppを中間ファイル.oに
$(OBJROOT)/%.o: $(SRCROOT)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(COMPILER) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

run: $(TARGETS)
	cd $(TARGETDIR); ./$(TARGETS); cd -

clean:
	rm -f $(OBJECTS) $(TARGETDIR)/$(TARGETS)