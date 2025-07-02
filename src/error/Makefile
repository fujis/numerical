# コンパイラ
COMPILER = g++
CXXFLAGS = -O3 -std=c++11

# ライブラリ関係
LDFLAGS = 
LIBS    = 

# インクルードフォルダ
INCLUDE = -I../../shared/inc

# 生成バイナリファイル名＆ディレクトリ
TARGETS = error
TARGETDIR = ../../bin

# ソースファイルの場所(サブフォルダも走査して*.cppファイルをリスト化)
SRCROOT   = .
SRCDIRS  := $(shell find $(SRCROOT) -type d)
SOURCES   = $(foreach dir, $(SRCDIRS), $(wildcard $(dir)/*.cpp))

# 中間ファイル(*.o)を置く場所&ファイル名(cppファイルから決定)
OBJROOT   = .
OBJECTS   = $(addprefix $(OBJROOT)/, $(SOURCES:.cpp=.o)) 
OBJDIRS   = $(addprefix $(OBJROOT)/, $(SRCDIRS)) 

# 実行ファイルの作成
$(TARGETS): $(OBJECTS) $(LIBS)
	$(COMPILER) -o $(TARGETDIR)/$@ $^ $(LDFLAGS)

# 中間バイナリのディレクトリを掘りながら.cppを中間ファイル.oに
$(OBJROOT)/%.o: $(SRCROOT)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(COMPILER) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

run: $(TARGETS)
	cd $(TARGETDIR); ./$(TARGETS); cd -

clean:
	rm -f $(OBJECTS) $(TARGETDIR)/$(TARGETS)