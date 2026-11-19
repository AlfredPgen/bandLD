# bandLD - build:  make          (dynamic)
#                  make static   (single portable binary)
# Debian/Ubuntu:   apt install g++ make libdeflate-dev libzstd-dev libsqlite3-dev
CXX      ?= g++
CXXFLAGS ?= -O3 -march=x86-64-v2 -mtune=generic -std=c++17 -fopenmp -Wall -Wno-unused-result -Wno-psabi
LIBS     ?= -lzstd -ldeflate -lsqlite3 -lpthread -ldl -lm

SRCS = src/main.cpp src/util.cpp src/bgen.cpp src/bandfile.cpp src/build.cpp src/query.cpp
HDRS = $(wildcard src/*.hpp)

all: bandLD

bandLD: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS) $(LIBS)

static: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -static -o bandLD $(SRCS) -Wl,--start-group $(LIBS) -Wl,--end-group

asan: $(SRCS) $(HDRS)
	$(CXX) -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=undefined \
	  -std=c++17 -fopenmp -march=x86-64-v2 -o bandLD-asan $(SRCS) $(LIBS)

test: bandLD
	bash tests/run_tests.sh ./bandLD

clean:
	rm -f bandLD bandLD-asan

.PHONY: all static asan test clean
