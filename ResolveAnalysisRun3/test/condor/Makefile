CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall $(shell root-config --cflags)
LDLIBS   := $(shell root-config --libs)

runHeepEffi: runHeepEffi.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f runHeepEffi
