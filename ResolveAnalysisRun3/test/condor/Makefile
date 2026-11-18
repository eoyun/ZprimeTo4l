CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall $(shell root-config --cflags)
LDLIBS   := $(shell root-config --libs)

TARGETS  := fillHeepEffi plotHeepEffi plotCutflow plotCut plotCutMarginal plotFrac2ndtrk plotIsoDist runHeepEffi

all: $(TARGETS)

%: %.cc HeepEffiCore.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f $(TARGETS)
