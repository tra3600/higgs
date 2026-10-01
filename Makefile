CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

higgs: higgs.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test: higgs
	./higgs test

clean:
	rm -f higgs

.PHONY: test clean
