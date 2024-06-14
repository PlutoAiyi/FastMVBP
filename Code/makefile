# Makefile

CXX = g++

CXXFLAGS = -std=c++20 -O3

TARGET = FastMVBP

SRC = main.cpp

DEPS = Framework.hpp Graph.hpp RandList.hpp MBitSet.hpp

all: $(TARGET)

$(TARGET): $(SRC) $(DEPS)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all clean
