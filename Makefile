EXE = SCarborSNV

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CXX = g++
CXXFLAGS = -Wall -g -std=c++11
PPFLAGS = -Iinclude
LDFLAGS = -Llib
LDLIBS = -lhts -lpthread

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(PPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ) $(EXE)

.PHONY: all clean
