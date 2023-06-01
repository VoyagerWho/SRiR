CXX = UPCXX_GASNET_CONDUIT=udp upcxx
CXXFLAGS = -Wall -Wextra -Wpedantic -O2 -std=c++11

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS := $(wildcard $(SRC_DIR)/*.cpp) main.cpp
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(filter-out main.cpp,$(SRCS))) $(OBJ_DIR)/main.o
DEPS := $(OBJS:%.o=%.d)

$(BIN_DIR)/main: $(OBJS) | $(BIN_DIR)	
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/main.o: main.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR)
	
.PHONY: run
run:
	upcxx-run -shared-heap 512M -n 8 -N 8 -ssh-servers stud204-01,stud204-02,stud204-05,stud204-07,stud204-09,stud204-11,stud204-12,stud204-14  ./bin/main -m 8 -f 3D.txt --compareX
	
.PHONY: nodes
nodes:
	/opt/nfs/config/station204_name_list.sh 1 16 > nodes

.PHONY: source
source:
	source /opt/nfs/config/source_upcxx.sh

