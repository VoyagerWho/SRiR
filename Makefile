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
	upcxx-run -shared-heap 512M -n 16 -N 64 $$(upcxx-nodes nodes) ./bin/main -m 1600 -f .txt --compareX
	
.PHONY: nodes
nodes:
	/opt/nfs/config/station204_name_list.sh 1 16 > nodes

.PHONY: source
source:
	source /opt/nfs/config/source_upcxx_2023.3.sh

