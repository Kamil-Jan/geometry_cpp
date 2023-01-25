SRC_DIR := src
INC_DIR := include
OBJ_DIR := obj

CC := g++
INCLUDES := -I$(INC_DIR)
CFLAGS := -std=c++20 -Wshadow -Wall

SRCS := $(shell find $(SRC_DIR) -name *.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

TARGET := tests
all: release

debug: CFLAGS += -DDEBUG -g -fsanitize=address -fsanitize=undefined
debug: $(TARGET)

release: CFLAGS += -O2
release: $(TARGET)

clean:
	rm -r $(OBJ_DIR)

$(TARGET): $(OBJS)
	@echo "--- Linking target $@"
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo "--- Building $@"
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
	@echo ""
