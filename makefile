# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c11 -g -fsanitize=address -DDEBUG -Iinclude
LDFLAGS = -L$(LIB_DIR) -lmylib  # Linker flags to specify libraries

# Directories
SRC_DIR = src
BIN_DIR = bin

# Target executable
TARGET = $(BIN_DIR)/test_app

# Source and object files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BIN_DIR)/%.o)

# Default target
all: $(TARGET)

# Linking step
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Compilation step for each .c file
$(BIN_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean the build
clean:
	rm -rf $(BIN_DIR)

# Phony targets
.PHONY: all clean
