# Compiler and flags
CC = gcc
# CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c11 -g -fsanitize=address -DDEBUG -Iinclude  
# CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L  -g -fsanitize=address -Iinclude
# CFLAGS = -pg -Iinclude # the run program to generate the gmon.out and open it with command gprof program gmon.out > report.txt
 CFLAGS = -O3 -march=native -ftree-vectorize -fopt-info-vec -Wall -Wextra -fstack-protector-strong -D_FORTIFY_SOURCE=2 -Iinclude 
# LIB_DIR = 
LDFLAGS = -lm  
# LDFLAGS = -L$(LIB_DIR) -lmylib  # Linker flags to specify libraries

# Directories
SRC_DIR = src
BIN_DIR = bin

# Target executable
TARGET = $(BIN_DIR)/EXEC_CG

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
