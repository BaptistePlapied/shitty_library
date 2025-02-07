# Compiler and flags
CC = gcc
CFLAGS = -Wall -g `pkg-config --cflags sdl2`
LDFLAGS = `pkg-config --libs sdl2`

# Directories
SRC_DIR = src
BUILD_DIR = build
BUILD_DIR_BIN = $(BUILD_DIR)/bin
BUILD_DIR_DEBUG = $(BUILD_DIR)/debug
INC_DIR = include

# Automatically find all .c files in the src directory
SRCS = $(wildcard $(SRC_DIR)/*.c)

# Object and assembly files
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))
ASMS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR_BIN)/%.s, $(SRCS))

# Output executable
TARGET = $(BUILD_DIR)/program

# Default target (build the program and generate assembly files)
all: $(TARGET) $(ASMS)

# Link object files into an executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile .c files into .o files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Generate assembly (.s) files in the build/bin directory
$(BUILD_DIR_BIN)/%.s: $(SRC_DIR)/%.c | $(BUILD_DIR_BIN)
	$(CC) $(CFLAGS) -S -I$(INC_DIR) $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Create the build/bin directory if it doesn't exist
$(BUILD_DIR_BIN):
	mkdir -p $(BUILD_DIR_BIN)

$(BUILD_DIR_DEBUG):
	mkdir -p $(BUILD_DIR_DEBUG)

# Create a tags for VIM
tags:
	ctags -R .

# Run the compiled program
run: $(TARGET)
	./$(TARGET)

# Clean build artifacts (but keep the build directories)
clear:
	rm -f $(BUILD_DIR)/*.o $(TARGET)
	rm -f $(BUILD_DIR_BIN)/*.s
	rm -f $(BUILD_DIR_DEBUG)/*.o

#make clangd understand the path file
bear:
	bear -- make

# Phony targets (prevent conflicts with actual filenames)
.PHONY: all tags run clear bear

