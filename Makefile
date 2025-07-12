# Compiler and flags
CC ?= gcc
CFLAGS ?= -Wall -Wextra -std=c99 -O3 -g -pedantic -Werror
LFLAGS ?= -lz

SRC_DIR := src
BUILD_DIR := build
INSTALL_DIR ?= /usr/local/bin

# Target object file for main.c
MAIN_C := $(SRC_DIR)/main.c
OBJ := $(BUILD_DIR)/main.o

# Default target builds the executable in build/
all: fastx2bam

fastx2bam: $(OBJ)
	$(CC) $(CFLAGS) $(LFLAGS) -o fastx2bam $(OBJ)

# Compile main.c into object file in build/
$(BUILD_DIR)/main.o: $(MAIN_C)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

install: fastx2bam
	@echo "Installing fastx2bam to $(INSTALL_DIR)"
	install -m 755 fastx2bam $(INSTALL_DIR)

clean:
	rm -f fastx2bam $(BUILD_DIR)/*.o

.PHONY: all clean install
