###########################################################
#
#
#			Makefile for the math library
#
#
###########################################################

# compiler and flags
CC = gcc
CFLAGS = -Wall -g
LCUNIT = -lcunit

# Directories
DIR_SEQUENTIAL = sequential
DIR_MULTITHREAD = multithread
DIR_TESTS = tests/sequential
DIR_UTILS = utils
DIR_OBJECT = object_math_lib
DIR_OBJECT_SEQUENTIAL = $(DIR_OBJECT)/$(DIR_SEQUENTIAL)
DIR_OBJECT_MULTITHREAD = $(DIR_OBJECT)/$(DIR_MULTITHREAD)
DIR_OBJECT_TESTS = $(DIR_OBJECT)/$(DIR_TESTS)
DIR_OBJECT_UTILS = $(DIR_OBJECT)/$(DIR_UTILS)
DIR_HEADERS = headers

# Cunit flags for tests
CUNIT_PREFIX = $(shell brew --prefix cunit)
LDFLAGS = -L$(CUNIT_PREFIX)/lib
INCLUDES = -I$(CUNIT_PREFIX)/include -I$(DIR_HEADERS)

# all C file from seq, mult and tests
SEQUENTIAL = $(wildcard $(DIR_SEQUENTIAL)/*.c)
MULTITHREAD = $(wildcard $(DIR_MULTITHREAD)/*.c)
TESTS = $(wildcard $(DIR_TESTS)/*.c)
UTILS = $(wildcard $(DIR_UTILS)/*.c)

# object of seq and mult
OBJS_SEQUENTIAL = $(patsubst $(DIR_SEQUENTIAL)/%.c, $(DIR_OBJECT_SEQUENTIAL)/%.o, $(SEQUENTIAL))
OBJS_MULTITHREAD = $(patsubst $(DIR_MULTITHREAD)/%.c, $(DIR_OBJECT_MULTITHREAD)/%.o, $(MULTITHREAD))
OBJS_UTILS = $(patsubst $(DIR_UTILS)/%.c, $(DIR_OBJECT_UTILS)/%.o, $(UTILS))
BINS_TESTS = $(patsubst $(DIR_TESTS)/%.c, $(DIR_OBJECT_TESTS)/%, $(TESTS))

# target result (exec)
TARGET_MAIN = $(DIR_OBJECT)/main



all: $(TARGET_MAIN)
tests: $(BINS_TESTS)
# run all tests
run-tests: tests
	@for t in $(BINS_TESTS); do echo "Running $$t:"; ./$$t; echo ""; done



# Compile sequential sources
$(DIR_OBJECT_SEQUENTIAL)/%.o: $(DIR_SEQUENTIAL)/%.c | $(DIR_OBJECT_SEQUENTIAL)
	$(CC) $(CFLAGS) -I$(DIR_HEADERS) -c $< -o $@

# Compile multithread sources
$(DIR_OBJECT_MULTITHREAD)/%.o: $(DIR_MULTITHREAD)/%.c | $(DIR_OBJECT_MULTITHREAD)
	$(CC) $(CFLAGS) -I$(DIR_HEADERS) -c $< -o $@

# Compile utils sources
$(DIR_OBJECT_UTILS)/%.o: $(DIR_UTILS)/%.c | $(DIR_OBJECT_UTILS)
	$(CC) $(CFLAGS) -I$(DIR_HEADERS) -c $< -o $@

# Rule to build each test binary
$(DIR_OBJECT_TESTS)/%: $(DIR_TESTS)/%.c $(OBJS_SEQUENTIAL) $(OBJS_MULTITHREAD) $(OBJS_UTILS)| $(DIR_OBJECT_TESTS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS_SEQUENTIAL) $(OBJS_MULTITHREAD) $(OBJS_UTILS) $(LDFLAGS) $(LCUNIT) -o $@

# Link all object files with main.c
$(TARGET_MAIN): main.c $(OBJS_SEQUENTIAL) $(OBJS_MULTITHREAD) $(OBJS_UTILS) | $(DIR_OBJECT)
	$(CC) $(CFLAGS) -I$(DIR_HEADERS) $^ -o $@



# Create required object directories
$(DIR_OBJECT):
	mkdir -p $@

$(DIR_OBJECT_SEQUENTIAL):
	mkdir -p $@

$(DIR_OBJECT_MULTITHREAD):
	mkdir -p $@

$(DIR_OBJECT_UTILS):
	mkdir -p $@

$(DIR_OBJECT_TESTS):
	mkdir -p $@


clean:
	rm -rf $(DIR_OBJECT)
