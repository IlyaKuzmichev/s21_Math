CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -Werror -c

TESTFLAGS := -lcheck -lpthread -pthread

ifeq ($(shell uname -s),Linux)
	TESTFLAGS += -lsubunit -lm -lrt
endif

GCOVFLAGS = -fprofile-arcs -ftest-coverage

HEADERS=s21_math.h s21_utils_math.h
SOURCES= s21_math.c
OBJECTS=$(SOURCES:.c=.o)
TESTS_SOURCES=s21_math_test.c

GCDA_FILES=*.gcda
GCNO_FILES=*.gcno
GCOV_FILES=*.gcov

FILES_TO_CHECK=$(SOURCES) $(HEADERS) $(TESTS_SOURCES)


.PHONY: test rebuild clean stylecheck fmt
.SECONDARY: $(OBJECTS)


all: s21_math.a


%.out: %.o $(OBJECTS)
	$(CC) $^ -o $@
	
%.o: %.c
	$(CC) $(CFLAGS) $^ -o $@ 

s21_math.a:  $(OBJECTS)
	ar rcs s21_math.a $^
	ranlib s21_math.a

gcov_report: test
	lcov -t test.out -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

report_nolcov:
	gcov *gcno > coverage.txt
	open coverage.txt

test: test_s21_math.a
	$(CC) $(TESTS_SOURCES) $(TESTFLAGS) $(GCOVFLAGS) -L. s21_math.a -o test.out
	./test.out

test_s21_math.a: all_functions
	ar rcs s21_math.a *.o
	ranlib s21_math.a

all_functions: $(SOURCES)
	$(CC) $(CFLAGS) $^ $(GCOVFLAGS) -g

rebuild: clean all

clean: 
	@rm -rf  $(OBJECTS) $(GCDA_FILES) $(GCNO_FILES) $(GCOV_FILES) 
	@rm -rf report/ *.a *.out *.info *.o *.txt
.clang-format:
	@cp ../materials/linters/.clang-format .
	
stylecheck: .clang-format
	@clang-format -n $(FILES_TO_CHECK)
	
fmt: .clang-format
	@clang-format -i $(FILES_TO_CHECK)

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem $(FILES_TO_CHECK)

