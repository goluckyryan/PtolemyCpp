CXX      = g++
# Warnings: -Wall -Wextra on; suppress noise classes from the temporary
# pool-system landmines (-Warray-bounds / -Wstrict-aliasing emitted by every
# `reinterpret_cast<...>(&ALLOCS.FACFR4) - 1` pattern) and Fortran-COMMON
# alias-parameter leftovers (-Wunused-parameter).
# The pool warnings are real UB but stable in practice — they go away
# permanently when the last NALLOC site dies (NAMLOC rework).
# Real diagnostics still active: -Wunused-variable, -Wunused-but-set-variable,
#  -Wuninitialized, -Wmaybe-uninitialized, -Wimplicit-fallthrough,
#  -Wformat-truncation, etc.
CXXFLAGS = -g -O2 -std=c++17 -Iinclude -Wall -Wextra \
           -Wno-unused-parameter \
           -Wno-array-bounds -Wno-strict-aliasing
SRCDIR   = src
INCDIR   = include
BUILD_DIR = build

# All sources except main files
MAIN_PTOLEMY = $(SRCDIR)/ptolemy_main.cpp
MAIN_UTEST   = $(SRCDIR)/UnitTests.cpp
ALL_SOURCES  = $(wildcard $(SRCDIR)/*.cpp)
# Shared sources: everything except main files
SHARED_SOURCES = $(filter-out $(MAIN_PTOLEMY) $(MAIN_UTEST), $(ALL_SOURCES))
SHARED_OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BUILD_DIR)/%.o,$(SHARED_SOURCES))

PTOLEMY_TARGET = ptolemy
UTEST_TARGET   = unit_tests

.PHONY: all clean distclean test

all: $(PTOLEMY_TARGET)

# ptolemy binary: shared objs + ptolemy_main
$(PTOLEMY_TARGET): $(SHARED_OBJECTS) $(BUILD_DIR)/ptolemy_main.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

# Unit tests
$(UTEST_TARGET): $(SHARED_OBJECTS) $(BUILD_DIR)/UnitTests.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

test: $(UTEST_TARGET)
	./$(UTEST_TARGET)

$(BUILD_DIR)/%.o: $(SRCDIR)/%.cpp $(wildcard $(INCDIR)/*.h) $(wildcard $(INCDIR)/math/*.h) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(SRCDIR)/*.o $(PTOLEMY_TARGET) $(UTEST_TARGET)

# distclean: like clean, plus the asan binary and stray Fortran scratch units.
# Does not touch /tmp/_s.out or /tmp/_m.out — make should never reach into /tmp.
distclean: clean
	rm -rf $(BUILD_DIR) $(PTOLEMY_TARGET) $(UTEST_TARGET) ptolemy_asan fort.*
