# Compiler settings
FC := ifx
FFLAGS := -fopenmp -O3 -xHost
LDFLAGS := -qmkl -fopenmp  # Removed -llapack -lblas as they're included in -qmkl

# Directories
SRC_DIR := source
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/

# Specific source files (in compilation order)
MOD_GLOBAL := $(SRC_DIR)/mod_global.f90
BOGOLIUBOV := $(SRC_DIR)/bogoliubov_transf.f90
MAIN_PROGRAM := $(SRC_DIR)/spinwavesNonCol.f90

# Object files
OBJ_GLOBAL := $(OBJ_DIR)/mod_global.o
OBJ_BOGOLIUBOV := $(OBJ_DIR)/bogoliubov_transf.o
OBJ_MAIN := $(OBJ_DIR)/spinwavesNonCol.o

# Module files (.mod) - these will be created in OBJ_DIR
MOD_GLOBAL_MOD := $(OBJ_DIR)/mod_global.mod
MOD_BOGOLIUBOV_MOD := $(OBJ_DIR)/bogoliubov_transf.mod

# All object files
OBJS := $(OBJ_GLOBAL) $(OBJ_BOGOLIUBOV) $(OBJ_MAIN)

# Target executable - named main.exe
TARGET := $(BIN_DIR)/main.exe

# Default target
all: $(BUILD_DIR) $(TARGET)

# Create build directories first
$(BUILD_DIR):
	@echo "Creating build directories..."
	mkdir -p $(OBJ_DIR)
	mkdir -p $(BIN_DIR)

# Link all object files to create executable
$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) -o $@
	@echo "Build complete! Executable created at: $@"
	@echo "Compiler: Intel Fortran (ifx) with OpenMP and MKL"

# Main program depends on both modules
$(OBJ_MAIN): $(MAIN_PROGRAM) $(MOD_GLOBAL_MOD) $(MOD_BOGOLIUBOV_MOD)
	$(FC) $(FFLAGS) -module $(OBJ_DIR) -c $< -o $@

# bogoliubov_transf depends on mod_global
$(OBJ_BOGOLIUBOV): $(BOGOLIUBOV) $(MOD_GLOBAL_MOD)
	$(FC) $(FFLAGS) -module $(OBJ_DIR) -c $< -o $@

# mod_global (independent module)
$(OBJ_GLOBAL): $(MOD_GLOBAL) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -module $(OBJ_DIR) -c $< -o $@

# Module files are created during compilation
# These rules ensure dependencies work correctly
$(MOD_GLOBAL_MOD): $(OBJ_GLOBAL)
	@touch $@  # Module file is created during compilation, just update timestamp

$(MOD_BOGOLIUBOV_MOD): $(OBJ_BOGOLIUBOV)
	@touch $@  # Module file is created during compilation, just update timestamp

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -rf $(BUILD_DIR)

# Clean and rebuild
rebuild: clean all

# Run the program
run: all
	@echo "Running main.exe..."
	@echo "==================="
	./$(TARGET)

# Show build info
info:
	@echo "Intel Fortran Build Information"
	@echo "==============================="
	@echo "Compiler: $(FC)"
	@echo "Compile flags: $(FFLAGS)"
	@echo "Link flags: $(LDFLAGS)"
	@echo "Module directory: $(OBJ_DIR)"
	@echo ""
	@echo "Source files:"
	@echo "  1. $(MOD_GLOBAL) (module)"
	@echo "  2. $(BOGOLIUBOV) (module, depends on mod_global)"
	@echo "  3. $(MAIN_PROGRAM) (main program)"
	@echo ""
	@echo "Target executable: $(TARGET)"
	@echo "Build directory: $(BUILD_DIR)"

# Test compilation (compile only)
compile: $(BUILD_DIR) $(OBJS)
	@echo "Compilation successful! Object files in $(OBJ_DIR)"

# Install (copy executable to current directory)
install: all
	@echo "Installing main.exe to current directory..."
	cp $(TARGET) ./main.exe
	@echo "Executable available as ./main.exe"

# Debug: Show what commands will be run
debug:
	@echo "FFLAGS: $(FFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "OBJ_DIR: $(OBJ_DIR)"
	@echo "Compilation commands:"
	@echo "  mod_global: $(FC) $(FFLAGS) -module $(OBJ_DIR) -c $(MOD_GLOBAL) -o $(OBJ_GLOBAL)"
	@echo "  bogoliubov: $(FC) $(FFLAGS) -module $(OBJ_DIR) -c $(BOGOLIUBOV) -o $(OBJ_BOGOLIUBOV)"
	@echo "  main: $(FC) $(FFLAGS) -module $(OBJ_DIR) -c $(MAIN_PROGRAM) -o $(OBJ_MAIN)"
	@echo "  link: $(FC) $(LDFLAGS) $(OBJS) -o $(TARGET)"

# Phony targets
.PHONY: all clean rebuild run info compile install debug help

# Help target
help:
	@echo "Intel Fortran Makefile for Spinwaves Program"
	@echo "============================================="
	@echo ""
	@echo "Project Structure:"
	@echo "  source/mod_global.f90          - Global definitions module"
	@echo "  source/bogoliubov_transf.f90   - Bogoliubov transformation module"
	@echo "  source/spinwavesNonCol.f90     - Main program"
	@echo ""
	@echo "Available targets:"
	@echo "  all      : Build the program (default)"
	@echo "  clean    : Remove all build files"
	@echo "  rebuild  : Clean and rebuild"
	@echo "  run      : Build and run the program"
	@echo "  compile  : Compile only (create object files)"
	@echo "  install  : Build and copy main.exe to current directory"
	@echo "  debug    : Show compilation commands for debugging"
	@echo "  info     : Show build information"
	@echo "  help     : Show this help message"
	@echo ""
	@echo "Compiler: Intel Fortran (ifx)"
	@echo "Features: OpenMP, Intel MKL"
	@echo "Optimization: -O3 -xHost"
	@echo ""
	@echo "Output: $(TARGET)"