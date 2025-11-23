#!/bin/bash

# Build script for Ancient Gene Predictor

set -e  # Exit on error

echo "========================================"
echo "  Ancient Gene Predictor Build Script"
echo "========================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
BUILD_TYPE="Release"
CLEAN=false
TESTS=ON
BENCHMARKS=ON
OPENMP=ON
INSTALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --no-tests)
            TESTS=OFF
            shift
            ;;
        --no-benchmarks)
            BENCHMARKS=OFF
            shift
            ;;
        --no-openmp)
            OPENMP=OFF
            shift
            ;;
        --install)
            INSTALL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --debug          Build in debug mode"
            echo "  --clean          Clean before building"
            echo "  --no-tests       Don't build tests"
            echo "  --no-benchmarks  Don't build benchmarks"
            echo "  --no-openmp      Disable OpenMP"
            echo "  --install        Install after building"
            echo "  -h, --help       Show this help"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Check for required tools
echo "Checking dependencies..."

if ! command -v cmake &> /dev/null; then
    echo -e "${RED}Error: cmake not found${NC}"
    echo "Please install cmake version 3.18 or higher"
    exit 1
fi

if ! command -v g++ &> /dev/null && ! command -v clang++ &> /dev/null; then
    echo -e "${RED}Error: No C++ compiler found${NC}"
    echo "Please install g++ or clang++"
    exit 1
fi

# Check zlib
if ! ldconfig -p | grep -q libz.so; then
    echo -e "${YELLOW}Warning: zlib not found${NC}"
    echo "You may need to install zlib-dev or zlib-devel"
fi

echo -e "${GREEN}✓ Dependencies OK${NC}"
echo ""

# Clean if requested
if [ "$CLEAN" = true ]; then
    echo "Cleaning build directory..."
    rm -rf build
    echo -e "${GREEN}✓ Cleaned${NC}"
    echo ""
fi

# Create build directory
mkdir -p build
cd build

# Configure
echo "Configuring with CMake..."
echo "  Build type: $BUILD_TYPE"
echo "  Tests: $TESTS"
echo "  Benchmarks: $BENCHMARKS"
echo "  OpenMP: $OPENMP"
echo ""

cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DBUILD_TESTS=$TESTS \
      -DBUILD_BENCHMARKS=$BENCHMARKS \
      -DUSE_OPENMP=$OPENMP \
      ..

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Configuration failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Configuration successful${NC}"
echo ""

# Build
echo "Building..."
make -j$(nproc)

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Build successful${NC}"
echo ""

# Run tests if built
if [ "$TESTS" = "ON" ] && [ -f "tests/agp_tests" ]; then
    echo "Running tests..."
    make test
    echo ""
fi

# Install if requested
if [ "$INSTALL" = true ]; then
    echo "Installing..."
    sudo make install
    echo -e "${GREEN}✓ Installation successful${NC}"
    echo ""
fi

# Summary
echo "========================================"
echo "  Build Complete!"
echo "========================================"
echo ""
echo "Executable: $(pwd)/agp"
echo ""
echo "To run:"
echo "  ./agp --help"
echo "  ./agp ../examples/test.fasta"
echo ""

if [ "$INSTALL" != true ]; then
    echo "To install system-wide:"
    echo "  sudo make install"
    echo ""
fi

echo "For more information, see README.md"
