#!/bin/bash

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Directory for the release
RELEASE_DIR="../../releases"
VERSION_DIR="$RELEASE_DIR/singer_pair-$VERSION-alpha-linux-x86_64"

# Create version directory
mkdir -p $VERSION_DIR

# Compile the program with optimizations and debugging information
g++ -std=c++17 -O3 -g -static *.cpp -o $VERSION_DIR/singer_pair

# Compile the debug version of the program
g++ -std=c++17 -g -static *.cpp -o $VERSION_DIR/singer_pair_debug


