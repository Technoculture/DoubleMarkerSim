#!/bin/bash

# Check if Python 3.12 is installed
if ! command -v python3.12 &> /dev/null
then
    echo "Python 3.12 could not be found. Please install it and make
sure it's in your PATH."
    exit 1
fi

# Check if Poetry is installed
if ! command -v poetry &> /dev/null
then
    echo "Poetry could not be found."
    pip install poetry
    exit 1
fi

# Define the command function for Poetry run Python
py() {
    poetry run python "$@"
}

py src/down_syndrome_screening/analysis.py
