#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

export PYTHONPATH=$(pwd)
echo "PYTHONPATH is set to $(pwd)"

echo "Installing required packages..."
pip install -r requirements.txt

if which g16 > /dev/null; then
    echo "Gaussian 16 found at $(which g16)"
else
    echo "Gaussian 16 not found in PATH. Please ensure Gaussian 16 is installed and available in PATH."
    exit 1
fi

echo "Setup completed successfully."