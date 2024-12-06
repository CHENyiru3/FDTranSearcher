#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <mode> <args...>"
    echo "mode: r (reference-based) or d (de novo)"
    exit 1
fi

MODE=$1
shift  


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


if [ "$MODE" = "r" ]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/reference_based_module/reference_based_main.py"
elif [ "$MODE" = "d" ]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/de_novo_module/transposon_analyzer.py"
else
    echo "Invalid mode. Use 'r' for reference-based or 'd' for de novo."
    exit 1
fi

python3 "$PYTHON_SCRIPT" "$@"