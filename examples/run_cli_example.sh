#!/bin/bash
# Basic usage of the metaclonotypist command line interface
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TCRPATH="$SCRIPT_DIR/data/tcrdata.csv"
HLAPATH="$SCRIPT_DIR/data/hladata.csv"
OUTPUTDIR="$SCRIPT_DIR/out/"

metaclonotypist --tcrpath "$TCRPATH" --hlapath "$HLAPATH" -o "$OUTPUTDIR" 
