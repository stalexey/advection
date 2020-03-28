#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

rm "$DIR"/interpolation*.dat
"$DIR"/runInterpolation
gnuplot "$DIR"/displayInterpolation.gnu