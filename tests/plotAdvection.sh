#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

rm "$DIR"/advection*.dat
"$DIR"/runAdvection
gnuplot "$DIR"/displayAdvection.gnu