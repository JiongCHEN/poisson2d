#!/bin/bash

exe=build/main

echo "== Solve =="
time $exe $1

echo "== Visualize =="
gnuplot -e "set terminal png; splot '$1' with lines" --persist > func.png