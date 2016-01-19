#!/bin/bash

gnuplot -e "set terminal png; splot '$1' with lines" --persist > func.png