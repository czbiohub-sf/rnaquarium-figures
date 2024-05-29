#!/bin/bash

FILES=(data/stats-*.csv)
# print header
head -n 1 ${FILES[0]}
# print contents without header and sort/merge
tail -n+2 -q ${FILES[@]} | (sort -us -t , -k 1,1)
