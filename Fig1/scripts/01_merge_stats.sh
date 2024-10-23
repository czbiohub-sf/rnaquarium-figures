#!/bin/bash

# print header
head -n 1 "$1"
# print contents without header and sort/merge
tail -n+2 -q "$@" | (sort -us -t , -k 1,1)
