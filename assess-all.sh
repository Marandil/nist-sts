#!/bin/bash

# Bit length of streams to be tested:
N=$((32*1024*1024*8))  # 32 MB

# List of generators to be tested
# [0] Input File (unsupported)
# [1] Linear Congruential
# [2] Quadratic Congruential I
# [3] Quadratic Congruential II
# [4] Cubic Congruential
# [5] XOR
# [6] Modular Exponentiation
# [7] Blum-Blum-Shub
# [8] Micali-Schnorr
# [9] G Using SHA-1
GENS="1 2 3 4 5 6 7 8 9"

# Set 0 if you only want to apply Gambler test
# Set 1 if you want to apply all of them
TESTS=1

# (modify if you want other tests if TESTS=0)
IF0="0000000000000001"

# Parameter Adjustments (0 to continue, number + value to set)
# e.g. "7 16384" will set number of runs per starting point to 16384
# if all tests are selected.
# The numbers below are tied to the tests selected and may vary,
# e.g. in TESTS=0 case Gambler has number 0.

#  [1] Block Frequency Test - block length(M):         128
#  [2] NonOverlapping Template Test - block length(m): 9
#  [3] Overlapping Template Test - block length(m):    9
#  [4] Approximate Entropy Test - block length(m):     10
#  [5] Serial Test - block length(m):                  16
#  [6] Linear Complexity Test - block length(M):       500
#  [7] Gambler - runs per starting point(M):           1000 

PARAMS="0"

# How many bitstreams to test?
BITSTREAMS=1024

TESTS_ARGS="$TESTS"
if (( $TESTS == 0 )); then
    TESTS_ARGS="$TESTS_ARGS $IF0"
fi

for GEN in $GENS; do
    CALL_PARAMS="$GEN $TESTS_ARGS $PARAMS $BITSTREAMS"
    time ./assess $N 2>&1 1> assess-$GEN.log <<< "$CALL_PARAMS" &
done

wait
