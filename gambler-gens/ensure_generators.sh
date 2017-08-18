#!/bin/bash

# Ensure that all generators are there (compiled)

mkdir -p bin

function compileSmallC # $1: generator name
{
	if [ ! -f bin/$1 ]; then
		gcc --std=c99 generators/$1.c -o bin/$1 -O3
	fi
}

function compileEstream # $1: generator name
{
	if [ ! -f bin/$1 ]; then
		gcc --std=c99 generators/common/ecrypt-sync.c generators/$1/$1.c generators/common/gen.c -o bin/$1 -iquote generators/common -iquote generators/$1 -O3
	fi
}

compileSmallC "spritz"
compileSmallC "vmpc"
compileSmallC "rc4p"
compileSmallC "c_rand"
compileSmallC "randu"
compileSmallC "hc128"
compileSmallC "mt19937ar"

if [ ! -f bin/los-rng ]; then
	g++ --std=c++11 generators/los-rng.cpp -o bin/los-rng -O3
fi

compileEstream "rabbit"
compileEstream "trivium"
compileEstream "sosemanuk"
compileEstream "salsa20"
compileEstream "grain"
compileEstream "mickey"
compileEstream "ffcsr"

