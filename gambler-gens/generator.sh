#!/bin/bash

#>&2 echo "start - $1 $2 $3 $4"

#
# Usage: generator.sh [kdf: sha|hex|echo] [generator] [length] [keybase]

DATE=`date`
echo "$DATE $1 $2 $3 $4" >> generator.log

# set the key derivation function [possible are sha, echo, rel, wep and hex]
KDF=$1

# set the default IV for the ciphers
IHEX=00000000000000000000000000000000

source ensure_generators.sh

# base to be added to all sha kdfs
if [ -z $ADD_BASE ]; then
	ADD_BASE="GAMBLER_0001"
fi

function sha # $1 keybase # hashes keybase and returns hex
{
	printf "$ADD_BASE#$1" | sha256sum | cut -c1-64
}

function hex # $1 keybase # interprets key as integer and returns hex repr.
{
	printf "%064x\n" $1
}

function rel # $1 keybase 
{
	PREFIX=`sha $ADD_BASE | cut -c1-48`
	SUFFIX=`hex $1 | cut -c49-64`
	echo "$PREFIX$SUFFIX"
}

function wep # $1 keybase 
{
	PREFIX=`hex $1 | cut -c27-64`          #128 + 24 bits (first 128 gets cut with kdf128)
	SUFFIX=`sha $ADD_BASE | cut -c1-26`    #104 bits
	echo "$PREFIX$SUFFIX"
}

function kdf32 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c57-64
}

function kdf32d # $1 keybase
{
	KEY=`kdf32 $1`
	echo $((16#$KEY))
}

function kdf40 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c55-64
}

function kdf64 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c49-64
}

function kdf64d # $1 keybase
{
	KEY=`kdf64 $1`
	echo $((16#$KEY))
}

function kdf80 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c45-64
}

function kdf128 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c33-64
}

function kdf192 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c17-64
}

function kdf256 # $1 keybase
{
	KEY=`$KDF $1`
	echo $KEY | cut -c1-64
}


# generator functions:
# [generator] [length in bytes] [keybase in hex]

function urandom # [length] [keybase] # keybase ignored for urandom
{
	BLLEN=$(($1/1024))
	dd if=/dev/urandom of=/dev/stdout count=$BLLEN bs=1KiB iflag=fullblock status=none
}

function openssl-rnd # [length] [keybase] # keybase ignored for openssl rand
{
	openssl rand $1
}

function rc4 # [length] [keybase]
{
	KEY128=`kdf128 $2`
	head -c $1 /dev/zero | openssl rc4 -K $KEY128
}
	
function rc4-40 # [length] [keybase]
{
	KEY40=`kdf40 $2`
	head -c $1 /dev/zero | openssl rc4 -K $KEY40
}
	
function spritz # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/spritz $1 $KEY128
}

function vmpc # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/vmpc $1 $KEY128
}
	
function rc4p # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/rc4p $1 $KEY128
}

function aes128ctr # [length] [keybase]
{
	KEY128=`kdf128 $2`
	head -c $1 /dev/zero | openssl enc -aes-128-ctr -iv $IHEX -K $KEY128 2>tmp.log
}

function aes192ctr # [length] [keybase]
{
	KEY192=`kdf192 $2`
	head -c $1 /dev/zero | openssl enc -aes-192-ctr -iv $IHEX -K $KEY192 2>tmp.log
}

function aes256ctr # [length] [keybase]
{
	KEY256=`kdf256 $2`
	head -c $1 /dev/zero | openssl enc -aes-256-ctr -iv $IHEX -K $KEY256 2>tmp.log
}

function hc128 # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/hc128 $1 $KEY128
}

function rabbit # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/rabbit $1 $KEY128
}

function trivium # [length] [keybase]
{
	KEY80=`kdf80 $2`
	bin/trivium $1 $KEY80
}

function sosemanuk # [length] [keybase]
{
	KEY256=`kdf256 $2`
	bin/sosemanuk $1 $KEY256
}
	
function salsa20 # [length] [keybase]
{
	KEY256=`kdf256 $2`
	bin/salsa20 $1 $KEY256
}
	
function grain # [length] [keybase]
{
	KEY128=`kdf128s $2`
	bin/grain $1 $KEY128
}
	
function mickey # [length] [keybase]
{
	KEY128=`kdf128 $2`
	bin/mickey $1 $KEY128
}
	
function ffcsr # [length] [keybase]
{
	KEY80=`kdf80 $2`
	bin/ffcsr $1 $KEY80
}

function mersenne # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng Mersenne $1 $DKEY64
}

function mersenne_ar # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/mt19937ar $1 $DKEY64
}

function minstd # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng Minstd $1 $DKEY64
}

function oldbsd # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng Rand $1 $DKEY64
}

function borland # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng Borland $1 $DKEY64
}

function vs # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng VS $1 $DKEY64
}

function cmrg # [length] [keybase]
{
	DKEY64=`kdf64d $2`
	bin/los-rng CMRG $1 $DKEY64
}

function randu # [length] [keybase]
{
	DKEY32=`kdf32d $2`
	bin/los-rng RANDU $1 $DKEY32
	#bin/randu $1 $DKEY32
}

function c_rand # [length] [keybase]
{
	DKEY32=`kdf32d $2`
	bin/los-rng C_PRG $1 $DKEY32
	#bin/c_rand $1 $DKEY32
}

$2 $3 $4

# close the file descriptor:
exec 1>&-

#>&2 echo "stop  - $1 $2 $3 $4"

