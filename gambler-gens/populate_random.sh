#!/bin/bash


# make NSEQ random sequences
if [ -z "$NSEQ" ]; then
	NSEQ=4096
fi

# the sequences are BLEN bytes long
if [ -z "$BLEN" ]; then
	BLEN=$((256*1024))
fi

# list of supported generators
if [ -z "$GENS" ]; then
	# GENS="urandom openssl-rnd rc4 rc4-40 spritz vmpc rc4p aes128ctr aes192ctr aes256ctr hc128 rabbit trivium sosemanuk salsa20 grain mickey ffcsr mersenne mersenne_ar minstd oldbsd borland vs cmrg randu c_rand"
	GENS="rc4 rc4-40 spritz vmpc rc4p randu"
fi

# list of supported KDFs
if [ -z "$KDFS" ]; then
	KDFS="sha rel wep"
fi

for GEN in $GENS; do
	for KDF in $KDFS; do
		mkdir -p seq/$KDF/$GEN
	done
done

#
# Generate NSEQ random sequences from various sources:
# (the sequences have length of $BLEN, 256 KB by default)
#
# seq/N - the key/seed of each sequence is serial, i.e. 1, 2, 3, ... NSEQ
# seq/R - the key/seed of each sequence is psudorandom 
#
# seq/*/urandom/*		: sequence from /dev/urandom
# seq/*/openssl-rnd/*		: sequence from openssl rand function
# seq/*/rc4/*			: sequence from RC4 (128 bit key)
# seq/*/rc4-40/*		: sequence from RC4 (40 bit key)
# seq/*/spritz/*		: sequence from Spritz
# seq/*/vmpc/*			: sequence from VMPC (KSA variant)
# seq/*/rc4p/*			: sequence from RC4+
# seq/*/aes128ctr/*		: sequence from AES-128-CTR
# seq/*/aes192ctr/*		: sequence from AES-192-CTR
# seq/*/aes256ctr/*		: sequence from AES-256-CTR
# seq/*/hc128/*			: sequence from hc128
# seq/*/rabbit/*		: sequence from rabbit (eSTREAM implementation)
# seq/*/trivium/*		: sequence from trivium (eSTREAM implementation)
# seq/*/sosemanuk/*		: sequence from sosemanuk (eSTREAM implementation)
# seq/*/salsa20/*		: sequence from salsa20 (eSTREAM implementation)
# seq/*/grain/*			: sequence from grain (eSTREAM implementation)
# seq/*/mickey/*		: sequence from mickey (eSTREAM implementation)
# seq/*/ffcsr/*			: sequence from ffcsr (eSTREAM implementation 80bit key)
# seq/*/mersenne/*		: sequence from Mersenne Twister
# seq/*/mersenne_ar/*		: sequence from Mersenne Twister AR
# seq/*/c_rand/*		: sequence of first bytes of C rand() function
# seq/*/randu/*			: sequence of first two bytes of RANDU function

source ensure_generators.sh

function generate # $1: index number
{
	for GEN in $GENS; do
		for KDF in $KDFS; do
			FNAME="seq/$KDF/$GEN/$1"
			./generator.sh $KDF $GEN $BLEN $1 > $FNAME
		done
	done
}

function s_to_mm_ss # $1 seconds
{
	MIN=$(($1/60))
	SEC=$(($1%60))
	printf "%d:%02d" $MIN $SEC
}

for i in $(seq 1 $NSEQ); do
	generate $i
	if ! (($i % 8)) && (($SECONDS)); then
		ELAPSED=`s_to_mm_ss $SECONDS`
		CPS=$(($i/$SECONDS))
		if (($CPS)); then
			ETA=$((($NSEQ-$i)/$CPS))
			ETA=`s_to_mm_ss $ETA`
		else
			ETA="INF"
		fi
		printf "\r$i / $NSEQ   $CPS cps in $ELAPSED eta $ETA"
	fi
done
echo ""

echo "Merging single sequence files into collections..."

for GEN in $GENS; do
	for KDF in $KDFS; do
		cat seq/$KDF/$GEN/* > seq/$KDF/$GEN-sequence
	done
done
