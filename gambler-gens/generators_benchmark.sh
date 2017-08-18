FKEY="00000000000000000000000000000000"
IHEX="00000000000000000000000000000000" # just for IVs here

KEY32=`echo $FKEY | cut -c1-8`
DKEY32=`echo $((16#$KEY32))`
KEY80=`echo $FKEY | cut -c1-20`
KEY128=`echo $FKEY | cut -c1-32`
KEY192=`echo $FKEY | cut -c1-48`
KEY256=$FKEY

BLEN=65536
BLLEN=$(($BLEN/1024))

FNAME="tmp.out"

printf "\n/dev/urandom"
time dd if=/dev/urandom of=$FNAME count=$BLLEN bs=1KiB iflag=fullblock status=none
printf "\nopenssl rand"
time openssl rand -out $FNAME $BLEN
printf "\nopenssl rc4"
time head -c $BLEN /dev/zero | openssl rc4 -out $FNAME -K $KEY128
printf "\nspritz"
time bin/spritz $BLEN $KEY128 > $FNAME
printf "\nvmpc"
time bin/vmpc $BLEN $KEY128 > $FNAME
printf "\nrc4+"
time bin/rc4p $BLEN $KEY128 > $FNAME
printf "\naes128ctr"
time head -c $BLEN /dev/zero | openssl enc -aes-128-ctr -out $FNAME -K $KEY128 -iv $IHEX
printf "\naes192ctr"
time head -c $BLEN /dev/zero | openssl enc -aes-192-ctr -out $FNAME -K $KEY192 -iv $IHEX
printf "\naes256ctr"
time head -c $BLEN /dev/zero | openssl enc -aes-256-ctr -out $FNAME -K $KEY256 -iv $IHEX
printf "\nc_rand"
time bin/c_rand $BLEN $DKEY32 > $FNAME
printf "\nrandu"
time bin/randu $BLEN $DKEY32 > $FNAME
printf "\nhc128"
time bin/hc128 $BLEN $KEY128 > $FNAME
printf "\nrabbit"
time bin/rabbit $BLEN $KEY128 > $FNAME
printf "\ntrivium"
time bin/trivium $BLEN $KEY80 > $FNAME
printf "\nsosemanuk"
time bin/sosemanuk $BLEN $KEY256 > $FNAME
printf "\nsalsa20"
time bin/salsa20 $BLEN $KEY256 > $FNAME
printf "\ngrain"
time bin/grain $BLEN $KEY128 > $FNAME
printf "\nmickey"
time bin/mickey $BLEN $KEY128 > $FNAME
printf "\nf-fcsr"
time bin/ffcsr $BLEN $KEY80 > $FNAME

rm $FNAME
