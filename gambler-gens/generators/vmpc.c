

/*---------------------------------------------------------------------------------------------------

                 Implementation of the VMPC Stream Cipher
                                 and 
              the VMPC-MAC Authenticated Encryption Scheme
                                 in C

                Author of the algorithms: Bartosz Zoltak
              Author of the implementation: Bartosz Zoltak

                         www.vmpcfunction.com 

-----------------------------------------------------------------------------------------------------
----------------------- Usage of the algorithms: ----------------------------------------------------
-----------------------------------------------------------------------------------------------------

unsigned char Key[64], Vector[64]; Message[1000]; MessageMAC[20];

Encryption:

   VMPCInitKey(Key, Vector, 16, 16);
   VMPCEncrypt(Message, 1000);

Decryption:

   VMPCInitKey(Key, Vector, 16, 16);
   VMPCEncrypt(Message, 1000);

   (the VMPCEncrypt function is used for both encryption and decryption). 

Authenticated Encryption (with the MAC tag):

   VMPCInitKey(Key, Vector, 16, 16);
   VMPCInitMAC();
   VMPCEncryptMAC(Message, 1000);
   VMPCOutputMAC();                  //The MAC tag is saved in the 20-byte "MAC" table
   Move(MAC, MessageMAC, 20);        //Save the generated MAC tag in the "MessageMAC" table

Decryption and verification of the MAC tag:

   VMPCInitKey(Key, Vector, 16, 16);
   VMPCInitMAC();
   VMPCDecryptMAC(Message, 1000);
   VMPCOutputMAC();                  //The MAC tag is saved in the 20-byte "MAC" table

   If the 20-byte tables "MAC" and "MessageMAC" are identical, the message was correctly decrypted -
   the correct key was used and the message was not corrupted.

----------------------------------------------------------------------------------------------------
The VMPCInitKey / VMPCInitKey16 functions (employing the VMPC-KSA3 key initialization algorithm)
provide higher security level but about 1/3 lower efficiency.
than the basic VMPCInitKeyBASIC / VMPCInitKey16BASIC functions.

If only the system efficiency allows, the author recommends to use the VMPCInitKey / VMPCInitKey16 functions.
At the same time the VMPCInitKeyBASIC / VMPCInitKey16BASIC functions also remain secure. 
----------------------------------------------------------------------------------------------------
CAUTION! 
A DIFFERENT value of the initialization vector ("Vector")
should be used for each encryption with the same key ("Key").

Encrypting two messages with THE SAME key and THE SAME initialization vector
drastically reduces security level!

The key is a secret value.
The initialization vector is not secret - it can be passed in plain form along with the encrypted message.
-----------------------------------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//---------------------------------------------------------------------------------------------------
//----------------------------------------- IMPLEMENTATION: -----------------------------------------
//---------------------------------------------------------------------------------------------------

//--- The "int" type denotes a 32-bit integer


//----------- VMPC Stream Cipher variables: -----------
unsigned char P[256];
unsigned char s, n;

//----------- VMPC-MAC Authenticated Encryption Scheme variables: -----------
unsigned char MAC[32];
unsigned char m1, m2, m3, m4, mn;

//----------------- Test data: -----------------
unsigned char TestKey[16]         = {0x96, 0x61, 0x41, 0x0A, 0xB7, 0x97, 0xD8, 0xA9, 0xEB, 0x76, 0x7C, 0x21, 0x17, 0x2D, 0xF6, 0xC7};
unsigned char TestVector[16]      = {0x4B, 0x5C, 0x2F, 0x00, 0x3E, 0x67, 0xF3, 0x95, 0x57, 0xA8, 0xD2, 0x6F, 0x3D, 0xA2, 0xB1, 0x55};

unsigned char TestOutPInd[8]      = {0, 1, 2, 3, 252, 253, 254, 255};
unsigned int  TestOutInd[16]      = {0, 1, 2, 3, 252, 253, 254, 255, 1020, 1021, 1022, 1023, 102396, 102397, 102398, 102399};


unsigned char TestOutPBASIC[8]    = {0x3F, 0xA5, 0x22, 0x67, 0x75, 0xB3, 0xD2, 0xC3};
unsigned char TestOutBASIC[16]    = {0xA8, 0x24, 0x79, 0xF5, 0xB8, 0xFC, 0x66, 0xA4, 0xE0, 0x56, 0x40, 0xA5, 0x81, 0xCA, 0x49, 0x9A};
      //VMPCInitKeyBASIC(TestKey, TestVector, 16, 16);  OR  VMPCInitKey16BASIC(TestKey, TestVector);
      //P[TestOutPInd[x]]==TestOutPBASIC[x];  x=0,1,...,7
      //Table[x]=0;  x=0,1,...,102399
      //VMPCEncrypt(Table, 102400);  OR  VMPCEncryptMAC(Table, 102400);
      //Table[TestOutInd[x]]==TestOutBASIC[x];  x=0,1,...,15


unsigned char TestOutP[8]         = {0x1F, 0x00, 0xE2, 0x03, 0x5C, 0xEE, 0xC2, 0x2B};
unsigned char TestOut[16]         = {0xB6, 0xEB, 0xAE, 0xFE, 0x48, 0x17, 0x24, 0x73, 0x1D, 0xAE, 0xC3, 0x5A, 0x1D, 0xA7, 0xE1, 0xDC};
      //VMPCInitKey(TestKey, TestVector, 16, 16);  OR  VMPCInitKey16(TestKey, TestVector);
      //P[TestOutPInd[x]]==TestOutP[x];  x=0,1,...,7
      //Table[x]=0;  x=0,1,...,102399
      //VMPCEncrypt(Table, 102400);  OR  VMPCEncryptMAC(Table, 102400);
      //Table[TestOutInd[x]]==TestOut[x];  x=0,1,...,15


unsigned char TestOutMACBASIC[20] = {0x9B, 0xDA, 0x16, 0xE2, 0xAD, 0x0E, 0x28, 0x47, 0x74, 0xA3, 0xAC, 0xBC, 0x88, 0x35, 0xA8, 0x32, 0x6C, 0x11, 0xFA, 0xAD};
      //Table[x]=x;  x=0,1,2,...,254,255
      //VMPCInitKeyBASIC(TestKey, TestVector, 16, 16);  OR  VMPCInitKey16BASIC(TestKey, TestVector);
      //VMPCInitMAC();
      //VMPCEncryptMAC(Table, 256);
      //VMPCOutputMAC();
      //MAC[x]==TestOutMACBASIC[x];  x=0,1,...,19

unsigned char TestOutMAC[20]      = {0xA2, 0xB6, 0x0D, 0xB7, 0xB3, 0x90, 0x1D, 0x5C, 0x99, 0x61, 0x7C, 0xE2, 0xA3, 0x95, 0x02, 0x81, 0x75, 0x3A, 0x0C, 0x98};
      //Table[x]=x & 255;  x=0,1,2,...,999998,999999;  (Table[0]=0; Table[1]=1; ...; Table[999998]=62; Table[999999]=63)
      //VMPCInitKey(TestKey, TestVector, 16, 16);  OR  VMPCInitKey16(TestKey, TestVector);
      //VMPCInitMAC();
      //VMPCEncryptMAC(Table, 1000000);
      //VMPCOutputMAC();
      //MAC[x]==TestOutMAC[x];  x=0,1,...,19

//-----------------------------------------------------------------------------------------------------------

unsigned char Permut123[256]=   //Permut123[x]=x
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,
72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,
109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,
139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,
169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,
199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,
228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255};

//-----------------------------------------------------------------------------------------------------------


//---------- VMPC Stream Cipher: ----------


void VMPCInitKeyRound(unsigned char Data[], unsigned char Len, unsigned char Src)
/*Data: key or initialization vector
  Len=1,2,3,...,64: key/initialization vector length (in bytes)
  Src=0: first initialization of the key (the P table and the s variable will be restored to their initial values first)
  Src=1: re-initialization of the key, e.g. with the initialization vector*/
{
  if (Src==0) {
    memcpy(P, Permut123, 256);
    s=0;
  }
  unsigned char k=0;
  n=0;

  for (int x=0; x<768; x++)
  {
    s=P[ (s + P[n] + Data[k]) & 255 ];

    unsigned char t=P[n];  P[n]=P[s];  P[s]=t;

    k++;  if (k==Len) k=0;
    n++;
  }
}


void VMPCInitKey(unsigned char Key[], unsigned char Vec[], unsigned char KeyLen, unsigned char VecLen)   //KeyLen, VecLen = 1,2,3,...,64
{
  VMPCInitKeyRound(Key, KeyLen, 0);
  VMPCInitKeyRound(Vec, VecLen, 1);
  VMPCInitKeyRound(Key, KeyLen, 1);
}


void VMPCInitKeyBASIC(unsigned char Key[], unsigned char Vec[], unsigned char KeyLen, unsigned char VecLen)   //KeyLen, VecLen = 1,2,3,...,64
{
  VMPCInitKeyRound(Key, KeyLen, 0);
  VMPCInitKeyRound(Vec, VecLen, 1);
}



void VMPCEncrypt(unsigned int Len)
{
  for (unsigned int x=0; x<Len; x++)
  {
    s=P[ (s + P[n]) & 255 ];

    putchar(P[(P[P[ s ]]+1) & 255]);

    unsigned char t=P[n];  P[n]=P[s];  P[s]=t;

    n++;
  }
}

//---------------------------------------------------------------------------------------------------
//----------------------------------------------- END -----------------------------------------------
//---------------------------------------------------------------------------------------------------


int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: vmpc [length] [key]\n");
		return -1;
	}
	int length = atoi(argv[1]);
	size_t len = strlen(argv[2]);
	size_t even_len = len + len % 2;
	char keystr[even_len];
	keystr[0] = '0';
	strcpy(&(keystr[len % 2]), argv[2]);
	
	char *pos = keystr;
	unsigned char key[even_len / 2], vec[16];
	for(int i = 0; i < 16; ++i) vec[i] = 0;
	for(int i = 0; i < even_len; ++i)
	{
		sscanf(pos, "%2hhx", &key[i]);
		pos += 2;
	}
	
	//* //for KSA version
	VMPCInitKeyBASIC(key, vec, even_len / 2, 16);
	/*/  //for KSA3 version
	VMPCInitKey(key, vec, even_len / 2, 16);
	//*/
	VMPCEncrypt(length);
	
	return 0;
}