/*
 * A simple RC4A Implementation based on
 * "A New Weakness in the RC4 Keystream
 * Generator and an Approach to Improve the
 * Security of the Cipher"
 * by Souradyuti Paul and Bart Preneel
 */

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint8_t byte;

typedef byte RC4State[256];

typedef struct {
  RC4State S1;
  RC4State S2;
  int i, j1, j2;
} RC4A;

void swap(byte *S, int a, int b)
{
  uint8_t tmp = S[a];
  S[a] = S[b];
  S[b] = tmp;
}

void initialize(RC4A *rc4a, uint8_t *key, size_t keyLength)
{
  int i, j;

  /* Initialize RC4 state 1 */
  for(i=0; i<256; ++i)
  {
    rc4a->S1[i] = i;
  }
  j = 0;
  for(i=0; i<256; ++i)
  {
    j = (j + rc4a->S1[i] + key[i % keyLength]) % 256;
    swap(rc4a->S1, i, j);
  }

  /* Derive K2 from key, in this case we use K2 = S1 as "derived" */
  for(i=0; i<256; ++i)
  {
    rc4a->S2[i] = i;
  }
  j = 0;
  for(i=0; i<256; ++i)
  {
    j = (j + rc4a->S2[i] + rc4a->S1[i]) % 256;
    swap(rc4a->S2, i, j);
  }

  rc4a->i = 0;
  rc4a->j1 = 0;
  rc4a->j2 = 0;
}

void generate(RC4A *rc4a, size_t length)
{
  byte *S1 = rc4a->S1;
  byte *S2 = rc4a->S2;
  int i = rc4a->i, j1 = rc4a->j1, j2 = rc4a->j2;
  for(;length-->0;)
  {
    i++;
    j1 += S1[i];
    swap(S1, i, j1);
    putchar(S2[S1[i] + S1[j1]]);
    j2 += S2[i];
    swap(S2, i, j2);
    putchar(S1[S2[i] + S2[j2]]);
  }
  rc4a->i = i;
  rc4a->j1 = j1;
  rc4a->j2 = j2;
}

int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: rc4a [length] [key]\n");
		return -1;
	}
	int length = atoi(argv[1]);
	size_t len = strlen(argv[2]);
	size_t even_len = len + len % 2;
	char keystr[even_len];
	keystr[0] = '0';
	strcpy(&(keystr[len % 2]), argv[2]);

	char *pos = keystr;
	unsigned char key[even_len / 2];
	for(int i = 0; i < even_len; ++i)
	{
		sscanf(pos, "%2hhx", &key[i]);
		pos += 2;
	}

	RC4A state;
	initialize(&state, key, even_len / 2);
	generate(&state, length);
	return 0;
}
