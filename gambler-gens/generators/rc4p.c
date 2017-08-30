/********************
 * Implementation by Bartlomiej Surma based on:
 *
 * "Analysis of RC4 and Proposal of Additional Layers
 * for Better Security Margin" (2008)
 * by
 * Subhamoy Maitra and Goutam Paul
 *
 * Usage:
 * RC4p_State state;
 * initialization(&state, key, 16, iv, 16);
 * stream(&state, 256);	//256 bytes will be outputed to stdout
 *
 * U can use stream function consecutively as many times as you like.
 *
 * Note that (according to Subhamoy and Goutam)
 * iv should be of the same length as key
 * to avoid Time Memory Trade-Off attacks.
 ********************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void swap(unsigned char* arr, int a, int b)
{
	unsigned char tmp = arr[a];
	arr[a] = arr[b];
	arr[b] = tmp;
	/*
	if(a != b)
	{
		arr[a] ^= arr[b];
		arr[b] ^= arr[a];
		arr[a] ^= arr[b];
	}
	*/
}

typedef struct {
  unsigned char S[256];
  unsigned char i;
  unsigned char j;
} RC4p_State;

void layer1(RC4p_State* state, unsigned char* key, int key_len)
{
	for(int i = 0; i < 256; ++i) state->S[i] = i;
	state->j = 0;
	for(int i = 0; i < 256; ++i)
	{
		state->j = (state->j + state->S[i] + key[i % key_len]) & 255;
		swap(state->S, i, state->j);
	}
}

void layer2(RC4p_State* state, unsigned char* key, int key_len, unsigned char* iv, int iv_len)
{
	for(int i = 127; i >= 0; --i)
	{
		state->j = ((state->j + state->S[i]) ^ (key[i % key_len] + iv[i % iv_len])) & 255;
		swap(state->S, i, state->j);
	}

	for(int i = 128; i < 256; ++i)
	{
		state->j = ((state->j + state->S[i]) ^ (key[i % key_len] + iv[i % iv_len])) & 255;
		swap(state->S, i, state->j);
	}
}

void layer3(RC4p_State* state, unsigned char* key, int key_len)
{
	int i;
	for(int y=0; y<256; ++y)
	{
		i = y & 1 ? 256 - ((y + 1) >> 1) : y >> 1;
		state->j = (state->j + state->S[i] + key[i % key_len]) & 255;
		swap(state->S, i, state->j);
	}
}

void initialization(RC4p_State* state, unsigned char* key, int key_len, unsigned char* iv, int iv_len)
{
	layer1(state, key, key_len);
	layer2(state, key, key_len, iv, iv_len);
	layer3(state, key, key_len);
	state->i = 0;
	state->j = 0;
}

void stream(RC4p_State* state, size_t length)
{
	unsigned char t1, t2, t3;
	for(int k=0; k < length; ++k){
		state->i += 1;
		state->j += state->S[state->i];
		swap(state->S, state->i, state->j);
		t1 = (state->S[state->i] + state->S[state->j]) & 255;
		t2 = ((
			state->S[((state->i >> 3) ^ (state->j << 5)) & 255]
			+
			state->S[((state->i << 5) ^ (state->j >> 3)) & 255]
		) & 255) ^ 0xAA;
		t3 = (state->j + state->S[state->j]) & 255;
		putchar(((state->S[t1] + state->S[t2]) & 255) ^ state->S[t3]);
	}
}

int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: rc4a [length] [key]\n");
		return -1;
	}
	size_t length = strtoull(argv[1], NULL, 10);

  /* Key: hex to binary */
  size_t strLength = strlen(argv[2]);
  size_t keyLength = (strLength / 2) + (strLength % 2);

  char *keystr = malloc(keyLength * 2 + 1);
  unsigned char *key = malloc(keyLength);

  keystr[0] = '0';
	memcpy(keystr + (strLength % 2), argv[2], strLength);
  keystr[strLength] = 0; /* null-terminate string */

	char *pos = keystr;
  unsigned char *outpos = key;
	for(int i = 0; i < keyLength; ++i)
	{
		sscanf(pos, "%2hhx", outpos);
		pos += 2;
    outpos++;
	}

	unsigned char vec[16] = { 0 };

	RC4p_State state;
	initialization(&state, key, keyLength, vec, 16);
	stream(&state, length);

  free(keystr);
  free(key);

	return 0;
}
