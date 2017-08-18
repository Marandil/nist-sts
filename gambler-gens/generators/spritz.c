/**********************
 * Code taken from https://github.com/jedisct1/spritz
 * Written in 2014-2016 by Frank Denis.
 * see the licence here: http://creativecommons.org/publicdomain/zero/1.0/
 * 
 * Modified by Bartlomiej Surma
 * ********************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
//#include "spritz.h"

#define N 256

#if defined(_MSC_VER)
# define ALIGNED(S) __declspec(align(S))
#elif defined(__GNUC__)
# define ALIGNED(S) __attribute__((aligned(S)))
#else
# define ALIGNED(S)
#endif

ALIGNED(64) typedef struct State_ {
    unsigned char s[N];
    unsigned char a;
    unsigned char i;
    unsigned char j;
    unsigned char k;
    unsigned char w;
    unsigned char z;
} State;

#define LOW(B)  ((B) & 0xf)
#define HIGH(B) ((B) >> 4)

static void
initialize_state(State *state)
{
    unsigned int v;

    for (v = 0; v < N; v++) {
        state->s[v] = (unsigned char) v;
    }
    state->a = 0;
    state->i = 0;
    state->j = 0;
    state->k = 0;
    state->w = 1;
    state->z = 0;
}

static void
update(State *state)
{
    unsigned char t;
    unsigned char y;

    state->i += state->w;
    y = state->j + state->s[state->i];
    state->j = state->k + state->s[y];
    state->k = state->i + state->k + state->s[state->j];
    t = state->s[state->i];
    state->s[state->i] = state->s[state->j];
    state->s[state->j] = t;
}

static unsigned char
output(State *state)
{
    const unsigned char y1 = state->z + state->k;
    const unsigned char x1 = state->i + state->s[y1];
    const unsigned char y2 = state->j + state->s[x1];

    state->z = state->s[y2];

    return state->z;
}

static void
crush(State *state)
{
    unsigned char v;
    unsigned char x1;
    unsigned char x2;
    unsigned char y;

    for (v = 0; v < N / 2; v++) {
        y = (N - 1) - v;
        x1 = state->s[v];
        x2 = state->s[y];
        if (x1 > x2) {
            state->s[v] = x2;
            state->s[y] = x1;
        } else {
            state->s[v] = x1;
            state->s[y] = x2;
        }
    }
}

static void
whip(State *state)
{
    const unsigned int r = N * 2;
    unsigned int       v;

    for (v = 0; v < r; v++) {
        update(state);
    }
    state->w += 2;
}

static void
shuffle(State *state)
{
    whip(state);
    crush(state);
    whip(state);
    crush(state);
    whip(state);
    state->a = 0;
}

static void
absorb_nibble(State *state, const unsigned char x)
{
    unsigned char t;
    unsigned char y;

    if (state->a == N / 2) {
        shuffle(state);
    }
    y = N / 2 + x;
    t = state->s[state->a];
    state->s[state->a] = state->s[y];
    state->s[y] = t;
    state->a++;
}

static void
absorb_byte(State *state, const unsigned char b)
{
    absorb_nibble(state, LOW(b));
    absorb_nibble(state, HIGH(b));
}

static void
absorb(State *state, const unsigned char *msg, size_t length)
{
    size_t v;

    for (v = 0; v < length; v++) {
        absorb_byte(state, msg[v]);
    }
}

static unsigned char
drip(State *state)
{
    if (state->a > 0) {
        shuffle(state);
    }
    update(state);

    return output(state);
}

static void
squeeze(State *state, size_t outlen)
{
    size_t v;

    if (state->a > 0) {
        shuffle(state);
    }
    for (v = 0; v < outlen; v++) {
        putchar(drip(state));
    }
}

void
spritz_stream(size_t outlen, const unsigned char *key, size_t keylen)
{
    State state;

    initialize_state(&state);
    absorb(&state, key, keylen);
    squeeze(&state, outlen);
}

int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: spritz [length] [key]\n");
		return -1;
	}
	int length = atoi(argv[1]);
	size_t len = strlen(argv[2]);
	size_t even_len = len + len % 2;
	char keystr[even_len];
	keystr[0] = '0';
	strcpy(&(keystr[len % 2]), argv[2]);
	
	char *pos = keystr;
	const unsigned char key[even_len / 2];
	for(int i = 0; i < even_len; ++i)
	{
		sscanf(pos, "%2hhx", &key[i]);
		pos += 2;
	}
	
	spritz_stream(length, key, even_len / 2);
	
	return 0;
}