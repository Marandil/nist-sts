#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hc128_ref.h"

int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: hc128 [length] [128bit key]\n");
		return -1;
	}
	char *pEnd;
	uint64 length = strtoull(argv[1], &pEnd, 10);
	char keystr[] = "00000000000000000000000000000000";
	size_t len = strlen(argv[2]);
	strcpy(&(keystr[32-len]), argv[2]);
	
	char *pos = keystr;
	uint8 key[16], iv[16];
	for(int i = 0; i < 16; ++i)
	{
		sscanf(pos, "%2hhx", &key[i]);
		iv[i] = 0;
		pos += 2;
	}
	
	HC128_State state;
	Initialization(&state, key, iv);
	for(int i=0; i<length; i+=4)
	{
		OneStep(&state);
		putchar(state.keystreamword & 0xff);
		putchar((state.keystreamword >> 8) & 0xff);
		putchar((state.keystreamword >> 16) & 0xff);
		putchar(state.keystreamword >> 24);
	}
	
	return 0;
}
