#include <stdlib.h>
#include <stdio.h>

unsigned long long V;

void srandu(int seed)
{
	V = seed;
}

int randu()
{
	V = (65539ull * V) % 2147483648ull;
	return V;
}


int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: c_rand [length] [seed]\n");
		return -1;
	}
	int length = atoi(argv[1]);
	int seed = atoi(argv[2]);
	
	srandu(seed);
	
	for(int i=0; i<length; i+=2)
	{
		int r = randu();
		putchar(r & 0xff);
		putchar(r >> 8);
		putchar(r >> 16);
		putchar(r >> 24);
	}
	
	return 0;
}
