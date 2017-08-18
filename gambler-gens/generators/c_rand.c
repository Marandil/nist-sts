#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "Usage: c_rand [length] [seed]\n");
		return -1;
	}
	int length = atoi(argv[1]);
	int seed = atoi(argv[2]);
	
	srand(seed);
	
	for(int i=0; i<length; ++i)
	{
		int r = rand();
		putchar(r & 0xff);
	}
	
	return 0;
}
