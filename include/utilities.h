/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
              U T I L I T Y  F U N C T I O N  P R O T O T Y P E S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdint.h>

int		displayGeneratorOptions();
int		generatorOptions(char** streamFile);
void	chooseTests();
void	fixParameters();
void	fileBasedBitStreams(char *streamFile);
void	readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile);
void	readHexDigitsInBinaryFormat(FILE *fp);
int		convertToBits(uint8_t *x, size_t xBitLength, size_t bitsNeeded, size_t *num_0s, size_t *num_1s, size_t *bitsRead);
void	openOutputStreams(int option);
void	invokeTestSuite(int option, char *streamFile);
void	nist_test_suite();

#define safe_scanf(pattern, value)  do { } while ( 1 != scanf((pattern), (value)) )
