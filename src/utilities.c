/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
U T I L I T I E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/generators.h"
#include "../include/stat_fncs.h"

#include "../include/memutils.h"

int
displayGeneratorOptions()
{
	int		option = 0;

	printf("           G E N E R A T O R    S E L E C T I O N \n");
	printf("           ______________________________________\n\n");
	printf("    [0] Input File                 [1] Linear Congruential\n");
	printf("    [2] Quadratic Congruential I   [3] Quadratic Congruential II\n");
	printf("    [4] Cubic Congruential         [5] XOR\n");
	printf("    [6] Modular Exponentiation     [7] Blum-Blum-Shub\n");
	printf("    [8] Micali-Schnorr             [9] G Using SHA-1\n\n");
	printf("   Enter Choice: ");
	safe_scanf("%d", &option);
	printf("\n\n");

	return option;
}


int
generatorOptions(char** streamFile)
{
	char	file[200];
	int		option = NUMOFGENERATORS+1;
	FILE	*fp;

	while ( (option < 0) || (option > NUMOFGENERATORS) ) {
		option = displayGeneratorOptions();
		switch( option ) {
			case 0:
				printf("\t\tUser Prescribed Input File: ");
				safe_scanf("%s", file);
				*streamFile = (char*)calloc(200, sizeof(char));
				sprintf(*streamFile, "%s", file);
				printf("\n");
				if ( (fp = fopen(*streamFile, "r")) == NULL ) {
					printf("File Error:  file %s could not be opened.\n",  *streamFile);
					exit(-1);
				}
				else
					fclose(fp);
				break;
			case 1:
				*streamFile = "Linear-Congruential";
				break;
			case 2:
				*streamFile = "Quadratic-Congruential-1";
				break;
			case 3:
				*streamFile = "Quadratic-Congruential-2";
				break;
			case 4:
				*streamFile = "Cubic-Congruential";
				break;
			case 5:
				*streamFile = "XOR";
				break;
			case 6:
				*streamFile = "Modular-Exponentiation";
				break;
			case 7:
				*streamFile = "Blum-Blum-Shub";
				break;
			case 8:
				*streamFile = "Micali-Schnorr";
				break;
			case 9:
				*streamFile = "G using SHA-1";
				break;

			/* INTRODUCE NEW PRNG NAMES HERE */
			/*
			case 10:  *streamFile = "myNewPRNG";
				break;
			*/
			default:
				printf("Error:  Out of range - Try again!\n");
				break;
		}
	}
	return option;
}


void
chooseTests()
{
	int		i;

	printf("                S T A T I S T I C A L   T E S T S\n");
	printf("                _________________________________\n\n");
	printf("    [01] Frequency                       [02] Block Frequency\n");
	printf("    [03] Cumulative Sums                 [04] Runs\n");
	printf("    [05] Longest Run of Ones             [06] Rank\n");
	printf("    [07] Discrete Fourier Transform      [08] Nonperiodic Template Matchings\n");
	printf("    [09] Overlapping Template Matchings  [10] Universal Statistical\n");
	printf("    [11] Approximate Entropy             [12] Random Excursions\n");
	printf("    [13] Random Excursions Variant       [14] Serial\n");
	printf("    [15] Linear Complexity               [16] Gambler\n\n");
	printf("         INSTRUCTIONS\n");
	printf("            Enter 0 if you DO NOT want to apply all of the\n");
	printf("            statistical tests to each sequence and 1 if you DO.\n\n");
	printf("   Enter Choice: ");
	safe_scanf("%d", &testVector[0]);
	printf("\n");
	if ( testVector[0] == 1 )
		for( i=1; i<=NUMOFTESTS; i++ )
			testVector[i] = 1;
	else {
		printf("         INSTRUCTIONS\n");
		printf("            Enter a 0 or 1 to indicate whether or not the numbered statistical\n");
		printf("            test should be applied to each sequence.\n\n");
		printf("      1234567891111111\n");
		printf("               0123456\n");
		printf("      ");
		for ( i=1; i<=NUMOFTESTS; i++ )
			safe_scanf("%1d", &testVector[i]);
		printf("\n\n");
	}
}


void
fixParameters()
{
	int		counter, testid;

	//  Check to see if any parameterized tests are selected
	if ( (testVector[TEST_BLOCK_FREQUENCY] != 1) && (testVector[TEST_NONPERIODIC] != 1) &&
		 (testVector[TEST_OVERLAPPING] != 1) && (testVector[TEST_APEN] != 1) &&
		 (testVector[TEST_SERIAL] != 1) && (testVector[TEST_LINEARCOMPLEXITY] != 1)  &&
		 (testVector[TEST_GAMBLER] != 1) )
			return;

	do {
		counter = 1;
		printf("        P a r a m e t e r   A d j u s t m e n t s\n");
		printf("        -----------------------------------------\n");
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 )
			printf("    [%d] Block Frequency Test - block length(M):         %zd\n", counter++, tp.blockFrequencyBlockLength);
		if ( testVector[TEST_NONPERIODIC] == 1 )
			printf("    [%d] NonOverlapping Template Test - block length(m): %zd\n", counter++, tp.nonOverlappingTemplateBlockLength);
		if ( testVector[TEST_OVERLAPPING] == 1 )
			printf("    [%d] Overlapping Template Test - block length(m):    %zd\n", counter++, tp.overlappingTemplateBlockLength);
		if ( testVector[TEST_APEN] == 1 )
			printf("    [%d] Approximate Entropy Test - block length(m):     %zd\n", counter++, tp.approximateEntropyBlockLength);
		if ( testVector[TEST_SERIAL] == 1 )
			printf("    [%d] Serial Test - block length(m):                  %zd\n", counter++, tp.serialBlockLength);
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 )
			printf("    [%d] Linear Complexity Test - block length(M):       %zd\n", counter++, tp.linearComplexitySequenceLength);
		if ( testVector[TEST_GAMBLER] == 1 )
		{
			printf("    [%d] Gambler -\truns per starting point(M):          %zd\n", counter++, tp.gamblerRunsPerStartingPoint);
			printf("                  \tfirst starting point (s_start):      %zd\n", tp.gamblerStartStartingPoint);
			printf("                  \tlast starting point (s_end):         %zd\n", tp.gamblerEndStartingPoint);
			printf("                  \ttest mask (mask):                    %zd\n", tp.gamblerTestMask);
		}
		printf("\n");
		printf("   Select Test (0 to continue): ");
		safe_scanf("%1d", &testid);
		printf("\n");

		counter = 0;
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Block Frequency Test block length: ");
				safe_scanf("%zd", &tp.blockFrequencyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_NONPERIODIC] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter NonOverlapping Template Test block Length: ");
				safe_scanf("%zd", &tp.nonOverlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_OVERLAPPING] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Overlapping Template Test block Length: ");
				safe_scanf("%zd", &tp.overlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_APEN] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Approximate Entropy Test block Length: ");
				safe_scanf("%zd", &tp.approximateEntropyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_SERIAL] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Serial Test block Length: ");
				safe_scanf("%zd", &tp.serialBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Linear Complexity Test block Length: ");
				safe_scanf("%zd", &tp.linearComplexitySequenceLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_GAMBLER] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("         Last starting point == 0 equals the last possible\n");
				printf("         Test mask is a bit mask: T1 - 1\n");
				printf("                                  T2 - 2\n");
				printf("                                  T3 - 4\n");
				printf("                                  T4 - 8\n");

				printf("   Enter number of runs per Gambler's starting point: ");
				safe_scanf("%zd", &tp.gamblerRunsPerStartingPoint);
				printf("\n");
				printf("   Enter the first Gambler's starting point: ");
				safe_scanf("%zd", &tp.gamblerStartStartingPoint);
				printf("\n");
				printf("   Enter the last Gambler's starting point: ");
				safe_scanf("%zd", &tp.gamblerEndStartingPoint);
				printf("\n");
				printf("   Enter the Gambler's test bit mask: ");
				safe_scanf("%zd", &tp.gamblerTestMask);
				printf("\n");
				continue;
			}
		}
	} while ( testid != 0 );
}


void
fileBasedBitStreams(char *streamFile)
{
	FILE	*fp;
	int		mode;

	printf("   Input File Format:\n");
	printf("    [0] ASCII - A sequence of ASCII 0's and 1's\n");
	printf("    [1] Binary - Each byte in data file contains 8 bits of data\n\n");
	printf("   Select input mode:  ");
	safe_scanf("%1d", &mode);
	printf("\n");
	if ( mode == 0 ) {
		if ( (fp = fopen(streamFile, "r")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n",  streamFile);
			exit(-1);
		}
		readBinaryDigitsInASCIIFormat(fp, streamFile);
		fclose(fp);
	}
	else if ( mode == 1 ) {
		if ( (fp = fopen(streamFile, "rb")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n", streamFile);
			exit(-1);
		}
		readHexDigitsInBinaryFormat(fp);
		fclose(fp);
	}
}


void
readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile)
{
	size_t   i, j, num_0s, num_1s, bitsRead;
	unsigned bit;

	if ( (epsilon = (BitSequence *) mcalloc(tp.n, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	printf("     Statistical Testing In Progress.........\n\n");
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		for ( j=0; j<tp.n; j++ ) {
			if ( fscanf(fp, "%1d", &bit) == EOF ) {
				printf("ERROR:  Insufficient data in file %s.  %zd bits were read.\n", streamFile, bitsRead);
				fclose(fp);
				mfree(epsilon);
				return;
			}
			else {
				bitsRead++;
				if ( bit == 0 )
					num_0s++;
				else
					num_1s++;
				epsilon[j] = bit;
			}
		}
		fprintf(freqfp, "\t\tBITSREAD = %zd 0s = %zd 1s = %zd\n", bitsRead, num_0s, num_1s);
		nist_test_suite();
	}
	free(epsilon);
}


void
readHexDigitsInBinaryFormat(FILE *fp)
{
	size_t  i, done, num_0s, num_1s, bitsRead;
	uint8_t buffer[4];

	if ( (epsilon = (BitSequence *) mcalloc(tp.n,sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		return;
	}

	printf("     Statistical Testing In Progress.........\n\n");
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		done = 0;
		do {
			if ( fread(buffer, sizeof(unsigned char), 4, fp) != 4 ) {
				printf("READ ERROR:  Insufficient data in file.\n");
				free(epsilon);
				return;
			}
			done = convertToBits(buffer, 32, tp.n, &num_0s, &num_1s, &bitsRead);
		} while ( !done );
		fprintf(freqfp, "\t\tBITSREAD = %zd 0s = %zd 1s = %zd\n", bitsRead, num_0s, num_1s);

		nist_test_suite();

	}
	mfree(epsilon);
}


int
convertToBits(BYTE *x, size_t xBitLength, size_t bitsNeeded, size_t *num_0s, size_t *num_1s, size_t *bitsRead)
{
	size_t  i, j, count, bit;
	BYTE	  mask;
	size_t  zeros, ones;

	count = 0;
	zeros = ones = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
				(*num_1s)++;
				ones++;
			}
			else {
				bit = 0;
				(*num_0s)++;
				zeros++;
			}
			mask >>= 1;
			epsilon[*bitsRead] = bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}

	return 0;
}


void
openOutputStreams(int option)
{
	int		i, numOfBitStreams, numOfOpenFiles = 0;
	char	freqfn[200], summaryfn[200], statsDir[200], resultsDir[200];

	sprintf(freqfn, "experiments/%s/freq.txt", generatorDir[option]);
	if ( (freqfp = fopen(freqfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open freq file: <%s>", freqfn);
		exit(-1);
	}
	sprintf(summaryfn, "experiments/%s/finalAnalysisReport.txt", generatorDir[option]);
	if ( (summary = fopen(summaryfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open stats file: <%s>", summaryfn);
		exit(-1);
	}

	for( i=1; i<=NUMOFTESTS; i++ ) {
		if ( testVector[i] == 1 ) {
			sprintf(statsDir, "experiments/%s/%s/stats.txt", generatorDir[option], testNames[i]);
			sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
			if ( (stats[i] = fopen(statsDir, "w")) == NULL ) {	/* STATISTICS LOG */
				printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				exit(-1);
			}
			else
				numOfOpenFiles++;
			if ( (results[i] = fopen(resultsDir, "w")) == NULL ) {	/* P_VALUES LOG   */
				 printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				 printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				 printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				 exit(-1);
			}
			else
				numOfOpenFiles++;
		}
	}
	printf("   How many bitstreams? ");
	safe_scanf("%d", &numOfBitStreams);
	tp.numOfBitStreams = numOfBitStreams;
	printf("\n");
}


void
invokeTestSuite(int option, char *streamFile)
{
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	fprintf(freqfp, "\t\tFILE = %s\t\tALPHA = %6.4f\n", streamFile, ALPHA);
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	if ( option != 0 )
		printf("     Statistical Testing In Progress.........\n\n");
	switch( option ) {
		case 0:
			fileBasedBitStreams(streamFile);
			break;
		case 1:
			lcg();
			break;
		case 2:
			quadRes1();
			break;
		case 3:
			quadRes2();
			break;
		case 4:
			cubicRes();
			break;
		case 5:
			exclusiveOR();
			break;
		case 6:
			modExp();
			break;
		case 7:
			bbs();
			break;
		case 8:
			micali_schnorr();
			break;
		case 9:
			SHA1();
			break;

		/* INTRODUCE NEW PSEUDO RANDOM NUMBER GENERATORS HERE */

		default:
			printf("Error in invokeTestSuite!\n");
			break;
	}
	printf("     Statistical Testing Complete!!!!!!!!!!!!\n\n");
}


void
nist_test_suite()
{
	if ( (testVector[0] == 1) || (testVector[TEST_FREQUENCY] == 1) )
		Frequency(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_BLOCK_FREQUENCY] == 1) )
		BlockFrequency(tp.blockFrequencyBlockLength, tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_CUSUM] == 1) )
		CumulativeSums(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_RUNS] == 1) )
		Runs(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_LONGEST_RUN] == 1) )
		LongestRunOfOnes(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_RANK] == 1) )
		Rank(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_FFT] == 1) )
		DiscreteFourierTransform(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_NONPERIODIC] == 1) )
		NonOverlappingTemplateMatchings(tp.nonOverlappingTemplateBlockLength, tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_OVERLAPPING] == 1) )
		OverlappingTemplateMatchings(tp.overlappingTemplateBlockLength, tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_UNIVERSAL] == 1) )
		Universal(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_APEN] == 1) )
		ApproximateEntropy(tp.approximateEntropyBlockLength, tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION] == 1) )
		RandomExcursions(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION_VAR] == 1) )
		RandomExcursionsVariant(tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_SERIAL] == 1) )
		Serial(tp.serialBlockLength,tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_LINEARCOMPLEXITY] == 1) )
		LinearComplexity(tp.linearComplexitySequenceLength, tp.n);

	if ( (testVector[0] == 1) || (testVector[TEST_GAMBLER] == 1) )
		Gambler(tp.gamblerRunsPerStartingPoint, tp.n, tp.gamblerStartStartingPoint,
			      tp.gamblerEndStartingPoint, tp.gamblerTestMask);
}
