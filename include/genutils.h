
#ifndef _GENUTILS_H_
#define _GENUTILS_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"

typedef struct _MP_struct {
	int		size;	/*  in bytes  */
	int		bitlen;	/*  in bits, duh  */
	uint8_t	*val;
	} MP;

#define	FREE(A)	if ( (A) ) { free((A)); (A) = NULL; }
#define	ASCII2BIN(ch)	( (((ch) >= '0') && ((ch) <= '9')) ? ((ch) - '0') : (((ch) >= 'A') && ((ch) <= 'F')) ? ((ch) - 'A' + 10) : ((ch) - 'a' + 10) )

#ifndef EXPWD
#define	EXPWD		((DBLWORD)1<<NUMLEN)
#endif

#define	sniff_bit(ptr,mask)		(*(ptr) & mask)

/*
 * Function Declarations
 */
int		greater(uint8_t *x, uint8_t *y, int l);
int		less(uint8_t *x, uint8_t *y, int l);
uint8_t	bshl(uint8_t *x, int l);
void	bshr(uint8_t *x, int l);
int		Mult(uint8_t *A, uint8_t *B, int LB, uint8_t *C, int LC);
void	ModSqr(uint8_t *A, uint8_t *B, int LB, uint8_t *M, int LM);
void	ModMult(uint8_t *A, uint8_t *B, int LB, uint8_t *C, int LC, uint8_t *M, int LM);
void	smult(uint8_t *A, uint8_t b, uint8_t *C, int L);
void	Square(uint8_t *A, uint8_t *B, int L);
void	ModExp(uint8_t *A, uint8_t *B, int LB, uint8_t *C, int LC, uint8_t *M, int LM);
int		DivMod(uint8_t *x, int lenx, uint8_t *n, int lenn, uint8_t *quot, uint8_t *rem);
void	Mod(uint8_t *x, int lenx, uint8_t *n, int lenn);
void	Div(uint8_t *x, int lenx, uint8_t *n, int lenn);
void	sub(uint8_t *A, int LA, uint8_t *B, int LB);
int		negate(uint8_t *A, int L);
uint8_t	add(uint8_t *A, int LA, uint8_t *B, int LB);
void	prettyprintBstr(char *S, uint8_t *A, int L);
void	byteReverse(uint32_t *buffer, int byteCount);
void	ahtopb (char *ascii_hex, uint8_t *p_binary, int bin_len);

#endif  /* _GENUTILS_H_ */
