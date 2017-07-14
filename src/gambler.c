#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include "../include/cephes.h"
#include "../include/externs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                            G A M B L E R  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

 // Number of states/total capital
#define N 33

#ifndef USE_BIT_TRACKER
#define USE_BIT_TRACKER 1
#endif

typedef struct
{
	BitSequence* seq;
	size_t pos;
	size_t max;
	unsigned wraps;
} stream_state;

bool
next_bit(stream_state* state)
{
	bool bit = state->seq[state->pos++];

	//fprintf(stderr, "%d", bit);

	if ( state->pos == state->max )
	{
		state->pos = 0;
		state->wraps++;
		//fprintf(stderr, "Gambler next_bit wrap no. %d\n", state->wraps);
	}

	return bit;
}

#if USE_BIT_TRACKER == 0
/* Read 64 bits from the stream_state as 64-bit little endian, unsigned integer
 */
uint64_t
next_u64(stream_state* state)
{
	uint64_t s = 0; uint64_t m = 1;
	int c;
	#pragma unroll 64
	for ( c = 0; c < 64; ++c )
	{
		if ( next_bit(state) )
			s |= m;
		m <<= 1;
	}
	return s;
}
#else
/* A simplified version of BitTracker
   will only split on 1 value instead on any number of them
 */
bool
bit_track(uint64_t split, stream_state* state)
{
	uint64_t mask = 0x8000000000000000ull;
	uint64_t left = 0x0000000000000000ull;
	uint64_t right = 0xffffffffffffffffull;
	do {
		bool bit = next_bit(state);
		if ( bit )
			left |= mask;
		else
			right &= ~mask;

		mask >>= 1;
	} while ( left < split && split < right && mask );
	if ( !mask )
		fprintf(stderr, "bit_track aborted abnormally\n");

	return left >= split;
}
#endif

uint64_t
frac_to_fixed(uint64_t nom, uint64_t den)
{
	uint64_t r = 0;
	for ( unsigned i=64; i > 0; --i )
	{
		nom <<= 1;
		r <<= 1;
		if ( nom >= den )
		{
			nom -= den;
			r |= 1;
		}
	}
	return r;
}

// Transition probabilities
uint64_t T1[N]; // = 0.48
uint64_t T2[N]; // = i / (2i+1)
uint64_t T3[N]; // = (i^3)/(i^3+(i+1)^3)

// Expected winning outcomes
double T1e[N];
double T2e[N];
double T3e[N];

// Expected game lengths
double T1le[N];
double T2le[N];
double T3le[N];

// Expected game lengths
double T1lv[N];
double T2lv[N];
double T3lv[N];

void
compute_times(uint64_t *prob, double *e_out, double *var_out)
{
	#define p2d(P) ((double)(P) / (double)(1ull << 32) / (double)(1ull << 32))

	double p, q;
	double d[N]; // d[k] = prod[i=1, k] { q(i) / p(i) }
	double f[N]; // f[k] = sum[i=1, k] { 1 / (p(i)d[i]) }
	double g[N]; // g[m] = sum[k=1, m-1] { d[k]f[k] }
	double h[N]; // h[m] = sum[k=0, m-1] { d[k] }

	double z[N]; // z[k] = 2(ET[k] - 1) / q(k)
	double x[N]; // x[m] = sum[k=1, m-1] { z[k] / d[k-1] }
	double c[N]; // c[m] = sum[k=1, m-1] { z[k] / d[k-1] * h[k] }
	double v[N]; // v[m] = x[m] * h[m] - c[m]

	unsigned i;

	d[0] = 1;
	f[0] = 0;
	g[0] = 0;
	h[0] = 0;
	z[0] = 0;
	x[0] = 0;
	c[0] = 0;
	v[0] = 0;
	for ( i = 1; i < N; ++i )
	{
		p = p2d(prob[i]);
		q = p2d(-prob[i]);
		d[i] = d[i-1] / p * q;
		f[i] = f[i-1] + (1 / d[i-1] / q);
		g[i] = g[i-1] + d[i-1] * f[i-1];
		h[i] = h[i-1] + d[i-1];
	}
	// Wrap to accomodate for non-existing g[N] and h[N]
	g[0] = g[i-1] + d[i-1] * f[i-1];
	h[0] = h[i-1] + d[i-1];

	for ( i = 1; i < N; ++i )
	{
		e_out[i] = g[0] / h[0] * h[i] - g[i];
		z[i] = 2 * (e_out[i] - 1) / p2d(-prob[i]);
		// the variables below have slightly different indexing han g and h, i.e. v[N-1] corresponds to the v[N] from python code.
		x[i] = x[i-1] + z[i] / d[i-1];
		c[i] = c[i-1] + z[i] / d[i-1] * h[i];
		v[i] = x[i] * h[(i+1)%N] - c[i];
	}
	for ( i = 1; i < N; ++i )
	{
		var_out[i] = (g[0] + v[N-1]) / h[0] * h[i] - g[i] - v[i-1] - e_out[i]*e_out[i];
	}
	#undef p2d
}

void
init_arrays ()
{
	static bool done = false;
	if ( done ) return;

	// For all possible points i:
	for ( unsigned i = 1; i < N; ++i )
	{
		// Compute different fracs for winning probabilities
		T1[i] = frac_to_fixed(48, 100);
		T2[i] = frac_to_fixed(i, 2*i + 1);
		T3[i] = frac_to_fixed(i*i*i, i*i*i + (i+1)*(i+1)*(i+1));
	}

	// For all starting points, compute expected prob.
	// T1
	{
		double p = 0.48;
		double q = 1-p;
		double r = q / p;

		double d = (1 - pow(r, N));
		for ( unsigned s = 1; s < N; ++s )
		{
			T1e[s] = (1 - pow(r, s)) / d;
		}

		compute_times(T1, T1le, T1lv);
	}
	// T2
	{
		double d = N * (N + 1);
		for ( unsigned s = 1; s < N; ++s )
		{
			T2e[s] = (double)(s * (s+1)) / d;
		}

		compute_times(T2, T2le, T2lv);
	}
	// T3
	{
		double d = N * N * (N + 1) * (N + 1);
		for ( unsigned s = 1; s < N; ++s )
		{
			T3e[s] = (double)(s * s * (s+1) * (s+1)) / d;
		}

		compute_times(T3, T3le, T3lv);
	}

	done = true;
}

void
run_gambler(stream_state* st, unsigned s, uint64_t* p, uint64_t* w_out, uint64_t* l_out)
{
	uint64_t l = 0;
	unsigned i = s;

	while ( i > 0 && i < N )
	{
		#if USE_BIT_TRACKER == 1
		// Negate the result of bit_track : result 0 suggests value < p;
		bool step = !bit_track(p[i], st);
		#else
		bool step = (next_u64(st) < p[i]);
		#endif

		if ( step )
		{
			++i;
		}
		else
		{
			--i;
		}
		++l;
	}

	//fprintf(stderr, "Game from starting point %d was %s. Game length: %zd\n", s, i == N ? "won" : "lost", l);

	if ( i == N )
		++(*w_out);

	(*l_out) += l;
}

// Total wins per starting point
uint64_t W1[N];
uint64_t W2[N];
uint64_t W3[N];

// Total game lengths per stering point
uint64_t L1[N];
uint64_t L2[N];
uint64_t L3[N];

void
zero_arrays()
{
	for ( unsigned i=0; i < N; ++i )
	{
		W1[i] = 0;
		W2[i] = 0;
		W3[i] = 0;
		L1[i] = 0;
		L2[i] = 0;
		L3[i] = 0;
	}
}

void
print_statistics(unsigned M, double *EXs, uint64_t *Xs, double *ETs, double *VTs, uint64_t *Ts)
{
	// For each starting point s
	unsigned s;
	double X, T, Y, sigma2, cdf, cdfc, p_value;

	fprintf(stats[TEST_GAMBLER], "\t\tWin probabilities:\n");
	// Win probability
	for ( s = 1; s < N; ++s )
	{
		// Compute variance of the distribution and the actual variable Y
		sigma2 = EXs[s] * (1 - EXs[s]);
		X = (double)(Xs[s]);
		Y = (X - M * EXs[s]) / sqrt(M * sigma2);
		// Compute CDF and P-Values
		cdf = cephes_normal(Y);
		cdfc = 1 - cdf;
		p_value = cdf < cdfc ? cdf : cdfc;

		fprintf(stats[TEST_GAMBLER], "%s\t\tY = %f\tp_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", Y, p_value); fflush(stats[TEST_GAMBLER]);
		fprintf(results[TEST_GAMBLER], "%f\n", p_value); fflush(results[TEST_GAMBLER]);
	}

	fprintf(stats[TEST_GAMBLER], "\t\tGame durations:\n");
	// Game duration
	for ( s = 1; s < N; ++s )
	{
		// Compute the actual variable Y
		sigma2 = VTs[s];
		T = (double)(Ts[s]);
		Y = (T - M * ETs[s]) / sqrt(M * sigma2);
		// Compute CDF and P-Values
		cdf = cephes_normal(Y);
		cdfc = 1 - cdf;
		p_value = cdf < cdfc ? cdf : cdfc;

		fprintf(stats[TEST_GAMBLER], "%s\t\tY = %f\tp_value = %f\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", Y, p_value); fflush(stats[TEST_GAMBLER]);
		fprintf(results[TEST_GAMBLER], "%f\n", p_value); fflush(results[TEST_GAMBLER]);
	}
}


void
Gambler(int M, int n)
{
	stream_state st = { epsilon, 0, n, 0 };
	unsigned s, c;

	init_arrays();
	zero_arrays();

	// For each starting point, run gambler simulator for stream_state and T table
	for ( s = 1; s < N; ++s )
	{
		for ( c = M; c > 0; --c )
		{
			run_gambler(&st, s, T1, &W1[s], &L1[s]);
			run_gambler(&st, s, T2, &W2[s], &L2[s]);
			run_gambler(&st, s, T3, &W3[s], &L3[s]);
		}
	}

	if ( st.wraps )
		fprintf(stderr, "WARNING: Gambler next_bit wrapped %d times.\n", st.wraps);

	fprintf(stats[TEST_GAMBLER], "\t\t\t        GAMBLER TEST\n");
	fprintf(stats[TEST_GAMBLER], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_GAMBLER], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_GAMBLER], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_GAMBLER], "\t\tNo. of next_bit wraps: %d                    \n", st.wraps);
	fprintf(stats[TEST_GAMBLER], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_GAMBLER], "\t\tT1 :                                         \n");
	for ( s = 1; s < N; ++s ) {
		fprintf(stats[TEST_GAMBLER], "\t\t\tStarting point %d:                         \n", s);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(a) Games won     (/M) (expected) = %ld (%lf) (%lf)\n", W1[s], (double)W1[s]/M, T1e[s]);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(b) Total length  (/M) (expected) = %ld (%lf) (%lf)\n", L1[s], (double)L1[s]/M, T1le[s]);
	}
	fprintf(stats[TEST_GAMBLER], "\t\tT2 :                                         \n");
	for ( s = 1; s < N; ++s ) {
		fprintf(stats[TEST_GAMBLER], "\t\t\tStarting point %d:                         \n", s);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(a) Games won     (/M) (expected) = %ld (%lf) (%lf)\n", W2[s], (double)W2[s]/M, T2e[s]);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(b) Total length  (/M) (expected) = %ld (%lf) (%lf)\n", L2[s], (double)L2[s]/M, T2le[s]);
	}
	fprintf(stats[TEST_GAMBLER], "\t\tT3 :                                         \n");
	for ( s = 1; s < N; ++s ) {
		fprintf(stats[TEST_GAMBLER], "\t\t\tStarting point %d:                         \n", s);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(a) Games won     (/M) (expected) = %ld (%lf) (%lf)\n", W3[s], (double)W3[s]/M, T3e[s]);
		fprintf(stats[TEST_GAMBLER], "\t\t\t\t(b) Total length  (/M) (expected) = %ld (%lf) (%lf)\n", L3[s], (double)L3[s]/M, T3le[s]);
	}
	fprintf(stats[TEST_GAMBLER], "\t\t---------------------------------------------\n");

	print_statistics(M, T1e, W1, T1le, T1lv, L1);
	print_statistics(M, T2e, W2, T2le, T2lv, L2);
	print_statistics(M, T3e, W3, T3le, T3lv, L3);
}
