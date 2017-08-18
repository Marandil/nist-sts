#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <set>
#include <random>

constexpr uint64_t pow2m1(int n) { return (1ull << n) - 1; }

class PRNG
{
public:
    virtual void setSeed(uint64_t seed) = 0;
    virtual uint64_t nextInt() = 0;
    virtual uint32_t getNrOfBits() = 0;
};

class LCG : public PRNG
{
public:
    const uint64_t M;
    const uint64_t a, b;
    const uint32_t nrOfBits;
    
    LCG(uint64_t M_, uint64_t a_, uint64_t b_, const uint32_t nrOfBits_)
        : M(M_)
        , a(a_)
        , b(b_)
        , nrOfBits(nrOfBits_)
    {
    }
    
    void setSeed(uint64_t seed)
    {
        this->seed = seed;
    }
    
    uint64_t nextInt()
    {
        seed = (a * seed + b) % M;
        return seed & pow2m1(nrOfBits);
    }
    
    uint32_t getNrOfBits()
    {
        return nrOfBits;
    }
    
private:
    uint64_t seed = 1;
};

class SomeBits : public PRNG
{
public:
    SomeBits(std::shared_ptr<PRNG> prng_, uint32_t mostSig_, uint32_t leastSig_)
        : prng(prng_)
        , mostSig(mostSig_)
        , leastSig(leastSig_)
        , nrOfBits(mostSig_ + 1 - leastSig_)
    {
    }
    
    void setSeed(uint64_t seed)
    {
        prng->setSeed(seed);
    }
    
    uint64_t nextInt()
    {
        return getBits(prng->nextInt());
    }
    
    uint32_t getNrOfBits()
    {
        return nrOfBits;
    }
    
private:
    const std::shared_ptr<PRNG> prng;
    const uint32_t mostSig;
    const uint32_t leastSig;
    const uint32_t nrOfBits;
    
    uint64_t getBits(uint64_t n)
    {
        n = n >> leastSig;
        return n & pow2m1(nrOfBits);
    }
};

class CMRG : public PRNG
{
public:
    CMRG()
    {
        reset();
    }
    
    void setSeed(uint64_t seed)
    {
        x[0] = 0;
        x[1] = seed;
        x[2] = 0;
        y[0] = 0;
        y[1] = 0;
        y[2] = seed;
        n = 0;
    }
    
    uint64_t nextInt()
    {
        int64_t nextx = xa *  x[(n+1)%3]  -  xb * x[n];
        int64_t nexty = ya *  y[(n+2)%3]  -  yb * y[n];
        nextx = mymod(nextx, xm);
        nexty = mymod(nexty, ym);
        x[n] = nextx;
        y[n] = nexty;
        n = (n + 1)%3;
        return static_cast<uint64_t>( (zm + nextx - nexty) % zm );
    }
    
    uint32_t getNrOfBits()
    {
        return 31;
    }
    
private:
    int64_t x[3], y[3];
    const int64_t xa = 63308;
    const int64_t xb = 183326;
    const int64_t xm = 2147483647LL;
    const int64_t ya = 86098;
    const int64_t yb = 539608;
    const int64_t ym = 2145483479LL;
    const int64_t zm = 2147483647LL;
    int n;
    
    void reset()
    {
        setSeed(1);
    }
    
    int64_t mymod(int64_t a, int64_t m)
    {
        if (m == 0)
            return a;
        if (m < 0)
            return mymod(-a, -m);
        int64_t res = a % m;
        if (res < 0)
            res += m;
        return res;
    }
};

class C_PRG : public PRNG
{
public:
    void setSeed(uint64_t seed)
    {
        srand(seed);
    }
    
    uint64_t nextInt()
    {
        return static_cast<uint64_t>(rand());
    }
    
    uint32_t getNrOfBits()
    {
        return 31;
    }
};

class BorlandPRNG : public PRNG
{
    uint64_t nextInt()
    {
        myseed = myseed * 0x015A4E35 + 1;
        return static_cast<uint64_t>( (myseed >> 16) & 0x7FFF );
    }
    
    void setSeed(uint64_t seed)
    {
        myseed = seed;
        //myrand();
    }
    
    uint32_t getNrOfBits()
    {
        return 15;
    }
    
private:
    uint32_t myseed = 0x015A4E36;
};

class VisualPRNG : public PRNG
{
    uint64_t nextInt()
    {
        myseed = myseed * 0x343FDu + 0x269EC3u;
        return static_cast<uint64_t>( (myseed >> 16) & 0x7FFF );
    }
    
    void setSeed(uint64_t seed)
    {
        myseed = seed;
        //myrand();
    }
    
    uint32_t getNrOfBits()
    {
        return 15;
    }
    
private:
    uint32_t myseed = 1;
};

class Mersenne : public PRNG
{
public:
    void setSeed(uint64_t seed)
    {
        eng.seed(seed);
    }
    
    uint64_t nextInt()
    {
//         //uint64_t r = static_cast<uint64_t>(eng());
        //fprintf(stderr, "%llu\n", r & pow2m1[63]);
        //return r & pow2m1[63];
        return static_cast<uint64_t>(eng());
    }
    
    uint32_t getNrOfBits()
    {
        return 64;
    }
    
    std::mt19937_64 eng;
};

class RandU : public PRNG
{
public:
    void setSeed(uint64_t seed)
    {
//        s = seed + (seed % 2 == 0 ? 1 : 0);
		s = seed;
    }
    
    uint64_t nextInt()
    {
        s = (65539llu * s) & pow2m1(31);
        return s;
    }
    
    uint32_t getNrOfBits()
    {
        return 31;
    }
    
    uint64_t s;
};

std::shared_ptr<PRNG> getPRNG(char* name)
{
	if (strcmp(name, "Rand") == 0)
    {
        return std::make_shared<LCG>(2147483648LL, 1103515245, 12345, 31);
    }
    else if (strcmp(name, "Rand0") == 0)
    {
        return std::make_shared<LCG>(2147483648LL, 1103515245, 12345, 8);
    }
    else if (strcmp(name, "Rand1") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<LCG>(2147483648LL, 1103515245, 12345, 31),
                15, 8
            );
    }
    else if (strcmp(name, "Rand3") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<LCG>(2147483648LL, 1103515245, 12345, 31),
                30, 23
            );
    }
    else if (strcmp(name, "Minstd") == 0)
    {
        return std::make_shared<LCG>(2147483647, 16807, 0, 31);
    }
    else if (strcmp(name, "Minstd0") == 0)
    {
        return std::make_shared<LCG>(2147483647, 16807, 0, 8);
    }
    else if (strcmp(name, "Minstd1") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<LCG>(2147483647, 16807, 0, 31),
                15, 8
            );
    }
    else if (strcmp(name, "NewMinstd") == 0)
    {
        return std::make_shared<LCG>(2147483647, 48271, 0, 31);
    }
    else if (strcmp(name, "NewMinstd0") == 0)
    {
        return std::make_shared<LCG>(2147483647, 48271, 0, 8);
    }
    else if (strcmp(name, "NewMinstd1") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<LCG>(2147483647, 48271, 0, 31),
                15, 8
            );
    }
    else if (strcmp(name, "NewMinstd3") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<LCG>(2147483647, 48271, 0, 31),
                30, 23
            );
    }
    else if (strcmp(name, "CMRG") == 0)
    {
        return std::make_shared<CMRG>();
    }
    else if (strcmp(name, "CMRG0") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<CMRG>(),
                7, 0
            );
    }
    else if (strcmp(name, "CMRG1") == 0)
    {
        return std::make_shared<SomeBits>(
                std::make_shared<CMRG>(),
                15, 8
            );
    }
    else if (strcmp(name, "Borland") == 0)
    {
        return std::make_shared<BorlandPRNG>();
    }
    else if (strcmp(name, "C_PRG") == 0)
    {
        return std::make_shared<C_PRG>();
    }
    else if (strcmp(name, "VS") == 0)
    {
        return std::make_shared<VisualPRNG>();
    }
    else if (strcmp(name, "Mersenne") == 0)
    {
        return std::make_shared<Mersenne>();
    }
    else if (strcmp(name, "RANDU") == 0)
    {
        return std::make_shared<RandU>();
    }
	return nullptr;
}



uint64_t curr;
int filled;

uint64_t nextChunk(std::shared_ptr<PRNG> prng)
{
    uint64_t r = 0;
    int nrOfBits = prng->getNrOfBits();
    while (filled < 64)
    {
        r = prng->nextInt();
        curr += (r << filled);
        filled += nrOfBits;
    }
    int used = nrOfBits + 64 - filled;
    uint64_t res = curr;
    curr = used < 64 ? (r >> used) : 0;
    filled = nrOfBits - used;
    return res;
}

void generateString(uint64_t nrOfBits, std::shared_ptr<PRNG> prng)
{
    uint64_t nrOfChunks = nrOfBits / 64;
    curr = 0;
    filled = 0;
    for (uint64_t i = 0; i < nrOfChunks; ++i)
    {
        uint64_t chunk = nextChunk(prng);
        //fprintf(stderr, "filled = %d\n", filled);
		static_assert(sizeof(uint64_t) == 8, "uint64_t not 64 bit long!");
        fwrite(&chunk, sizeof(uint64_t), 1, stdout);
    }
}

int main(int argc, char** argv)
{   
    if (argc < 4)
	{
		printf("Usage: %s [prng name] [length] [seed]\n", argv[0]);
        exit(0);
	}
    
    std::shared_ptr<PRNG> prng = getPRNG(argv[1]);
    if (!prng)
    {
        printf("Unknown prng: %s\n", argv[1]);
        exit(0);
    }

    char *pEnd;
	uint64_t length = strtoull(argv[2], &pEnd, 10);
	uint64_t seed = strtoull(argv[3], &pEnd, 10);
	
	prng->setSeed(seed);
	generateString(length * 8, prng);
}

