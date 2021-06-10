#ifndef RIPEMD160_CL
#define RIPEMD160_CL

#ifndef endian
#define endian(x) ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((x) >> 24)
#endif

__constant unsigned int RIPEMD160_IV[5] = {
    0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0,
};

__constant unsigned int K[8] = {
    0x5a827999, 0x6ed9eba1, 0x8f1bbcdc, 0xa953fd4e, 0x7a6d76e9, 0x6d703ef3, 0x5c4dd124, 0x50a28be6
};

#define rotl(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

#define F(x, y, z) ((x) ^ (y) ^ (z))

#define G(x, y, z) (((x) & (y)) | (~(x) & (z)))

#define H(x, y, z) (((x) | ~(y)) ^ (z))

#define I(x, y, z) (((x) & (z)) | ((y) & ~(z)))

#define J(x, y, z) ((x) ^ ((y) | ~(z)))

#define FF(a, b, c, d, e, m, s)\
    a += (F((b), (c), (d)) + (m));\
    a = (rotl((a), (s)) + (e));\
    c = rotl((c), 10)

#define GG(a, b, c, d, e, x, s)\
    a += G((b), (c), (d)) + (x) + K[0];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define HH(a, b, c, d, e, x, s)\
    a += H((b), (c), (d)) + (x) + K[1];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define II(a, b, c, d, e, x, s)\
    a += I((b), (c), (d)) + (x) + K[2];\
    a = rotl((a), (s)) + e;\
    c = rotl((c), 10)

#define JJ(a, b, c, d, e, x, s)\
    a += J((b), (c), (d)) + (x) + K[3];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define FFF(a, b, c, d, e, x, s)\
    a += F((b), (c), (d)) + (x);\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define GGG(a, b, c, d, e, x, s)\
    a += G((b), (c), (d)) + x + K[4];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define HHH(a, b, c, d, e, x, s)\
    a += H((b), (c), (d)) + (x) + K[5];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define III(a, b, c, d, e, x, s)\
    a += I((b), (c), (d)) + (x) + K[6];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

#define JJJ(a, b, c, d, e, x, s)\
    a += J((b), (c), (d)) + (x) + K[7];\
    a = rotl((a), (s)) + (e);\
    c = rotl((c), 10)

void ripemd160p1(__private const unsigned int x[8], unsigned int digest[5])
{
    __private unsigned int a = RIPEMD160_IV[0];
    __private unsigned int b = RIPEMD160_IV[1];
    __private unsigned int c = RIPEMD160_IV[2];
    __private unsigned int d = RIPEMD160_IV[3];
    __private unsigned int e = RIPEMD160_IV[4];

    /* round 1 */
    FF(a, b, c, d, e, x[0], 11);
    FF(e, a, b, c, d, x[1], 14);
    FF(d, e, a, b, c, x[2], 15);
    FF(c, d, e, a, b, x[3], 12);
    FF(b, c, d, e, a, x[4], 5);
    FF(a, b, c, d, e, x[5], 8);
    FF(e, a, b, c, d, x[6], 7);
    FF(d, e, a, b, c, x[7], 9);
    FF(c, d, e, a, b, 128, 11);
    FF(b, c, d, e, a, 0, 13);
    FF(a, b, c, d, e, 0, 14);
    FF(e, a, b, c, d, 0, 15);
    FF(d, e, a, b, c, 0, 6);
    FF(c, d, e, a, b, 0, 7);
    FF(b, c, d, e, a, 256, 9);
    FF(a, b, c, d, e, 0, 8);

    /* round 2 */
    GG(e, a, b, c, d, x[7], 7);
    GG(d, e, a, b, c, x[4], 6);
    GG(c, d, e, a, b, 0, 8);
    GG(b, c, d, e, a, x[1], 13);
    GG(a, b, c, d, e, 0, 11);
    GG(e, a, b, c, d, x[6], 9);
    GG(d, e, a, b, c, 0, 7);
    GG(c, d, e, a, b, x[3], 15);
    GG(b, c, d, e, a, 0, 7);
    GG(a, b, c, d, e, x[0], 12);
    GG(e, a, b, c, d, 0, 15);
    GG(d, e, a, b, c, x[5], 9);
    GG(c, d, e, a, b, x[2], 11);
    GG(b, c, d, e, a, 256, 7);
    GG(a, b, c, d, e, 0, 13);
    GG(e, a, b, c, d, 0x80, 12);

    /* round 3 */
    HH(d, e, a, b, c, x[3], 11);
    HH(c, d, e, a, b, 0, 13);
    HH(b, c, d, e, a, 256, 6);
    HH(a, b, c, d, e, x[4], 7);
    HH(e, a, b, c, d, 0, 14);
    HH(d, e, a, b, c, 0, 9);
    HH(c, d, e, a, b, 0x80, 13);
    HH(b, c, d, e, a, x[1], 15);
    HH(a, b, c, d, e, x[2], 14);
    HH(e, a, b, c, d, x[7], 8);
    HH(d, e, a, b, c, x[0], 13);
    HH(c, d, e, a, b, x[6], 6);
    HH(b, c, d, e, a, 0, 5);
    HH(a, b, c, d, e, 0, 12);
    HH(e, a, b, c, d, x[5], 7);
    HH(d, e, a, b, c, 0, 5);

    /* round 4 */
    II(c, d, e, a, b, x[1], 11);
    II(b, c, d, e, a, 0, 12);
    II(a, b, c, d, e, 0, 14);
    II(e, a, b, c, d, 0, 15);
    II(d, e, a, b, c, x[0], 14);
    II(c, d, e, a, b, 0x80, 15);
    II(b, c, d, e, a, 0, 9);
    II(a, b, c, d, e, x[4], 8);
    II(e, a, b, c, d, 0, 9);
    II(d, e, a, b, c, x[3], 14);
    II(c, d, e, a, b, x[7], 5);
    II(b, c, d, e, a, 0, 6);
    II(a, b, c, d, e, 256, 8);
    II(e, a, b, c, d, x[5], 6);
    II(d, e, a, b, c, x[6], 5);
    II(c, d, e, a, b, x[2], 12);

    /* round 5 */
    JJ(b, c, d, e, a, x[4], 9);
    JJ(a, b, c, d, e, x[0], 15);
    JJ(e, a, b, c, d, x[5], 5);
    JJ(d, e, a, b, c, 0, 11);
    JJ(c, d, e, a, b, x[7], 6);
    JJ(b, c, d, e, a, 0, 8);
    JJ(a, b, c, d, e, x[2], 13);
    JJ(e, a, b, c, d, 0, 12);
    JJ(d, e, a, b, c, 256, 5);
    JJ(c, d, e, a, b, x[1], 12);
    JJ(b, c, d, e, a, x[3], 13);
    JJ(a, b, c, d, e, 0x80, 14);
    JJ(e, a, b, c, d, 0, 11);
    JJ(d, e, a, b, c, x[6], 8);
    JJ(c, d, e, a, b, 0, 5);
    JJ(b, c, d, e, a, 0, 6);

    digest[0] = c;
    digest[1] = d;
    digest[2] = e;
    digest[3] = a;
    digest[4] = b;
}

void ripemd160p2(const unsigned int x[8], unsigned int digest[5])
{
    __private unsigned int a = RIPEMD160_IV[0];
    __private unsigned int b = RIPEMD160_IV[1];
    __private unsigned int c = RIPEMD160_IV[2];
    __private unsigned int d = RIPEMD160_IV[3];
    __private unsigned int e = RIPEMD160_IV[4];

    /* parallel round 1 */
    JJJ(a, b, c, d, e, x[5], 8);
    JJJ(e, a, b, c, d, 256, 9);
    JJJ(d, e, a, b, c, x[7], 9);
    JJJ(c, d, e, a, b, x[0], 11);
    JJJ(b, c, d, e, a, 0, 13);
    JJJ(a, b, c, d, e, x[2], 15);
    JJJ(e, a, b, c, d, 0, 15);
    JJJ(d, e, a, b, c, x[4], 5);
    JJJ(c, d, e, a, b, 0, 7);
    JJJ(b, c, d, e, a, x[6], 7);
    JJJ(a, b, c, d, e, 0, 8);
    JJJ(e, a, b, c, d, 0x80, 11);
    JJJ(d, e, a, b, c, x[1], 14);
    JJJ(c, d, e, a, b, 0, 14);
    JJJ(b, c, d, e, a, x[3], 12);
    JJJ(a, b, c, d, e, 0, 6);

    /* parallel round 2 */
    III(e, a, b, c, d, x[6], 9);
    III(d, e, a, b, c, 0, 13);
    III(c, d, e, a, b, x[3], 15);
    III(b, c, d, e, a, x[7], 7);
    III(a, b, c, d, e, x[0], 12);
    III(e, a, b, c, d, 0, 8);
    III(d, e, a, b, c, x[5], 9);
    III(c, d, e, a, b, 0, 11);
    III(b, c, d, e, a, 256, 7);
    III(a, b, c, d, e, 0, 7);
    III(e, a, b, c, d, 0x80, 12);
    III(d, e, a, b, c, 0, 7);
    III(c, d, e, a, b, x[4], 6);
    III(b, c, d, e, a, 0, 15);
    III(a, b, c, d, e, x[1], 13);
    III(e, a, b, c, d, x[2], 11);

    /* parallel round 3 */
    HHH(d, e, a, b, c, 0, 9);
    HHH(c, d, e, a, b, x[5], 7);
    HHH(b, c, d, e, a, x[1], 15);
    HHH(a, b, c, d, e, x[3], 11);
    HHH(e, a, b, c, d, x[7], 8);
    HHH(d, e, a, b, c, 256, 6);
    HHH(c, d, e, a, b, x[6], 6);
    HHH(b, c, d, e, a, 0, 14);
    HHH(a, b, c, d, e, 0, 12);
    HHH(e, a, b, c, d, 0x80, 13);
    HHH(d, e, a, b, c, 0, 5);
    HHH(c, d, e, a, b, x[2], 14);
    HHH(b, c, d, e, a, 0, 13);
    HHH(a, b, c, d, e, x[0], 13);
    HHH(e, a, b, c, d, x[4], 7);
    HHH(d, e, a, b, c, 0, 5);

    /* parallel round 4 */
    GGG(c, d, e, a, b, 0x80, 15);
    GGG(b, c, d, e, a, x[6], 5);
    GGG(a, b, c, d, e, x[4], 8);
    GGG(e, a, b, c, d, x[1], 11);
    GGG(d, e, a, b, c, x[3], 14);
    GGG(c, d, e, a, b, 0, 14);
    GGG(b, c, d, e, a, 0, 6);
    GGG(a, b, c, d, e, x[0], 14);
    GGG(e, a, b, c, d, x[5], 6);
    GGG(d, e, a, b, c, 0, 9);
    GGG(c, d, e, a, b, x[2], 12);
    GGG(b, c, d, e, a, 0, 9);
    GGG(a, b, c, d, e, 0, 12);
    GGG(e, a, b, c, d, x[7], 5);
    GGG(d, e, a, b, c, 0, 15);
    GGG(c, d, e, a, b, 256, 8);

    /* parallel round 5 */
    FFF(b, c, d, e, a, 0, 8);
    FFF(a, b, c, d, e, 0, 5);
    FFF(e, a, b, c, d, 0, 12);
    FFF(d, e, a, b, c, x[4], 9);
    FFF(c, d, e, a, b, x[1], 12);
    FFF(b, c, d, e, a, x[5], 5);
    FFF(a, b, c, d, e, 0x80, 14);
    FFF(e, a, b, c, d, x[7], 6);
    FFF(d, e, a, b, c, x[6], 8);
    FFF(c, d, e, a, b, x[2], 13);
    FFF(b, c, d, e, a, 0, 6);
    FFF(a, b, c, d, e, 256, 5);
    FFF(e, a, b, c, d, x[0], 15);
    FFF(d, e, a, b, c, x[3], 13);
    FFF(c, d, e, a, b, 0, 11);
    FFF(b, c, d, e, a, 0, 11);

    digest[0] = d;
    digest[1] = e;
    digest[2] = a;
    digest[3] = b;
    digest[4] = c;
}

void ripemd160sha256NoFinal(const unsigned int x[8], unsigned int digest[5])
{
    __private unsigned int digest1[5];
    __private unsigned int digest2[5];

    ripemd160p1(x, digest1);
    ripemd160p2(x, digest2);

    digest[0] = digest1[0] + digest2[0];
    digest[1] = digest1[1] + digest2[1];
    digest[2] = digest1[2] + digest2[2];
    digest[3] = digest1[3] + digest2[3];
    digest[4] = digest1[4] + digest2[4];
}

void ripemd160FinalRound(const unsigned int hIn[5], unsigned int hOut[5])
{
    hOut[0] = endian(hIn[0] + RIPEMD160_IV[1]);
    hOut[1] = endian(hIn[1] + RIPEMD160_IV[2]);
    hOut[2] = endian(hIn[2] + RIPEMD160_IV[3]);
    hOut[3] = endian(hIn[3] + RIPEMD160_IV[4]);
    hOut[4] = endian(hIn[4] + RIPEMD160_IV[0]);
}

#endif
#ifndef SECP256K1_CL
#define SECP256K1_CL

typedef struct uint256_t {
    unsigned int v[8];
} uint256_t;

/**
 * Prime modulus 2^256 - 2^32 - 977
 */
__constant unsigned int P[8] = {
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFC2F
};

#ifdef DEVICE_VENDOR_INTEL
// Intel devices have a mul_hi bug
inline unsigned int mul_hi977(unsigned int x)
{
    unsigned int high = x >> 16;
    unsigned int low = x & 0xffff;

    return (((low * 977) >> 16) + (high * 977)) >> 16;
}

// 32 x 32 multiply-add
inline void madd977(unsigned int *high, unsigned int *low, unsigned int *a, unsigned int *c)
{
    *low = *a * 977;
    unsigned int tmp = *low + *c;
    unsigned int carry = tmp < *low ? 1 : 0;
    *low = tmp;
    *high = mul_hi977(*a) + carry;
}
#else

inline void madd977(unsigned int *high, unsigned int *low, unsigned int *a, unsigned int *c)
{
    *low = *a * 977;
    __private unsigned int tmp = *low + *c;
    __private unsigned int carry = tmp < *low ? 1 : 0;
    *low = tmp;
    *high = mad_hi(*a, (unsigned int)977, carry);
}

#endif

// Add with carry
#define addc(a, b, sum, carry, tmp)      \
    sum = (a) + (carry);                 \
    tmp = ((sum) < (a)) * 1;             \
    sum = (sum) + (b);                   \
    carry = (tmp) | (((sum) < (b)) * 1);

// subtract with borrow
#define subc(a, b, diff, borrow, tmp)    \
    tmp = (a) - (borrow);                \
    borrow = ((tmp) > (a)) * 1;          \
    diff = (tmp) - (b);                  \
    borrow |= ((diff) > (tmp)) ? 1 : 0;

#define add256k(a, b, c, carry, tmp)    \
    carry = 0;                          \
    addc(a[7], b[7], c[7], carry, tmp); \
    addc(a[6], b[6], c[6], carry, tmp); \
    addc(a[5], b[5], c[5], carry, tmp); \
    addc(a[4], b[4], c[4], carry, tmp); \
    addc(a[3], b[3], c[3], carry, tmp); \
    addc(a[2], b[2], c[2], carry, tmp); \
    addc(a[1], b[1], c[1], carry, tmp); \
    addc(a[0], b[0], c[0], carry, tmp);

#define sub256k( a, b, c, borrow, tmp)   \
    borrow = 0;                          \
    subc(a[7], b[7], c[7], borrow, tmp); \
    subc(a[6], b[6], c[6], borrow, tmp); \
    subc(a[5], b[5], c[5], borrow, tmp); \
    subc(a[4], b[4], c[4], borrow, tmp); \
    subc(a[3], b[3], c[3], borrow, tmp); \
    subc(a[2], b[2], c[2], borrow, tmp); \
    subc(a[1], b[1], c[1], borrow, tmp); \
    subc(a[0], b[0], c[0], borrow, tmp);

#define isInfinity256k(a)        \
    (                           \
        (a[0] == 0xffffffff) && \
        (a[1] == 0xffffffff) && \
        (a[2] == 0xffffffff) && \
        (a[3] == 0xffffffff) && \
        (a[4] == 0xffffffff) && \
        (a[5] == 0xffffffff) && \
        (a[6] == 0xffffffff) && \
        (a[7] == 0xffffffff)    \
    )

#define greaterOrEqualToP(a)    \
    (                           \
        (a[0] == 0xffffffff) && \
        (a[1] == 0xffffffff) && \
        (a[2] == 0xffffffff) && \
        (a[3] == 0xffffffff) && \
        (a[4] == 0xffffffff) && \
        (a[5] == 0xffffffff) && \
        (a[6] >= 0xfffffffe) && \
        (a[7] >= 0xfffffc2f)    \
    )

#define equal256k(a, b)   \
    (                     \
        (a[0] == b[0]) && \
        (a[1] == b[1]) && \
        (a[2] == b[2]) && \
        (a[3] == b[3]) && \
        (a[4] == b[4]) && \
        (a[5] == b[5]) && \
        (a[6] == b[6]) && \
        (a[7] == b[7])    \
    )

#define modP256k(a, carry, tmp)        \
    if (greaterOrEqualToP(a)) {        \
        sub256k(a, P, a, carry, tmp);  \
    }

void multiply256k(__private const unsigned int x[8], __private const unsigned int y[8], __private unsigned int out_high[8], __private unsigned int out_low[8])
{
    __private unsigned long product;

    // First round, overwrite z
    product = (unsigned long)x[7] * y[7];
    out_low[7] = (unsigned int)product;
    
    product = (unsigned long)x[7] * y[6] + (unsigned int)(product >> 32);
    out_low[6] = (unsigned int)product;
    
    product = (unsigned long)x[7] * y[5] + (unsigned int)(product >> 32);
    out_low[5] = (unsigned int)product;
    
    product = (unsigned long)x[7] * y[4] + (unsigned int)(product >> 32);
    out_low[4] = (unsigned int)product;
    
    product = (unsigned long)x[7] * y[3] + (unsigned int)(product >> 32);
    out_low[3] = (unsigned int)product;
    
    product = (unsigned long)x[7] * y[2] + (unsigned int)(product >> 32);
    out_low[2] = (unsigned int)product;
        
    product = (unsigned long)x[7] * y[1] + (unsigned int)(product >> 32);
    out_low[1] = (unsigned int)product;
        
    product = (unsigned long)x[7] * y[0] + (unsigned int)(product >> 32);
    out_low[0] = (unsigned int)product;
    out_high[7] = (unsigned int)(product >> 32);

    product = (unsigned long)x[6] * y[7] + out_low[6];
    out_low[6] = (unsigned int)product;

    /** round6 */
    product = (unsigned long)x[6] * y[6] + out_low[5] + (product >> 32);
    out_low[5] = (unsigned int)product;

    product = (unsigned long)x[6] * y[5] + out_low[4] + (product >> 32);
    out_low[4] = (unsigned int)product;

    product = (unsigned long)x[6] * y[4] + out_low[3] + (product >> 32);
    out_low[3] = (unsigned int)product;

    product = (unsigned long)x[6] * y[3] + out_low[2] + (product >> 32);
    out_low[2] = (unsigned int)product;

    product = (unsigned long)x[6] * y[2] + out_low[1] + (product >> 32);
    out_low[1] = (unsigned int)product;
    
    product = (unsigned long)x[6] * y[1] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;
    
    product = (unsigned long)x[6] * y[0] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;
    out_high[6] = product >> 32;

    /** round 5 */
    product = (unsigned long)x[5] * y[7] + out_low[5];
    out_low[5] = (unsigned int)product;

    product = (unsigned long)x[5] * y[6] + out_low[4] + (product >> 32);
    out_low[4] = (unsigned int)product;

    product = (unsigned long)x[5] * y[5] + out_low[3] + (product >> 32);
    out_low[3] = (unsigned int)product;

    product = (unsigned long)x[5] * y[4] + out_low[2] + (product >> 32);
    out_low[2] = (unsigned int)product;

    product = (unsigned long)x[5] * y[3] + out_low[1] + (product >> 32);
    out_low[1] = (unsigned int)product;

    product = (unsigned long)x[5] * y[2] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;
    
    product = (unsigned long)x[5] * y[1] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;
    
    product = (unsigned long)x[5] * y[0] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;
    out_high[5] = product >> 32;

    /** round 4 */
    product = (unsigned long)x[4] * y[7] + out_low[4];
    out_low[4] = (unsigned int)product;

    product = (unsigned long)x[4] * y[6] + out_low[3] + (product >> 32);
    out_low[3] = (unsigned int)product;

    product = (unsigned long)x[4] * y[5] + out_low[2] + (product >> 32);
    out_low[2] = (unsigned int)product;

    product = (unsigned long)x[4] * y[4] + out_low[1] + (product >> 32);
    out_low[1] = (unsigned int)product;

    product = (unsigned long)x[4] * y[3] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;

    product = (unsigned long)x[4] * y[2] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;
    
    product = (unsigned long)x[4] * y[1] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;
    
    product = (unsigned long)x[4] * y[0] + out_high[5] + (product >> 32);
    out_high[5] = (unsigned int)product;
    out_high[4] = product >> 32;

    /** round 3 */
    product = (unsigned long)x[3] * y[7] + out_low[3];
    out_low[3] = (unsigned int)product;

    product = (unsigned long)x[3] * y[6] + out_low[2] + (product >> 32);
    out_low[2] = (unsigned int)product;

    product = (unsigned long)x[3] * y[5] + out_low[1] + (product >> 32);
    out_low[1] = (unsigned int)product;

    product = (unsigned long)x[3] * y[4] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;

    product = (unsigned long)x[3] * y[3] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;

    product = (unsigned long)x[3] * y[2] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;
    
    product = (unsigned long)x[3] * y[1] + out_high[5] + (product >> 32);
    out_high[5] = (unsigned int)product;
    
    product = (unsigned long)x[3] * y[0] + out_high[4] + (product >> 32);
    out_high[4] = (unsigned int)product;
    out_high[3] = product >> 32;

    /** round 2 */
    product = (unsigned long)x[2] * y[7] + out_low[2];
    out_low[2] = (unsigned int)product;

    product = (unsigned long)x[2] * y[6] + out_low[1] + (product >> 32);
    out_low[1] = (unsigned int)product;

    product = (unsigned long)x[2] * y[5] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;

    product = (unsigned long)x[2] * y[4] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;

    product = (unsigned long)x[2] * y[3] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;

    product = (unsigned long)x[2] * y[2] + out_high[5] + (product >> 32);
    out_high[5] = (unsigned int)product;
    
    product = (unsigned long)x[2] * y[1] + out_high[4] + (product >> 32);
    out_high[4] = (unsigned int)product;
    
    product = (unsigned long)x[2] * y[0] + out_high[3] + (product >> 32);
    out_high[3] = (unsigned int)product;
    out_high[2] = product >> 32;
    
    /** round 1 */
    product = (unsigned long)x[1] * y[7] + out_low[1];
    out_low[1] = (unsigned int)product;

    product = (unsigned long)x[1] * y[6] + out_low[0] + (product >> 32);
    out_low[0] = (unsigned int)product;

    product = (unsigned long)x[1] * y[5] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;

    product = (unsigned long)x[1] * y[4] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;

    product = (unsigned long)x[1] * y[3] + out_high[5] + (product >> 32);
    out_high[5] = (unsigned int)product;

    product = (unsigned long)x[1] * y[2] + out_high[4] + (product >> 32);
    out_high[4] = (unsigned int)product;
    
    product = (unsigned long)x[1] * y[1] + out_high[3] + (product >> 32);
    out_high[3] = (unsigned int)product;
    
    product = (unsigned long)x[1] * y[0] + out_high[2] + (product >> 32);
    out_high[2] = (unsigned int)product;
    out_high[1] = product >> 32;

    /** round 0 */
    product = (unsigned long)x[0] * y[7] + out_low[0];
    out_low[0] = (unsigned int)product;

    product = (unsigned long)x[0] * y[6] + out_high[7] + (product >> 32);
    out_high[7] = (unsigned int)product;

    product = (unsigned long)x[0] * y[5] + out_high[6] + (product >> 32);
    out_high[6] = (unsigned int)product;

    product = (unsigned long)x[0] * y[4] + out_high[5] + (product >> 32);
    out_high[5] = (unsigned int)product;

    product = (unsigned long)x[0] * y[3] + out_high[4] + (product >> 32);
    out_high[4] = (unsigned int)product;

    product = (unsigned long)x[0] * y[2] + out_high[3] + (product >> 32);
    out_high[3] = (unsigned int)product;
    
    product = (unsigned long)x[0] * y[1] + out_high[2] + (product >> 32);
    out_high[2] = (unsigned int)product;
    
    product = (unsigned long)x[0] * y[0] + out_high[1] + (product >> 32);
    out_high[1] = (unsigned int)product;
    out_high[0] = product >> 32;
}

void mulModP256k(__private const unsigned int a[8], __private const unsigned int b[8], __private unsigned int product_low[8])
{
    __private unsigned int high[8];
    __private unsigned int low[8];

    __private unsigned int hWord = 0;
    __private unsigned int carry = 0;
    __private unsigned int t = 0;
    __private unsigned int product6 = 0;
    __private unsigned int product7 = 0;
    __private unsigned int tmp;

    // 256 x 256 multiply
    multiply256k(a, b, high, low);
    product_low[7] = low[7];
    product_low[6] = low[6];
    product_low[5] = low[5];
    product_low[4] = low[4];
    product_low[3] = low[3];
    product_low[2] = low[2];
    product_low[1] = low[1];
    product_low[0] = low[0];

    // Add 2^32 * high to the low 256 bits (shift left 1 word and add)
    // Affects product[14] to product[6]
    addc(product_low[6], high[7], product_low[6], carry, tmp);
    addc(product_low[5], high[6], product_low[5], carry, tmp);
    addc(product_low[4], high[5], product_low[4], carry, tmp);
    addc(product_low[3], high[4], product_low[3], carry, tmp);
    addc(product_low[2], high[3], product_low[2], carry, tmp);
    addc(product_low[1], high[2], product_low[1], carry, tmp);
    addc(product_low[0], high[1], product_low[0], carry, tmp);

    addc(high[0], 0, product7, carry, tmp);
    product6 = carry;

    carry = 0;

    // Multiply high by 977 and add to low
    // Affects product[15] to product[5]
    for(int i = 7; i >= 0; i--) {
        madd977(&hWord, &t, &high[i], &hWord);
        addc(product_low[i], t, product_low[i], carry, tmp);
        t = 0;
    }
    addc(product7, hWord, high[7], carry, tmp);
    addc(product6, 0, high[6], carry, tmp);

    // Multiply high 2 words by 2^32 and add to low
    // Affects product[14] to product[7]
    carry = 0;

    addc(product_low[6], high[7], product_low[6], carry, tmp);
    addc(product_low[5], high[6], product_low[5], carry, tmp);

    addc(product_low[4], 0, product_low[4], carry, tmp);
    addc(product_low[3], 0, product_low[3], carry, tmp);
    addc(product_low[2], 0, product_low[2], carry, tmp);
    addc(product_low[1], 0, product_low[1], carry, tmp);
    addc(product_low[0], 0, product_low[0], carry, tmp);

    // Multiply top 2 words by 977 and add to low
    // Affects product[15] to product[7]
    carry = 0;
    hWord = 0;
    madd977(&hWord, &t, &high[7], &hWord);
    addc(product_low[7], t, product_low[7], carry, tmp);
    madd977(&hWord, &t, &high[6], &hWord);
    addc(product_low[6], t,  product_low[6], carry, tmp);
    addc(product_low[5], hWord,  product_low[5], carry, tmp);
    // Propagate carry
    addc(product_low[4], 0, product_low[4], carry, tmp);
    addc(product_low[3], 0, product_low[3], carry, tmp);
    addc(product_low[2], 0, product_low[2], carry, tmp);
    addc(product_low[1], 0, product_low[1], carry, tmp);
    addc(product_low[0], 0, product_low[0], carry, tmp);

    // Reduce if >= P
    if(carry || greaterOrEqualToP(product_low)) {
        carry = 0;
        sub256k(product_low, P, product_low, carry, tmp);
    }
}

/**
 * Subtraction mod p
 */
void subModP256k(__private const unsigned int a[8], __private const unsigned int b[8], __private unsigned int c[8])
{
    __private unsigned int borrow = 0;
    __private unsigned int tmp;
    
    sub256k(a, b, c, borrow, tmp);
    
    if (borrow) {
        add256k(c, P, c, borrow, tmp);
    }
}

/**
 * Multiplicative inverse mod P using Fermat's method of x^(p-2) mod p and addition chains
 */
void invModP256k(__private unsigned int x[8])
{
    __private unsigned int y[8] = {0, 0, 0, 0, 0, 0, 0, 1};

    mulModP256k(x, y, y);
    mulModP256k(x, x, x);
    mulModP256k(x, x, x);
    mulModP256k(x, y, y);
    mulModP256k(x, x, x);
    mulModP256k(x, y, y);
    mulModP256k(x, x, x);
    mulModP256k(x, x, x);
    mulModP256k(x, y, y);

    for(int i = 0; i < 5; i++) {
        mulModP256k(x, x, x);
    }

    for(int i = 0; i < 22; i++) {
        mulModP256k(x, y, y);
        mulModP256k(x, x, x);
    }

    mulModP256k(x, x, x);

    for(int i = 0; i < 222; i++) {
        mulModP256k(x, y, y);
        mulModP256k(x, x, x);
    }

    mulModP256k(x, y, x);
}

void addModP256k(__private const unsigned int a[8], __private const unsigned int b[8], __private unsigned int c[8])
{
    __private unsigned int borrow = 0;
    __private unsigned int tmp = 0;

    add256k(a, b, c, borrow, tmp);

    if(borrow) { sub256k(c, P, c, borrow, tmp); }

    else if(c[0] > P[0]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[0] < P[0]) {  }

    else if(c[1] > P[1]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[1] < P[1]) {  }

    else if(c[2] > P[2]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[2] < P[2]) {  }
    
    else if(c[3] > P[3]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[3] < P[3]) {  }
    
    else if(c[4] > P[4]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[4] < P[4]) {  }
    
    else if(c[5] > P[5]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[5] < P[5]) {  }
    
    else if(c[6] > P[6]) { sub256k(c, P, c, borrow, tmp); } 
    else if(c[6] < P[6]) {  }

    else if(c[7] > P[7]) { sub256k(c, P, c, borrow, tmp); } 
}

void doublePoint(__private unsigned int x[8], __private unsigned int y[8]) {
        unsigned int yInv[8];
        unsigned int x3[8];
        unsigned int rx[8];
        unsigned int ry[8];
        unsigned int s[8];

        addModP256k(y, y, yInv);
        invModP256k(yInv);
        
	    // s = 3x^2 / 2y
        mulModP256k(x, x, x3);      // x^2

        addModP256k(x3, x3, s);     // 2x^2
        addModP256k(s, x3, s);      // 3x^2
        mulModP256k(s, yInv, s);    // 3x^2 / 2y

	    //rx = s^2 - 2x
        mulModP256k(s, s, rx);
        subModP256k(rx, x, rx);
        subModP256k(rx, x, rx);

	    //ry = s * (px - rx) - py
        subModP256k(x, rx, ry);
        mulModP256k(s, ry, ry);
        subModP256k(ry, y, ry);

        x[0] = rx[0];
        x[1] = rx[1];
        x[2] = rx[2];
        x[3] = rx[3];
        x[4] = rx[4];
        x[5] = rx[5];
        x[6] = rx[6];
        x[7] = rx[7];

        y[0] = ry[0];
        y[1] = ry[1];
        y[2] = ry[2];
        y[3] = ry[3];
        y[4] = ry[4];
        y[5] = ry[5];
        y[6] = ry[6];
        y[7] = ry[7];
}

void addPoints(__private unsigned int x1[8], __private unsigned int y1[8], __private unsigned int x2[8], __private unsigned int y2[8]) {

	if(equal256k(x1,x2) && equal256k(y1,y2)) {
		doublePoint(x1, y1);
        return;
	}
    
	if(equal256k(x1, x2)) {
        x1[0] = 0xFFFFFFFF;
        x1[1] = 0xFFFFFFFF;
        x1[2] = 0xFFFFFFFF;
        x1[3] = 0xFFFFFFFF;
        x1[4] = 0xFFFFFFFF;
        x1[5] = 0xFFFFFFFF;
        x1[6] = 0xFFFFFFFF;
        x1[7] = 0xFFFFFFFF;
        y1[0] = 0xFFFFFFFF;
        y1[1] = 0xFFFFFFFF;
        y1[2] = 0xFFFFFFFF;
        y1[3] = 0xFFFFFFFF;
        y1[4] = 0xFFFFFFFF;
        y1[5] = 0xFFFFFFFF;
        y1[6] = 0xFFFFFFFF;
        y1[7] = 0xFFFFFFFF;
        return;
	}

	if(isInfinity256k(x1) && isInfinity256k(y1)) {
        x1[0] = x2[0];
        x1[1] = x2[1];
        x1[2] = x2[2];
        x1[3] = x2[3];
        x1[4] = x2[4];
        x1[5] = x2[5];
        x1[6] = x2[6];
        x1[7] = x2[7];

        y1[0] = y2[0];
        y1[1] = y2[1];
        y1[2] = y2[2];
        y1[3] = y2[3];
        y1[4] = y2[4];
        y1[5] = y2[5];
        y1[6] = y2[6];
        y1[7] = y2[7];
        
        return;
	}

	if(isInfinity256k(x2) && isInfinity256k(y2)) {
        return;
	}

    unsigned int rx[8];
    unsigned int ry[8];
    unsigned int s[8];

	subModP256k(y1, y2, rx);

	subModP256k(x1, x2, ry);
    invModP256k(ry);
    mulModP256k(rx, ry, s);

	//rx = (s*s - px - qx) % _p;
    mulModP256k(s, s, rx);
    subModP256k(rx, x1, rx);
    subModP256k(rx, x2, rx);

	//ry = (s * (px - rx) - py) % _p;
    subModP256k(x1, rx, ry);
    mulModP256k(s, ry, ry);
    subModP256k(ry, y1, ry);
    
    x1[0] = rx[0];
    x1[1] = rx[1];
    x1[2] = rx[2];
    x1[3] = rx[3];
    x1[4] = rx[4];
    x1[5] = rx[5];
    x1[6] = rx[6];
    x1[7] = rx[7];

    y1[0] = ry[0];
    y1[1] = ry[1];
    y1[2] = ry[2];
    y1[3] = ry[3];
    y1[4] = ry[4];
    y1[5] = ry[5];
    y1[6] = ry[6];
    y1[7] = ry[7];
}

void generatePubKey(
    __private const unsigned int k[8],
    __private unsigned int x[8],
    __private unsigned int y[8],
    __private volatile const unsigned int gx[64][8],
    __private volatile const unsigned int gy[64][8]
) {

    x[0] = 0xFFFFFFFF;
    x[1] = 0xFFFFFFFF;
    x[2] = 0xFFFFFFFF;
    x[3] = 0xFFFFFFFF;
    x[4] = 0xFFFFFFFF;
    x[5] = 0xFFFFFFFF;
    x[6] = 0xFFFFFFFF;
    x[7] = 0xFFFFFFFF;

    y[0] = 0xFFFFFFFF;
    y[1] = 0xFFFFFFFF;
    y[2] = 0xFFFFFFFF;
    y[3] = 0xFFFFFFFF;
    y[4] = 0xFFFFFFFF;
    y[5] = 0xFFFFFFFF;
    y[6] = 0xFFFFFFFF;
    y[7] = 0xFFFFFFFF;
    
    __private unsigned int px[8];
    __private unsigned int py[8];

	for(int i = 63; i >= 0; i--) {
		if(k[i / 32] >> abs((i - 63) % 32) & 1 ) {
            px[0] = gx[63-i][0];
            px[1] = gx[63-i][1];
            px[2] = gx[63-i][2];
            px[3] = gx[63-i][3];
            px[4] = gx[63-i][4];
            px[5] = gx[63-i][5];
            px[6] = gx[63-i][6];
            px[7] = gx[63-i][7];

            py[0] = gy[63-i][0];
            py[1] = gy[63-i][1];
            py[2] = gy[63-i][2];
            py[3] = gy[63-i][3];
            py[4] = gy[63-i][4];
            py[5] = gy[63-i][5];
            py[6] = gy[63-i][6];
            py[7] = gy[63-i][7];
            addPoints(x, y, px, py);
        }
    }
}

#endif
#ifndef SHA256_CL
#define SHA256_CL

__constant unsigned int _K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

__constant unsigned int _IV[8] = {
    0x6a09e667,
    0xbb67ae85,
    0x3c6ef372,
    0xa54ff53a,
    0x510e527f,
    0x9b05688c,
    0x1f83d9ab,
    0x5be0cd19
};

#define rotr(x, n) ((x) >> (n)) ^ ((x) << (32 - (n)))

#define MAJ(a, b, c) (((a) & (b)) ^ ((a) & (c)) ^ ((b) & (c)))

#define CH(e, f, g) (((e) & (f)) ^ (~(e) & (g)))

#define s0(x) (rotr((x), 7) ^ rotr((x), 18) ^ ((x) >> 3))

#define s1(x) (rotr((x), 17) ^ rotr((x), 19) ^ ((x) >> 10))

#define roundSha(a, b, c, d, e, f, g, h, m, k)\
    t = CH((e), (f), (g)) + (rotr((e), 6) ^ rotr((e), 11) ^ rotr((e), 25)) + (k) + (m);\
    (d) += (t) + (h);\
    (h) += (t) + MAJ((a), (b), (c)) + (rotr((a), 2) ^ rotr((a), 13) ^ rotr((a), 22))

void sha256PublicKey(const unsigned int x[8], const unsigned int y[8], unsigned int digest[8])
{
    __private unsigned int a, b, c, d, e, f, g, h;
    __private unsigned int w[16];
    __private unsigned int t;

    a = _IV[0];
    b = _IV[1];
    c = _IV[2];
    d = _IV[3];
    e = _IV[4];
    f = _IV[5];
    g = _IV[6];
    h = _IV[7];

    w[0] = (x[0] >> 8) | 0x04000000;
    w[1] = (x[1] >> 8) | (x[0] << 24);
    w[2] = (x[2] >> 8) | (x[1] << 24);
    w[3] = (x[3] >> 8) | (x[2] << 24);
    w[4] = (x[4] >> 8) | (x[3] << 24);
    w[5] = (x[5] >> 8) | (x[4] << 24);
    w[6] = (x[6] >> 8) | (x[5] << 24);
    w[7] = (x[7] >> 8) | (x[6] << 24);
    w[8] = (y[0] >> 8) | (x[7] << 24);
    w[9] = (y[1] >> 8) | (y[0] << 24);
    w[10] = (y[2] >> 8) | (y[1] << 24);
    w[11] = (y[3] >> 8) | (y[2] << 24);
    w[12] = (y[4] >> 8) | (y[3] << 24);
    w[13] = (y[5] >> 8) | (y[4] << 24);
    w[14] = (y[6] >> 8) | (y[5] << 24);
    w[15] = (y[7] >> 8) | (y[6] << 24);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[0]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[1]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[2]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[3]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[4]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[5]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[6]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[7]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[8]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[9]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[10]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[11]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[12]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[13]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[14]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[15]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[16]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[17]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[18]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[19]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[20]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[21]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[22]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[23]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[24]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[25]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[26]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[27]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[28]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[29]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[30]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[31]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[32]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[33]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[34]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[35]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[36]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[37]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[38]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[39]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[40]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[41]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[42]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[43]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[44]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[45]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[46]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[47]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[48]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[49]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[50]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[51]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[52]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[53]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[54]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[55]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[56]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[57]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[58]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[59]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[60]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[61]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[62]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[63]);

    a += _IV[0];
    b += _IV[1];
    c += _IV[2];
    d += _IV[3];
    e += _IV[4];
    f += _IV[5];
    g += _IV[6];
    h += _IV[7];

    digest[0] = a;
    digest[1] = b;
    digest[2] = c;
    digest[3] = d;
    digest[4] = e;
    digest[5] = f;
    digest[6] = g;
    digest[7] = h;

    w[0] = (y[7] << 24) | 0x00800000;
    w[15] = 520; // 65 * 8

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[0]);
    roundSha(h, a, b, c, d, e, f, g, 0, _K[1]);
    roundSha(g, h, a, b, c, d, e, f, 0, _K[2]);
    roundSha(f, g, h, a, b, c, d, e, 0, _K[3]);
    roundSha(e, f, g, h, a, b, c, d, 0, _K[4]);
    roundSha(d, e, f, g, h, a, b, c, 0, _K[5]);
    roundSha(c, d, e, f, g, h, a, b, 0, _K[6]);
    roundSha(b, c, d, e, f, g, h, a, 0, _K[7]);
    roundSha(a, b, c, d, e, f, g, h, 0, _K[8]);
    roundSha(h, a, b, c, d, e, f, g, 0, _K[9]);
    roundSha(g, h, a, b, c, d, e, f, 0, _K[10]);
    roundSha(f, g, h, a, b, c, d, e, 0, _K[11]);
    roundSha(e, f, g, h, a, b, c, d, 0, _K[12]);
    roundSha(d, e, f, g, h, a, b, c, 0, _K[13]);
    roundSha(c, d, e, f, g, h, a, b, 0, _K[14]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[15]);

    w[0] = w[0] + s0(0) + 0 + s1(0);
    w[1] = 0 + s0(0) + 0 + s1(w[15]);
    w[2] = 0 + s0(0) + 0 + s1(w[0]);
    w[3] = 0 + s0(0) + 0 + s1(w[1]);
    w[4] = 0 + s0(0) + 0 + s1(w[2]);
    w[5] = 0 + s0(0) + 0 + s1(w[3]);
    w[6] = 0 + s0(0) + w[15] + s1(w[4]);
    w[7] = 0 + s0(0) + w[0] + s1(w[5]);
    w[8] = 0 + s0(0) + w[1] + s1(w[6]);
    w[9] = 0 + s0(0) + w[2] + s1(w[7]);
    w[10] = 0 + s0(0) + w[3] + s1(w[8]);
    w[11] = 0 + s0(0) + w[4] + s1(w[9]);
    w[12] = 0 + s0(0) + w[5] + s1(w[10]);
    w[13] = 0 + s0(0) + w[6] + s1(w[11]);
    w[14] = 0 + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[16]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[17]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[18]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[19]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[20]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[21]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[22]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[23]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[24]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[25]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[26]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[27]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[28]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[29]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[30]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[31]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[32]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[33]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[34]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[35]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[36]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[37]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[38]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[39]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[40]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[41]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[42]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[43]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[44]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[45]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[46]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[47]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[48]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[49]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[50]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[51]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[52]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[53]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[54]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[55]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[56]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[57]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[58]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[59]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[60]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[61]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[62]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[63]);

    digest[0] += a;
    digest[1] += b;
    digest[2] += c;
    digest[3] += d;
    digest[4] += e;
    digest[5] += f;
    digest[6] += g;
    digest[7] += h;
}

void sha256PublicKeyCompressed(const unsigned int x[8], unsigned int yParity, unsigned int digest[8])
{
    __private unsigned int a, b, c, d, e, f, g, h;
    __private unsigned int w[16];
    __private unsigned int t;

    w[0] = 0x02000000 | ((yParity & 1) << 24) | (x[0] >> 8);

    w[1] = (x[1] >> 8) | (x[0] << 24);
    w[2] = (x[2] >> 8) | (x[1] << 24);
    w[3] = (x[3] >> 8) | (x[2] << 24);
    w[4] = (x[4] >> 8) | (x[3] << 24);
    w[5] = (x[5] >> 8) | (x[4] << 24);
    w[6] = (x[6] >> 8) | (x[5] << 24);
    w[7] = (x[7] >> 8) | (x[6] << 24);
    w[8] = (x[7] << 24) | 0x00800000;
    w[15] = 264; // 33 * 8

    a = _IV[0];
    b = _IV[1];
    c = _IV[2];
    d = _IV[3];
    e = _IV[4];
    f = _IV[5];
    g = _IV[6];
    h = _IV[7];

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[0]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[1]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[2]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[3]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[4]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[5]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[6]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[7]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[8]);
    roundSha(h, a, b, c, d, e, f, g, 0, _K[9]);
    roundSha(g, h, a, b, c, d, e, f, 0, _K[10]);
    roundSha(f, g, h, a, b, c, d, e, 0, _K[11]);
    roundSha(e, f, g, h, a, b, c, d, 0, _K[12]);
    roundSha(d, e, f, g, h, a, b, c, 0, _K[13]);
    roundSha(c, d, e, f, g, h, a, b, 0, _K[14]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[15]);

    w[0] = w[0] + s0(w[1]) + 0 + s1(0);
    w[1] = w[1] + s0(w[2]) + 0 + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + 0 + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + 0 + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + 0 + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + 0 + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(0) + w[1] + s1(w[6]);
    w[9] = 0 + s0(0) + w[2] + s1(w[7]);
    w[10] = 0 + s0(0) + w[3] + s1(w[8]);
    w[11] = 0 + s0(0) + w[4] + s1(w[9]);
    w[12] = 0 + s0(0) + w[5] + s1(w[10]);
    w[13] = 0 + s0(0) + w[6] + s1(w[11]);
    w[14] = 0 + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[16]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[17]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[18]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[19]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[20]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[21]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[22]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[23]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[24]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[25]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[26]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[27]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[28]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[29]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[30]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[31]);

    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[32]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[33]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[34]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[35]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[36]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[37]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[38]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[39]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[40]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[41]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[42]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[43]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[44]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[45]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[46]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[47]);


    w[0] = w[0] + s0(w[1]) + w[9] + s1(w[14]);
    w[1] = w[1] + s0(w[2]) + w[10] + s1(w[15]);
    w[2] = w[2] + s0(w[3]) + w[11] + s1(w[0]);
    w[3] = w[3] + s0(w[4]) + w[12] + s1(w[1]);
    w[4] = w[4] + s0(w[5]) + w[13] + s1(w[2]);
    w[5] = w[5] + s0(w[6]) + w[14] + s1(w[3]);
    w[6] = w[6] + s0(w[7]) + w[15] + s1(w[4]);
    w[7] = w[7] + s0(w[8]) + w[0] + s1(w[5]);
    w[8] = w[8] + s0(w[9]) + w[1] + s1(w[6]);
    w[9] = w[9] + s0(w[10]) + w[2] + s1(w[7]);
    w[10] = w[10] + s0(w[11]) + w[3] + s1(w[8]);
    w[11] = w[11] + s0(w[12]) + w[4] + s1(w[9]);
    w[12] = w[12] + s0(w[13]) + w[5] + s1(w[10]);
    w[13] = w[13] + s0(w[14]) + w[6] + s1(w[11]);
    w[14] = w[14] + s0(w[15]) + w[7] + s1(w[12]);
    w[15] = w[15] + s0(w[0]) + w[8] + s1(w[13]);

    roundSha(a, b, c, d, e, f, g, h, w[0], _K[48]);
    roundSha(h, a, b, c, d, e, f, g, w[1], _K[49]);
    roundSha(g, h, a, b, c, d, e, f, w[2], _K[50]);
    roundSha(f, g, h, a, b, c, d, e, w[3], _K[51]);
    roundSha(e, f, g, h, a, b, c, d, w[4], _K[52]);
    roundSha(d, e, f, g, h, a, b, c, w[5], _K[53]);
    roundSha(c, d, e, f, g, h, a, b, w[6], _K[54]);
    roundSha(b, c, d, e, f, g, h, a, w[7], _K[55]);
    roundSha(a, b, c, d, e, f, g, h, w[8], _K[56]);
    roundSha(h, a, b, c, d, e, f, g, w[9], _K[57]);
    roundSha(g, h, a, b, c, d, e, f, w[10], _K[58]);
    roundSha(f, g, h, a, b, c, d, e, w[11], _K[59]);
    roundSha(e, f, g, h, a, b, c, d, w[12], _K[60]);
    roundSha(d, e, f, g, h, a, b, c, w[13], _K[61]);
    roundSha(c, d, e, f, g, h, a, b, w[14], _K[62]);
    roundSha(b, c, d, e, f, g, h, a, w[15], _K[63]);

    digest[0] = a + _IV[0];
    digest[1] = b + _IV[1];
    digest[2] = c + _IV[2];
    digest[3] = d + _IV[3];
    digest[4] = e + _IV[4];
    digest[5] = f + _IV[5];
    digest[6] = g + _IV[6];
    digest[7] = h + _IV[7];
}
#endif
#ifndef BITCOIN_CL
#define BITCOIN_CL

#ifndef endian
#define endian(x) ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((x) >> 24)
#endif

void hashPublicKeyCompressed(const uint256_t x, const unsigned int yParity, unsigned int digest[5])
{
    __private unsigned int hash[8];

    sha256PublicKeyCompressed(x.v, yParity, hash);

    // Swap to little-endian
    hash[0] = endian(hash[0]);
    hash[1] = endian(hash[1]);
    hash[2] = endian(hash[2]);
    hash[3] = endian(hash[3]);
    hash[4] = endian(hash[4]);
    hash[5] = endian(hash[5]);
    hash[6] = endian(hash[6]);
    hash[7] = endian(hash[7]);

    ripemd160sha256NoFinal(hash, digest);
}

void hashPublicKey(const uint256_t x, const uint256_t y, unsigned int digest[5])
{
    __private unsigned int hash[8];

    sha256PublicKey(x.v, y.v, hash);

    // Swap to little-endian
    hash[0] = endian(hash[0]);
    hash[1] = endian(hash[1]);
    hash[2] = endian(hash[2]);
    hash[3] = endian(hash[3]);
    hash[4] = endian(hash[4]);
    hash[5] = endian(hash[5]);
    hash[6] = endian(hash[6]);
    hash[7] = endian(hash[7]);

    ripemd160sha256NoFinal(hash, digest);
}

#endif
#ifndef BLOOMFILTER_CL
#define BLOOMFILTER_CL

bool isInBloomFilter(__private const unsigned int hash[5], __global const unsigned int *targetList, const ulong *mask)
{
    unsigned int h5 = hash[0] + hash[1] + hash[2] + hash[3] + hash[4];

    return (false == 
        (
            (targetList[(((hash[0] << 6) | (h5 & 0x3f)) & *mask) / 32] & (0x01 << ((((hash[0] << 6) | (h5 & 0x3f)) & *mask) % 32))) == 0 ||
            (targetList[(((hash[1] << 6) | ((h5 >> 6) & 0x3f)) & *mask) / 32] & (0x01 << ((((hash[1] << 6) | ((h5 >> 6) & 0x3f)) & *mask) % 32))) == 0 ||
            (targetList[(((hash[2] << 6) | ((h5 >> 12) & 0x3f)) & *mask) / 32] & (0x01 << ((((hash[2] << 6) | ((h5 >> 12) & 0x3f)) & *mask) % 32))) == 0 ||
            (targetList[(((hash[3] << 6) | ((h5 >> 18) & 0x3f)) & *mask) / 32] & (0x01 << ((((hash[3] << 6) | ((h5 >> 18) & 0x3f)) & *mask) % 32))) == 0 || 
            (targetList[ (((hash[4] << 6) | ((h5 >> 24) & 0x3f)) & *mask) / 32] & (0x01 << ( (((hash[4] << 6) | ((h5 >> 24) & 0x3f)) & *mask) % 32))) == 0
        )
    );
}

#endif
#define COMPRESSED 0
#define UNCOMPRESSED 1
#define BOTH 2

typedef struct {
    bool compressed;
    unsigned int privateKey[8];
    unsigned int x[8];
    unsigned int y[8];
    unsigned int digest[5];
} CLDeviceResult;

void addResult(
    __private const bool compressed,
    __private const uint256_t privateKey,
    __private const uint256_t x,
    __private const uint256_t y,
    __private const unsigned int digest[5],
    __global CLDeviceResult* results,
    __global unsigned int* numResults
) {
    CLDeviceResult result;

    result.compressed = compressed;

    result.privateKey[0] = privateKey.v[0];
    result.privateKey[1] = privateKey.v[1];
    result.privateKey[2] = privateKey.v[2];
    result.privateKey[3] = privateKey.v[3];
    result.privateKey[4] = privateKey.v[4];
    result.privateKey[5] = privateKey.v[5];
    result.privateKey[6] = privateKey.v[6];
    result.privateKey[7] = privateKey.v[7];

    result.x[0] = x.v[0];
    result.x[1] = x.v[1];
    result.x[2] = x.v[2];
    result.x[3] = x.v[3];
    result.x[4] = x.v[4];
    result.x[5] = x.v[5];
    result.x[6] = x.v[6];
    result.x[7] = x.v[7];

    result.y[0] = y.v[0];
    result.y[1] = y.v[1];
    result.y[2] = y.v[2];
    result.y[3] = y.v[3];
    result.y[4] = y.v[4];
    result.y[5] = y.v[5];
    result.y[6] = y.v[6];
    result.y[7] = y.v[7];

    ripemd160FinalRound(digest, result.digest);

    results[atomic_add(numResults, 1)] = result;
}

/*
void div2(unsigned int a[8], unsigned int c[8]) {
    unsigned int carry = 0;

    c[0] = (a[0] + carry) >> 1;
    carry = a[0] % 2;
    
    c[1] = (a[1] + carry) >> 1;
    carry = a[1] % 2;
    
    c[2] = (a[2] + carry) >> 1;
    carry = a[2] % 2;
    
    c[3] = (a[3] + carry) >> 1;
    carry = a[3] % 2;
    
    c[4] = (a[4] + carry) >> 1;
    carry = a[4] % 2;
    
    c[5] = (a[5] + carry) >> 1;
    carry = a[5] % 2;

    c[6] = (a[6] + carry) >> 1;
    carry = a[6] % 2;

    c[7] = (a[7] + carry) >> 1;
}
*/

__kernel void _stepKernel(
    const unsigned int i,
    const unsigned int totalPoints,
    __global uint256_t* incXPtr,
    __global uint256_t* incYPtr,
    __global const unsigned int* targetList,
    const ulong mask,
    __global CLDeviceResult* results,
    __global unsigned int* numResults) {

    __private uint256_t k = {{0,0,0,0,0,0,0,1}};
    k.v[7] = get_local_id(0);

__private unsigned int gx[64][8] = { 
    { 0x79be667e, 0xf9dcbbac, 0x55a06295, 0xce870b07, 0x29bfcdb, 0x2dce28d9, 0x59f2815b, 0x16f81798 },
    { 0xc6047f94, 0x41ed7d6d, 0x3045406e, 0x95c07cd8, 0x5c778e4b, 0x8cef3ca7, 0xabac09b9, 0x5c709ee5 },
    { 0xe493dbf1, 0xc10d80f3, 0x581e4904, 0x930b1404, 0xcc6c1390, 0xee07584, 0x74fa94ab, 0xe8c4cd13 },
    { 0x2f01e5e1, 0x5cca351d, 0xaff3843f, 0xb70f3c2f, 0xa1bdd05, 0xe5af888a, 0x67784ef3, 0xe10a2a01 },
    { 0xe60fce93, 0xb59e9ec5, 0x3011aabc, 0x21c23e97, 0xb2a31369, 0xb87a5ae9, 0xc44ee89e, 0x2a6dec0a },
    { 0xd30199d7, 0x4fb5a22d, 0x47b6e054, 0xe2f378ce, 0xdacffcb8, 0x9904a61d, 0x75d0dbd4, 0x7143e65 },
    { 0xbf23c154, 0x2d16eab7, 0xb1051ea, 0xf832823c, 0xfc4c6f1d, 0xcdbafd81, 0xe37918e6, 0xf874ef8b },
    { 0x34ff3be4, 0x33f7a06, 0x696c3d09, 0xf7d1671c, 0xbcf55cd7, 0x535655, 0x64707745, 0x6769a24e },
    { 0x82822632, 0x12c609d9, 0xea2a6e3e, 0x172de238, 0xd8c39cab, 0xd5ac1ca1, 0x646e23f, 0xd5f51508 },
    { 0x465370b2, 0x87a79ff3, 0x905a857a, 0x9cf918d5, 0xadbc968, 0xd9e159d0, 0x926e2c00, 0xef34a24d },
    { 0x241febb8, 0xe23cbd77, 0xd664a18f, 0x66ad6240, 0xaaec6ecd, 0xc813b088, 0xd5b901b2, 0xe285131f },
    { 0x5d1bdb4e, 0xa172fa79, 0xfce4cc29, 0x83d8f8d9, 0xfc318b85, 0xf423de0d, 0xedcb6306, 0x9b920471 },
    { 0x175e159f, 0x728b865a, 0x72f99cc6, 0xc6fc846d, 0xe0b93833, 0xfd2222ed, 0x73fce5b5, 0x51e5b739 },
    { 0x423a013f, 0x3ff32d7, 0xa5ffbcc8, 0xe139c621, 0x30fdfeb5, 0xc6da121b, 0xce78049e, 0x46bc47d6 },
    { 0x111d6a45, 0xac1fb905, 0x8907a7a, 0xbcd68776, 0x49df662f, 0x3b3e2741, 0x302df6f7, 0x8416824a },
    { 0x4a4a6dc9, 0x7ac7c8b8, 0xad795dbe, 0xbcb9dcff, 0x7290b68a, 0x5ef74e56, 0xab5edde0, 0x1bced775 },
    { 0x363d90d4, 0x47b00c9c, 0x99ceac05, 0xb6262ee0, 0x53441c7e, 0x55552ffe, 0x526bad8f, 0x83ff4640 },
    { 0x4c1b9866, 0xed9a7e9b, 0x553973c6, 0xc93b02bf, 0xb62fb01, 0x2edfb59d, 0xd2712a5c, 0xaf92c541 },
    { 0xa4083877, 0xba83b12b, 0x529a2f3c, 0x780b54e, 0x3233edbc, 0x1a28f135, 0xe0c8f28c, 0xbeaaf3d1 },
    { 0xa804c641, 0xd28cc0b5, 0x3a4e3e1a, 0x2f56c86f, 0x6e0d880a, 0x454203b9, 0x8cd3db5a, 0x7940d33a },
    { 0x8b4b5f16, 0x5df3c2be, 0x8c6244b5, 0xb7456388, 0x43e4a781, 0xa15bcd1b, 0x69f79a55, 0xdffdf80c },
    { 0xed0c5ce4, 0xe1329171, 0x8ce17c7e, 0xc83c6110, 0x71af64ee, 0x417c997a, 0xbb3f2671, 0x4755e4be },
    { 0xfaecb013, 0xc44ce694, 0xb3b15c3f, 0x83f1fae8, 0xe5325456, 0x6e0552ce, 0xd4b6e6c8, 0x7cec8ab },
    { 0x9bb8a13, 0x2dcad2f2, 0xc8731a0b, 0x37cbcafd, 0xb3b2dd82, 0x4f23cd3e, 0x7f64eae, 0x9ad1b1f7 },
    { 0x723cbaa6, 0xe5db996d, 0x6bf771c0, 0xbd548c7, 0xb700dbff, 0xa6c0e77b, 0xcb611592, 0x5232fcda },
    { 0x57efa786, 0x437b744d, 0x343d7dc4, 0x5773a3c6, 0x2d240a43, 0x7984907, 0x1fd383d6, 0xca030d5 },
    { 0x264bbd43, 0x6a28bc42, 0xa2df7e9c, 0xd5226cb9, 0x1080577e, 0x327b012a, 0x7fafc777, 0xc584dd5 },
    { 0xa94c6524, 0xbd40d2bb, 0xdac85c05, 0x6236a79d, 0xa78bc61f, 0xd5bdec9d, 0x2bf26bd8, 0x4b2438e8 },
    { 0xeebfa4d4, 0x93bebf98, 0xba5feec8, 0x12c2d3b5, 0x9479612, 0x37a91983, 0x9a533eca, 0xe7dd7fa },
    { 0x381c4ad7, 0xa7a97bfd, 0xa61c6031, 0xc118495f, 0xc4ea4bc0, 0x8f6766d6, 0x76bee908, 0x47d297fd },
    { 0xe1efb9cd, 0x5adc63b, 0xcce10831, 0xd9538c47, 0x9cf1d05f, 0xefdd08b2, 0x448d7042, 0x2ede454c },
    { 0x5318f9b1, 0xa2697010, 0xc5ac235e, 0x9af475a8, 0xc7e5419f, 0x33d47b18, 0xd33feeb3, 0x29eb99a4 },
    { 0x100f44da, 0x696e7167, 0x2791d0a0, 0x9b7bde45, 0x9f1215a2, 0x9b3c03bf, 0xefd7835b, 0x39a48db0 },
    { 0x8c0989f2, 0xceb5c771, 0xa8415dff, 0x2b4c4199, 0xd8d9c8f9, 0x237d0808, 0x4b05284f, 0x1e4df706 },
    { 0xfb8f153c, 0x5e266704, 0xc4a48174, 0x3262c025, 0x9c528539, 0xbc95bc1b, 0xb1e63c33, 0xdc47bffd },
    { 0xe747333f, 0xd75d5175, 0x5a0cc9f0, 0xa7287084, 0x65a02c58, 0x7737a7b8, 0xb8fa1b8b, 0x4bb2629a },
    { 0xe1031be2, 0x62c7ed1b, 0x1dc9227a, 0x4a04c017, 0xa77f8d44, 0x64f3b385, 0x2c8acde6, 0xe534fd2d },
    { 0xf4b93f22, 0x4c8089ea, 0xb9f95dcd, 0xf29b2c9, 0x28a6ac5, 0xde94d857, 0x84e27e36, 0xa95c8356 },
    { 0x9d1aca1, 0xfce55236, 0xb19622ea, 0x25b08b0, 0xd51e8512, 0xf97e696c, 0x20d62fe1, 0x7b160e8a },
    { 0xc66c59cc, 0x454c2b9e, 0x18a2ad79, 0x3821cde7, 0x518b3a93, 0xbfc39562, 0xe97d7d04, 0x75ba7fc2 },
    { 0xfeea6cae, 0x46d55b53, 0xac2839f, 0x143bd7ec, 0x5cf8b266, 0xa41d6af5, 0x2d5e688d, 0x9094696d },
    { 0x4d000b62, 0x1adb87e1, 0xc53261af, 0x9db2e179, 0x141ecae0, 0xb331a187, 0xaa4040a, 0xee752b08 },
    { 0x71f570ca, 0x203da05d, 0xd6aa2621, 0x14717128, 0xd657a040, 0x3e1f1b77, 0xf89962fd, 0x475c58ef },
    { 0xa2b7b362, 0x9f7bd253, 0xb7d282b5, 0xc21da014, 0x46b4821d, 0xc65e7651, 0x6048b060, 0x43ff8359 },
    { 0xda67a91d, 0x91049cdc, 0xb367be4b, 0xe6ffca3c, 0xfeed657d, 0x808583de, 0x33fa978b, 0xc1ec6cb1 },
    { 0x4dbacd36, 0x5fa1ef58, 0x7c0c0cfa, 0xaf00d871, 0x8bbd9f35, 0xccea5a83, 0x5ee3cc82, 0x1fe741c9 },
    { 0x13d1ffc4, 0x81509bee, 0xe68f17d8, 0xff41c259, 0xf4c85f1, 0x52686050, 0x87eda8ba, 0xb4e218da },
    { 0x219b4f9c, 0xef6c6000, 0x7659c79c, 0x45b0533b, 0x3cc9d916, 0xce29dbff, 0x133b40ca, 0xa2e96db8 },
    { 0x53904faa, 0xb334cdd, 0xa6e00093, 0x5ef22151, 0xec08d0f7, 0xbb11069f, 0x57545ccc, 0x1a37b7c0 },
    { 0x1a575af, 0x9d414675, 0x3cf99119, 0x6316995d, 0x2a6ee7aa, 0xad0f85ad, 0x57cd0f1f, 0x38a47ca9 },
    { 0xf5f0e043, 0x7621d439, 0xca71f5c1, 0xb76155d6, 0xd3a61a83, 0xd3c20c6e, 0xe309d755, 0xe315565b },
    { 0x8f506f0b, 0x6c0b6e9a, 0x57a7f36d, 0x970ca4e3, 0x47cbc921, 0x46227642, 0xcbe781d9, 0xf5362d33 },
    { 0x8e7bcd0b, 0xd35983a7, 0x719cca77, 0x64ca9067, 0x79b53a04, 0x3a9b8bca, 0xeff959f4, 0x3ad86047 },
    { 0x33b35baa, 0x195e729d, 0xc350f319, 0x996950df, 0x3bc15b8d, 0x3d0389e7, 0x77d2808b, 0xf13f0351 },
    { 0x374deeae, 0x22c93f95, 0x5cb83ad2, 0x71f7e22, 0x56f6e109, 0xcad7bca6, 0xd71dc7b2, 0x4414bb36 },
    { 0x2380c09c, 0x7f3aeae5, 0x7c46e073, 0x95aeb0dc, 0x944dbaf2, 0xb62a9f0c, 0x5e8a64ad, 0x6ae7d616 },
    { 0x385eed34, 0xc1cdff21, 0xe6d08186, 0x89b81bde, 0x71a7f4f1, 0x8397e669, 0xa841e15, 0x99c43862 },
    { 0xf6f62208, 0x3daf5480, 0x456be13, 0x4d5f67d1, 0x47c82642, 0xbefc1ce2, 0xdc83a270, 0x78f2827c },
    { 0xfb26e518, 0x8f953de2, 0xbd70cb3c, 0x3d1fc255, 0xcd91c3ce, 0x7d8c6f36, 0x9d893209, 0x715adcb6 },
    { 0x89912259, 0x11b9132d, 0x28f5c6bc, 0x763ceab7, 0xd18c3706, 0xe8bd1d7, 0xed44db75, 0x60788c1e },
    { 0x6f9d9b8, 0x3ecf191, 0x637c73a4, 0x413dfa18, 0xfddf84a, 0x5947fbc9, 0xc606ed86, 0xc3fac3a7 },
    { 0xae86eeea, 0x252b411c, 0x1cdc36c2, 0x84482939, 0xda1745e5, 0xa7e4da17, 0x5c9d2274, 0x4b7fd72d },
    { 0x2248c9f9, 0xbbfff55, 0xe61d2f8c, 0x56dc2c48, 0x8718be75, 0xcf36f2ee, 0x7a147426, 0x7c169290 },
    { 0xe11a6e16, 0xe05c4407, 0x4ac11b48, 0xd94085d0, 0xa99f0877, 0xdd1c6f76, 0xfd0dac4b, 0xb50964e3 }
};
__private unsigned int gy[64][8] = { 
    { 0x483ada77, 0x26a3c465, 0x5da4fbfc, 0xe1108a8, 0xfd17b448, 0xa6855419, 0x9c47d08f, 0xfb10d4b8 },
    { 0x1ae168fe, 0xa63dc339, 0xa3c58419, 0x466ceaee, 0xf7f63265, 0x3266d0e1, 0x236431a9, 0x50cfe52a },
    { 0x51ed993e, 0xa0d455b7, 0x5642e209, 0x8ea51448, 0xd967ae33, 0xbfbdfe40, 0xcfe97bdc, 0x47739922 },
    { 0x5c4da8a7, 0x41539949, 0x293d082a, 0x132d13b4, 0xc2e213d6, 0xba5b7617, 0xb5da2cb7, 0x6cbde904 },
    { 0xf7e35073, 0x99e59592, 0x9db99f34, 0xf5793710, 0x1296891e, 0x44d23f0b, 0xe1f32cce, 0x69616821 },
    { 0x95038d9d, 0xae3d5c3, 0xb3d6dec9, 0xe9838065, 0x1f760cc3, 0x64ed8196, 0x5b3ff1f, 0x24106ab9 },
    { 0x5cb3866f, 0xc3300373, 0x7ad928a0, 0xba5392e4, 0xc522fc54, 0x811e2f78, 0x4dc37efe, 0x66831d9f },
    { 0x5d9d1162, 0x3a236c55, 0x3f6619d8, 0x9832098c, 0x55df16c3, 0xe8f8b681, 0x8491067a, 0x73cc2f1a },
    { 0x11f8a809, 0x8557dfe4, 0x5e8256e8, 0x30b60ace, 0x62d613ac, 0x2f7b17be, 0xd31b6eaf, 0xf6e26caf },
    { 0x35e531b3, 0x8368c082, 0xa4af8bda, 0xfdeec2c1, 0x588e09b2, 0x15d37a10, 0xa2f8fb20, 0xb33887f4 },
    { 0x513378d9, 0xff94f8d3, 0xd6c420bd, 0x13981df8, 0xcd50fd0f, 0xbd0cb5af, 0xabb3e66f, 0x2750026d },
    { 0x28438267, 0x79379e2e, 0x794bb994, 0x38a22656, 0x79eb1e99, 0x96c56e7b, 0x70330666, 0xf7b83103 },
    { 0xd3506e0d, 0x9e3c79eb, 0xa4ef97a5, 0x1ff71f5e, 0xacb5955a, 0xdd24345c, 0x6efa6ffe, 0xe9fed695 },
    { 0xb91ae00f, 0xe1e1d970, 0xa1179f7b, 0xbaf6b3c7, 0x720d8ec3, 0x524f009e, 0xd1236e6d, 0x8b548a34 },
    { 0x696911c, 0x478eaffb, 0xb90d48db, 0xff065952, 0xf0700089, 0x96daca4c, 0xa9a111d4, 0x2108e9d0 },
    { 0x529911b0, 0x16631e72, 0x943ef9f7, 0x39c0f457, 0x1de90cdb, 0x424742ac, 0xb2bf8f68, 0xa78dd66d },
    { 0x4e273ad, 0xfc732221, 0x953b4453, 0x97f33631, 0x45b9a890, 0x8199ecb, 0x62003c7f, 0x3bee9de9 },
    { 0xc1f792d3, 0x20be8a0f, 0x7fbcb753, 0xce56e69c, 0xc652ead7, 0xe43eb1ad, 0x72c4f3fd, 0xc68fe020 },
    { 0x40e9f612, 0xfeefbc79, 0xb8bf83d6, 0x9361b3e2, 0x2001e757, 0x6ed1ef90, 0xb12b534d, 0xf0b254b9 },
    { 0x95be8325, 0x2b2fa6d0, 0x3dec2842, 0xc16047e8, 0x1af18ca8, 0x9cf736a9, 0x43ce95fa, 0x6d46967a },
    { 0x4aad0a6f, 0x68d308b4, 0xb3fbd781, 0x3ab0da04, 0xf9e33654, 0x6162ee56, 0xb3eff0c6, 0x5fd4fd36 },
    { 0x221a9fc7, 0xbc2345bd, 0xbf3dad7f, 0x5a7ea680, 0x49d93925, 0x763ddab1, 0x63f9fa6e, 0xa07bf42f },
    { 0xcc09b5e9, 0xe9ecb57, 0xfc2e02c6, 0xec2fb13d, 0x9c32b286, 0xb85e2e2e, 0x8981dfd9, 0xab155070 },
    { 0x945bb2b2, 0xafeee3b9, 0xb6f9dd28, 0x4f863e85, 0xf54a840, 0xf4752d53, 0x64130627, 0xc3811c80 },
    { 0x96e867b5, 0x595cc498, 0xa9211374, 0x88824d6e, 0x2660a065, 0x37794948, 0x1dc069d, 0x9eb39f5f },
    { 0xd712db0b, 0xd1b48518, 0x893627c9, 0x28de03ec, 0x689b6d2a, 0xe5e9974a, 0xb07ab442, 0x74b02f9e },
    { 0xd87c6fa9, 0x4ee093b4, 0xd4f75ce2, 0x4c33be22, 0x6a118243, 0x717b8d8d, 0xe6122793, 0x7704ab11 },
    { 0xb5201fd9, 0x92f96280, 0xfd792195, 0x5019e3a, 0x7e5d3c60, 0xa0e39b2b, 0xc2e2c8db, 0xf18661f4 },
    { 0x5d9a8ca3, 0x970ef0f2, 0x69ee7eda, 0xf178089d, 0x9ae4cdc3, 0xa711f712, 0xddfd4fda, 0xe1de8999 },
    { 0x936af53b, 0x238eeee4, 0x8f3e5fa7, 0x9915ecc, 0xf0451032, 0xdb939c00, 0x93ace318, 0x7d493fc5 },
    { 0xecb4530, 0xd8af9be7, 0xb0154c1f, 0xfe477123, 0x464e3244, 0xa7a2d4c6, 0xad9fd233, 0xa8913797 },
    { 0xf44ccfeb, 0x4beda419, 0x5772d93a, 0xebb405e8, 0xa41f2b40, 0xd1e3ec65, 0x2c726eee, 0xfe91f92d },
    { 0xcdd9e131, 0x92a00b77, 0x2ec8f330, 0xc090666, 0xb7ff4a18, 0xff5195ac, 0xfbd5cd6, 0x2bc65a09 },
    { 0xfb4dbd04, 0x4f432034, 0xffd2172c, 0xb9dc966c, 0x60de6bf5, 0x156511aa, 0x736ac5a3, 0x5d72fa98 },
    { 0x6ca27a9d, 0xc5e06218, 0x16fa11d9, 0xb4bccd53, 0x1dde1389, 0xac542613, 0x90a45dd, 0xd949b095 },
    { 0xf2affe01, 0x45070c11, 0x4cc43603, 0x804c2581, 0xc88376aa, 0x6e1a969a, 0x9f8d961a, 0x6946f6d6 },
    { 0x9d706192, 0x8940405e, 0x6bb6a417, 0x6597535a, 0xf292dd41, 0x9e1ced79, 0xa44f18f2, 0x9456a00d },
    { 0xa67a92ec, 0x62962df, 0xb0e5f6a7, 0xa40eee90, 0xc37ef134, 0x4915609a, 0xbd5861b9, 0xbe001fd3 },
    { 0x1153188f, 0x5101f0c6, 0x3e56692c, 0xe0d8c27e, 0x6fe9e0ee, 0x9212b5e5, 0x34e050c5, 0x7ca04c44 },
    { 0xd9592fe2, 0xbfb30fcf, 0xbea4f3ce, 0xaac10cb2, 0xf00a60dd, 0xb1595597, 0x7ec3c69c, 0xf75f5956 },
    { 0xe57c6b6c, 0x97dce1ba, 0xb06e4e12, 0xbf3ecd5c, 0x981c8957, 0xcc41442d, 0x3155debf, 0x18090088 },
    { 0x6a0d5b8f, 0x18e0d255, 0xcb6d8255, 0x82d972cc, 0xcb7df5f1, 0x19c7293a, 0x3e72851f, 0x48302cea },
    { 0xeb42415b, 0x95dc880d, 0xd2555734, 0x5bc95b8d, 0xf2445d00, 0xc3363e7d, 0xf8649a72, 0xd35d420e },
    { 0x69303894, 0x1695122d, 0x57a937a3, 0xf71e29c9, 0x10d10835, 0x46f3835, 0xa2397fec, 0xfe86fec2 },
    { 0x9bacaa35, 0x481642bc, 0x41f463f7, 0xec9780e5, 0xdec7adc5, 0x8f740a1, 0x7e9ea8e2, 0x7a68be1d },
    { 0x16c3540e, 0x8a51892e, 0x7fdcfd59, 0xe838299d, 0xcc384a0, 0x9fc0535f, 0x60be10f8, 0x338eb623 },
    { 0x6008391f, 0xa991961d, 0xcecb9337, 0xb1b758bd, 0xa4ad0120, 0x6d5bd127, 0xe0db419d, 0xdb191c19 },
    { 0x24d9c605, 0xd959efea, 0xf5a44180, 0xc0372a6e, 0x394f8ac5, 0x3e905765, 0x27df01a7, 0x8d3b6bc7 },
    { 0x5bc087d0, 0xbc80106d, 0x88c9ecca, 0xc20d3c1c, 0x13999981, 0xe1443469, 0x9dcb096b, 0x22771c8 },
    { 0x3038f1cb, 0x8ab20dc3, 0xcc55fc52, 0xe1bb8698, 0xbdb93c5d, 0x9f4d7ea6, 0x67c5df2e, 0x77ebcdb7 },
    { 0x6b9f4e62, 0xbe5a052b, 0xf6218916, 0xdf7101a, 0xa5bf61bf, 0x3ed7e40a, 0x678430af, 0xdd2ecc82 },
    { 0x469f955d, 0x2afa6171, 0x9530c542, 0x4f1c3368, 0x48cf925d, 0x43bb8eaf, 0x30487d0c, 0x87fa243f },
    { 0x10b7770b, 0x2a3da4b3, 0x94031042, 0xca95145, 0x79e88e2e, 0x47fd68b3, 0xea10047e, 0x8460372a },
    { 0xa58a0185, 0x640abf87, 0xf9464036, 0x248d52bc, 0xaa6560ef, 0xbc889b70, 0x2bc503cc, 0xcb8d7418 },
    { 0x171165b6, 0x4fcd4f99, 0x16032c06, 0xf806f729, 0x3828d663, 0xe54321, 0x7875bea9, 0x8daf734a },
    { 0x6f8e8619, 0x3464956a, 0xf1598aef, 0xd509b09a, 0x93af9214, 0x8f846756, 0x99be48, 0x161bbc1a },
    { 0x283bebc3, 0xe8ea23f5, 0x6701de19, 0xe9ebf457, 0x6b304eec, 0x2086dc8c, 0xc0458fe5, 0x542e5453 },
    { 0x1bcd4e81, 0x7de73a0f, 0xaf2c5715, 0xb367cee7, 0xe657ca74, 0x48321bf6, 0xd15b20b5, 0x20aaa102 },
    { 0xf3e12881, 0x1012a34d, 0x58e846a7, 0x19d01769, 0x16d2cb31, 0xb8b7ab54, 0x49dbca3b, 0x58ba68f3 },
    { 0xda8b4d98, 0x7cc9ac9b, 0x27b87635, 0x59b136fa, 0x36969c84, 0xfdef9e11, 0x635c4222, 0x8e8f0ef1 },
    { 0x7c80c68e, 0x603059ba, 0x69b8e2a3, 0xe45c4d4, 0x7ea4dd2f, 0x5c281002, 0xd8689060, 0x3a842160 },
    { 0x19e993c9, 0x707302f9, 0x62ab0ace, 0x589ff0e9, 0x8d921155, 0x1472f728, 0x2334cb7a, 0x4eee38bc },
    { 0xfa059469, 0x2d21eed7, 0xa506bb55, 0xb435ba18, 0xe1637502, 0x35da2be2, 0x369d8a12, 0x883ea257 },
    { 0x87d6065b, 0x87a2d430, 0xe1ad5e25, 0x96f0af24, 0x17adc6e1, 0x38318c6f, 0x767fbf8b, 0x682bfc8 },
};
    /** 
    k.v[7] = i * (get_local_id(0) + 1) + 1;

    __private unsigned int mult[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
    mult[7] = totalPoints;
    multiply256k(k.v, mult, mult, k.v);
    */
    __private unsigned int digest[5];
    __private uint256_t x;
    __private uint256_t y;

    generatePubKey(k.v, x.v,y.v, gx, gy);
    
    #if defined(COMPRESSION_UNCOMPRESSED) || defined(COMPRESSION_BOTH)
        hashPublicKey(x, y, digest);
        if(isInBloomFilter(digest, targetList, &mask)) {
            addResult(false, k, x, y, digest, results, numResults);
        }
    #endif
    #if defined(COMPRESSION_COMPRESSED) || defined(COMPRESSION_BOTH)
        hashPublicKeyCompressed(x, y.v[7], digest);
        if(isInBloomFilter(digest, targetList, &mask)) {
            addResult(true, k, x, y, digest, results, numResults);
        }
    #endif

    /**
    for(int i = 0; i < totalPoints; i++) {
        addPoints(x.v, y.v, incXPtr, incYPtr);

        #if defined(COMPRESSION_UNCOMPRESSED) || defined(COMPRESSION_BOTH)
            hashPublicKey(x, y, digest);
            if(isInBloomFilter(digest, targetList, &mask)) {
                addResult(false, k, x, y, digest, results, numResults);
            }
        #endif
        #if defined(COMPRESSION_COMPRESSED) || defined(COMPRESSION_BOTH)
            hashPublicKeyCompressed(x, y.v[7], digest);
            if(isInBloomFilter(digest, targetList, &mask)) {
                addResult(true, k, x, y, digest, results, numResults);
            }
        #endif
    }
    */
}
