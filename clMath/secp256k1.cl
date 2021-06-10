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

    unsigned int rise[8];
	subModP256k(y1, y2, rise);

    unsigned int run[8];
    unsigned int s[8];
	subModP256k(x1, x2, run);
    invModP256k(run);
    mulModP256k(rise, run, s);

	//rx = (s*s - px - qx) % _p;
    unsigned int rx[8];
    mulModP256k(s, s, rx);
    subModP256k(rx, x1, rx);
    subModP256k(rx, x2, rx);

	//ry = (s * (px - rx) - py) % _p;
    unsigned int ry[8];
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

void generatePubKey(__private const unsigned int k[8],__private unsigned int x[8],__private unsigned int y[8]) {

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
    
    unsigned int gx[8] = { 0x79BE667E, 0xF9DCBBAC, 0x55A06295, 0xCE870B07, 0x029BFCDB, 0x2DCE28D9, 0x59F2815B, 0x16F81798 };
    unsigned int gy[8] = { 0x483ADA77, 0x26A3C465, 0x5DA4FBFC, 0x0E1108A8, 0xFD17B448, 0xA6855419, 0x9C47D08F, 0xFB10D4B8 };

	for(int i = 255; i >= 0; i--) {
		if(k[i / 32] >> abs((i - 255) % 32) & 1 ) {
            addPoints(x, y, gx, gy);
        }
        doublePoint(gx, gy);
    }
}

#endif
