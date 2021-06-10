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

    __private unsigned int digest[5];
    __private uint256_t x;
    __private uint256_t y;

    generatePubKey(k.v, x.v,y.v);
    
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
        addPoints(x, y, incXPtr, incYPtr);

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
