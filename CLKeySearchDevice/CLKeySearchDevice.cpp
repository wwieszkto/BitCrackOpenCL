#include <cmath>
#include "Logger.h"
#include "util.h"
#include "CLKeySearchDevice.h"

// Defined in bitcrack_cl.cpp which gets build in the pre-build event
extern char _bitcrack_cl[];

typedef struct {
    bool compressed;
    unsigned int privateKey[8];
    unsigned int x[8];
    unsigned int y[8];
    unsigned int digest[5];
} CLDeviceResult;

static void undoRMD160FinalRound(const unsigned int hIn[5], unsigned int hOut[5])
{
    unsigned int iv[5] = {
        0x67452301,
        0xefcdab89,
        0x98badcfe,
        0x10325476,
        0xc3d2e1f0
    };

    for(int i = 0; i < 5; i++) {
        hOut[i] = util::endian(hIn[i]) - iv[(i + 1) % 5];
    }
}

CLKeySearchDevice::CLKeySearchDevice(uint64_t device, int threads, int pointsPerThread, int blocks, int compressionMode)
{
    _threads = threads;
    _blocks = blocks;
    _points = pointsPerThread * threads * blocks;
    _device = (cl_device_id)device;

    if(threads <= 0 || threads % 32 != 0) {
        throw KeySearchException("KEYSEARCH_THREAD_MULTIPLE_EXCEPTION", "The number of threads must be a multiple of 32");
    }

    if(pointsPerThread <= 0) {
        throw KeySearchException("KEYSEARCH_MINIMUM_POINT_EXCEPTION", "At least 1 point per thread required");
    }

    std::string options = "";

    switch (compressionMode) {
        case PointCompressionType::COMPRESSED:
            options += " -DCOMPRESSION_COMPRESSED";
        break;
        case PointCompressionType::UNCOMPRESSED:
            options += " -DCOMPRESSION_UNCOMPRESSED";
        break;
        case PointCompressionType::BOTH:
            options += " -DCOMPRESSION_BOTH";
        break;
    }
    try {
        // Create the context
        _clContext = new cl::CLContext(_device);
        Logger::log(LogLevel::Info, "Compiling OpenCL kernels...");
        _clProgram = new cl::CLProgram(*_clContext, _bitcrack_cl, options);

        // Load the kernel
        _stepKernel = new cl::CLKernel(*_clProgram, "_stepKernel");

        _globalMemSize = _clContext->getGlobalMemorySize();

        _deviceName = _clContext->getDeviceName();
    } catch(cl::CLException ex) {
        throw KeySearchException(ex.msg, ex.description);
    }

    _iterations = 0;
}

CLKeySearchDevice::~CLKeySearchDevice()
{
    _clContext->free(_xInc);
    _clContext->free(_yInc);
    _clContext->free(_deviceResults);
    _clContext->free(_deviceResultsCount);

    delete _stepKernel;
    delete _initKeysKernel;
    delete _clContext;
}

uint64_t CLKeySearchDevice::getOptimalBloomFilterMask(double p, size_t n)
{
    double m = 3.6 * ceil((n * std::log(p)) / std::log(1 / std::pow(2, std::log(2))));

    unsigned int bits = (unsigned int)std::ceil(std::log(m) / std::log(2));

    return ((uint64_t)1 << bits) - 1;
}

void CLKeySearchDevice::initializeBloomFilter(const std::vector<struct hash160> &targets, uint64_t mask)
{
    size_t sizeInWords = (mask + 1) / 32;
    _targetMemSize = sizeInWords * sizeof(uint32_t);

    Logger::log(LogLevel::Info, "Initializing BloomFilter (" + util::format("%.1f", (double)_targetMemSize / (double)(1024)) + "KB)");

    uint32_t *buf = new uint32_t[sizeInWords];

    for(size_t i = 0; i < sizeInWords; i++) {
        buf[i] = 0;
    }

    for(unsigned int k = 0; k < targets.size(); k++) {

        unsigned int hash[5];
        unsigned int h5 = 0;

        uint64_t idx[5];

        undoRMD160FinalRound(targets[k].h, hash);

        for(int i = 0; i < 5; i++) {
            h5 += hash[i];
        }

        idx[0] = ((hash[0] << 6) | (h5 & 0x3f)) & mask;
        idx[1] = ((hash[1] << 6) | ((h5 >> 6) & 0x3f)) & mask;
        idx[2] = ((hash[2] << 6) | ((h5 >> 12) & 0x3f)) & mask;
        idx[3] = ((hash[3] << 6) | ((h5 >> 18) & 0x3f)) & mask;
        idx[4] = ((hash[4] << 6) | ((h5 >> 24) & 0x3f)) & mask;

        for(int i = 0; i < 5; i++) {
            uint64_t j = idx[i];
            buf[j / 32] |= 1 << (j % 32);
        }
    }

    _deviceTargetList.mask = mask;
    _deviceTargetList.ptr = _clContext->malloc(_targetMemSize);
    _deviceTargetList.size = targets.size();
    _clContext->copyHostToDevice(buf, _deviceTargetList.ptr, _targetMemSize);

    delete[] buf;
}

void CLKeySearchDevice::allocateBuffers()
{
    size_t numKeys = (size_t)_points;
    size_t size = numKeys * 8 * sizeof(unsigned int);

    _bufferMemSize = 
        8 * sizeof(unsigned int) +       // _xInc
        8 * sizeof(unsigned int) +       // _yInc
        128 * sizeof(CLDeviceResult) +   // _deviceResults
        sizeof(unsigned int);            // _deviceResultsCount

    Logger::log(LogLevel::Info, "Allocating Memory for Buffers (" + util::format("%.1f", (double)_bufferMemSize / (double)(1024 * 1024)) + "MB)");


    // Value to increment by
    _xInc = _clContext->malloc(8 * sizeof(unsigned int), CL_MEM_READ_ONLY);
    _yInc = _clContext->malloc(8 * sizeof(unsigned int), CL_MEM_READ_ONLY);

    // Buffer for storing results
    _deviceResults = _clContext->malloc(128 * sizeof(CLDeviceResult));
    _deviceResultsCount = _clContext->malloc(sizeof(unsigned int));
}

void CLKeySearchDevice::setIncrementor(secp256k1::ecpoint &p)
{
    unsigned int buf[8];

    p.x.exportWords(buf, 8, secp256k1::uint256::BigEndian);
    _clContext->copyHostToDevice(buf, _xInc, 8 * sizeof(unsigned int));

    p.y.exportWords(buf, 8, secp256k1::uint256::BigEndian);
    _clContext->copyHostToDevice(buf, _yInc, 8 * sizeof(unsigned int));
}

void CLKeySearchDevice::init(const secp256k1::uint256 &start, int compression, const secp256k1::uint256 &stride)
{
    if(start.cmp(secp256k1::N) >= 0) {
        throw KeySearchException("KEYSEARCH_STARTINGKEY_OUT_OF_RANGE", "Starting key is out of range");
    }

    _start = start;

    _stride = stride;

    _compression = compression;

    try {
        allocateBuffers();

        // Set the incrementor
        secp256k1::ecpoint p = secp256k1::multiplyPoint(secp256k1::uint256((uint64_t)_points ) * _stride, secp256k1::G());
        setIncrementor(p);
    } catch(cl::CLException ex) {
        throw KeySearchException(ex.msg, ex.description);
    }
}

void CLKeySearchDevice::doStep()
{
    try {
        _stepKernel->set_args(
            _iterations,
            _points,
            _xInc,
            _yInc,
            _deviceTargetList.ptr,
            _deviceTargetList.mask,
            _deviceResults,
            _deviceResultsCount);
        _stepKernel->call(_blocks, _threads);
        fflush(stdout);

        getResultsInternal();

        _iterations++;
    } catch(cl::CLException ex) {
        throw KeySearchException(ex.msg, ex.description);
    }
}

void CLKeySearchDevice::setBloomFilter()
{
    uint64_t bloomFilterMask = getOptimalBloomFilterMask(1.0e-9, _targetList.size());

    initializeBloomFilter(_targetList, bloomFilterMask);
}

void CLKeySearchDevice::setTargetsInternal()
{
    // Clean up existing list
    if(_deviceTargetList.ptr != NULL) {
        _clContext->free(_deviceTargetList.ptr);
    }

    setBloomFilter();
}

void CLKeySearchDevice::setTargets(const std::set<KeySearchTarget> &targets)
{
    try {
        _targetList.clear();

        for(std::set<KeySearchTarget>::iterator i = targets.begin(); i != targets.end(); ++i) {
            hash160 h(i->value);
            _targetList.push_back(h);
        }

        setTargetsInternal();
    } catch(cl::CLException ex) {
        throw KeySearchException(ex.msg, ex.description);
    }
}

size_t CLKeySearchDevice::getResults(std::vector<KeySearchResult> &results)
{
    size_t count = _results.size();
    for(size_t i = 0; i < count; i++) {
        results.push_back(_results[i]);
    }
    _results.clear();

    return count;
}

uint64_t CLKeySearchDevice::keysPerStep()
{
    return (uint64_t)_points;
}

std::string CLKeySearchDevice::getDeviceName()
{
    return _deviceName;
}

void CLKeySearchDevice::getMemoryInfo(uint64_t &freeMem, uint64_t &totalMem)
{
    freeMem = _globalMemSize - _targetMemSize - _pointsMemSize - _bufferMemSize;
    totalMem = _globalMemSize;
}

void CLKeySearchDevice::splatBigInt(unsigned int *ptr, int idx, secp256k1::uint256 &k)
{
    unsigned int buf[8];

    k.exportWords(buf, 8, secp256k1::uint256::BigEndian);

    memcpy(ptr + idx * 8, buf, sizeof(unsigned int) * 8);

}

bool CLKeySearchDevice::isTargetInList(const unsigned int hash[5])
{
    size_t count = _targetList.size();

    while(count) {
        if(memcmp(hash, _targetList[count - 1].h, 20) == 0) {
            return true;
        }
        count--;
    }

    return false;
}

void CLKeySearchDevice::removeTargetFromList(const unsigned int hash[5])
{
    size_t count = _targetList.size();

    while(count) {
        if(memcmp(hash, _targetList[count - 1].h, 20) == 0) {
            _targetList.erase(_targetList.begin() + count - 1);
            return;
        }
        count--;
    }
}

void CLKeySearchDevice::getResultsInternal()
{
    unsigned int numResults = 0;

    _clContext->copyDeviceToHost(_deviceResultsCount, &numResults, sizeof(unsigned int));

    if(numResults > 0) {
        CLDeviceResult *ptr = new CLDeviceResult[numResults];

        _clContext->copyDeviceToHost(_deviceResults, ptr, sizeof(CLDeviceResult) * numResults);

        unsigned int actualCount = 0;

        for(unsigned int i = 0; i < numResults; i++) {

            // might be false-positive
            if(!isTargetInList(ptr[i].digest)) {
                continue;
            }
            actualCount++;

            KeySearchResult minerResult;

            minerResult.privateKey = secp256k1::uint256(ptr[i].privateKey, secp256k1::uint256::BigEndian);
            minerResult.compressed = ptr[i].compressed;

            memcpy(minerResult.hash, ptr[i].digest, 20);

            minerResult.publicKey = secp256k1::ecpoint(secp256k1::uint256(ptr[i].x, secp256k1::uint256::BigEndian), secp256k1::uint256(ptr[i].y, secp256k1::uint256::BigEndian));

            removeTargetFromList(ptr[i].digest);

            _results.push_back(minerResult);
        }
        
        delete[] ptr;

        // Reset device counter
        numResults = 0;
        _clContext->copyHostToDevice(&numResults, _deviceResultsCount, sizeof(unsigned int));
    }
}

secp256k1::uint256 CLKeySearchDevice::readBigInt(unsigned int *src, int idx)
{
    unsigned int value[8] = {0};

    for(int k = 0; k < 8; k++) {
        value[k] = src[idx * 8 + k];
    }

    secp256k1::uint256 v(value, secp256k1::uint256::BigEndian);

    return v;
}

secp256k1::uint256 CLKeySearchDevice::getNextKey()
{
    return _start + secp256k1::uint256((uint64_t)_points) * _iterations * _stride;
}