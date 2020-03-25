// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <cassert>
#include <stdexcept>

#include <FairMQLogger.h>

#include "PreClusterSerializer.h"

using namespace std;

namespace o2
{
namespace mch
{

//_________________________________________________________________________________________________
PreClusterSerializer::PreClusterSerializer(void* buffer, uint32_t size, bool write)
{
  /// constructor

  // initialization is the same as in the "reset" method
  reset(buffer, size, write);
}

//_________________________________________________________________________________________________
int PreClusterSerializer::reset(void* buffer, uint32_t size, bool write)
{
  /// initialize the internal structure in write mode

  assert(buffer != nullptr);

  int status = 0;

  // store buffer
  mBuffer = buffer;

  // store buffer size
  mSize = size;

  // store write mode
  mWriteMode = write;

  // reset
  const uint32_t minSizeOfCluster = SSizeOfUShort + SSizeOfDigit;
  mSize4PreCluster = (size > minSizeOfCluster) ? size - minSizeOfCluster : 0;
  mSize4Digit = (size > SSizeOfDigit) ? size - SSizeOfDigit : 0;
  mLastNDigits = nullptr;
  //mPreClusters.clear();

  if (size >= SSizeOfUShort) {

    // store number of precluster and increment buffer
    mNPreClusters = reinterpret_cast<uint16_t*>(mBuffer);
    mBuffer = (reinterpret_cast<uint16_t*>(mBuffer) + 1);
    mCurrentSize = SSizeOfUShort;

    printf("[PreClusterSerializer::reset] mNPreClusters=%p\n", mNPreClusters);

    if (mWriteMode) {
      // assign 0 clusters in write mode
      *mNPreClusters = 0;
    } else {
      // read buffer otherwise
      status = 0; //readBuffer();
    }

  } else {

    mCurrentSize = mSize + 1;
    mNPreClusters = nullptr;

    if (mWriteMode) {
      LOG(ERROR) << "The buffer is too small to store the data block.";
      status = -ENOBUFS;
    } else {
      LOG(ERROR) << "The buffer is too small to contain the data block.";
      status = -EILSEQ;
    }
  }

  return status;
}

//_________________________________________________________________________________________________
int PreClusterSerializer::startPreCluster(const Digit& digit)
{
  /// start a new precluster

  assert(mWriteMode);

  // could move back to constructor (or reset)
  if (mCurrentSize > mSize4PreCluster) {
    LOG(ERROR) << "The buffer is too small to store a new precluster.";
    mLastNDigits = nullptr;
    return -ENOBUFS;
  }

  // use current buffer position to store number of digits
  // and increment buffer
  mLastNDigits = reinterpret_cast<uint16_t*>(mBuffer);
  *mLastNDigits = 1;
  mBuffer = (reinterpret_cast<uint16_t*>(mBuffer) + 1);
  mCurrentSize += SSizeOfUShort;

  // store digit and increment buffer
  auto lastDigit = reinterpret_cast<Digit*>(mBuffer);
  *lastDigit = digit;
  printf("[PreClusterSerializer::startPreCluster] *mLastNDigits=%d -> %p\n", (int)(*mLastNDigits), mLastNDigits);
  mBuffer = (reinterpret_cast<Digit*>(mBuffer) + 1);
  mCurrentSize += SSizeOfDigit;

  // increment number of pre-clusters
  (*mNPreClusters) += 1;
  printf("[PreClusterSerializer::startPreCluster] *mNPreClusters=%d -> %p\n", (int)(*mNPreClusters), mNPreClusters);

  // insert in list
  //PreClusterStruct preCluster = {*mLastNDigits, lastDigit};
  //mPreClusters.push_back(preCluster);

  return 0;
}

//_________________________________________________________________________________________________
int PreClusterSerializer::addDigit(const Digit& digit)
{
  /// add a new digit to the current precluster

  assert(mWriteMode);
  if (mCurrentSize > mSize4Digit) {
    LOG(ERROR) << "The buffer is too small to store a new digit.";
    return -ENOBUFS;
  }

  if (!mLastNDigits) {
    LOG(ERROR) << "No precluster to attach the new digit to.";
    return -EPERM;
  }

  // increment number of digits
  *mLastNDigits += 1;
  printf("[PreClusterSerializer::addDigit] *mLastNDigits=%d -> %p\n", (int)(*mLastNDigits), mLastNDigits);

  // assign digit to the buffer and increment buffer
  *(reinterpret_cast<Digit*>(mBuffer)) = digit;
  mBuffer = (reinterpret_cast<Digit*>(mBuffer) + 1);
  mCurrentSize += SSizeOfDigit;

  // increment number of digits in the stored cluster
  //if (!mPreClusters.empty()) {
  //  ++mPreClusters.back().nDigits;
  //}

  return 0;
}

//_________________________________________________________________________________________________
int PreClusterSerializer::readBuffer(std::vector<PreClusterStruct>& preClusters)
{
  /// read the buffer to fill the internal structure
  /// fNPreclus,fNdig,dig1, dig2 ... dign, fNdig, dig1 ... digN, ... last_precluster

  assert(!mWriteMode);

  // make sure mNPreClusters was assigned
  if (mNPreClusters == nullptr) {
    return -EILSEQ;
  }

  // initialize status
  int status = 0;

  // rolling pre-cluster
  PreClusterStruct preCluster = {};

  // loop over
  printf("*mNPreClusters: %d\n", (int)(*mNPreClusters));
  for (int i = 0; i < *mNPreClusters; ++i) {

    // store number of digits from buffer and increment
    uint16_t* nDigits = nullptr;
    if (mCurrentSize < mSize) {
      nDigits = reinterpret_cast<uint16_t*>(mBuffer);
      printf("[PreClusterSerializer::readBuffer] pre-cluster %d  *nDigits=%d\n", i, (int)(*nDigits));
      mBuffer = (reinterpret_cast<uint16_t*>(mBuffer) + 1);
      mCurrentSize += SSizeOfUShort;
    } else {
      LOG(ERROR) << "Cannot read the expected number of preclusters.";
      status = -EILSEQ;
      break;
    }

    // read the digits
    if (nDigits && (*nDigits) > 0 && mCurrentSize + (*nDigits) * SSizeOfDigit <= mSize) {

      auto digit = reinterpret_cast<Digit*>(mBuffer);
      mBuffer = (reinterpret_cast<Digit*>(mBuffer) + (*nDigits));
      mCurrentSize += (*nDigits) * SSizeOfDigit;

      // store
      preCluster.nDigits = *nDigits;
      preCluster.digits = digit;
      preClusters.push_back(preCluster);

    } else {

      if (*nDigits == 0) {
        LOG(ERROR) << "The precluster cannot contain 0 digit.";
      } else {
        LOG(ERROR) << "Cannot read the expected number of digits.";
      }

      status = -EILSEQ;
      break;
    }
  }

  // sanity check on read size
  if (mCurrentSize != mSize && status >= 0) {
    LOG(ERROR) << "The number of bytes read differs from the buffer size: "<<mCurrentSize<<", "<<mSize;
    status = -EILSEQ;
  }

  // reset in case of negative status
  if (status < 0) {
    mCurrentSize = mSize + 1;
    mNPreClusters = nullptr;
    preClusters.clear();
  }

  return status;
}

//_________________________________________________________________
std::ostream& operator<<(std::ostream& stream, const PreClusterStruct& cluster)
{
  stream << "{nDigits= " << cluster.nDigits;
  for (int i = 0; i < cluster.nDigits; ++i) {
    stream << ", digit[" << i << "]= " << cluster.digits[i];
  }
  stream << "}";

  return stream;
}

//_________________________________________________________________
std::ostream& operator<<(std::ostream& stream, const PreClusterSerializer& clusterBlock)
{
  stream << "{fNClusters= " << clusterBlock.getNPreClusters() << std::endl;
  //const auto& clusters(clusterBlock.getPreClusters());
  //int i(0);
  //for (const auto& cluster : clusters) {
  //  stream << "  cluster[" << i++ << "]= " << cluster << std::endl;
  //}
  stream << "}";

  return stream;
}

//_________________________________________________________________________________________________
uint32_t PreClusterSerializer::getPreClustersBufferSize(PreClusterFinder& finder)
{
  const uint32_t SSizeOfInt = sizeof(int);
  const uint32_t SSizeOfUShort = sizeof(unsigned short int);
  const uint32_t SSizeOfDigit = sizeof(Digit);

  // number of DEs with preclusters and total number of pads used
  int nUsedDigits(0);
  int nDEWithPreClusters = finder.getNDEWithPreClusters(nUsedDigits);
  uint32_t nPreClusters = finder.getNumberOfPreClusters();

  uint32_t bufSize = SSizeOfInt + nDEWithPreClusters * 2 * SSizeOfInt +
      (nDEWithPreClusters + nPreClusters) * SSizeOfUShort + nUsedDigits * SSizeOfDigit;

  return bufSize;
}

//_________________________________________________________________________________________________
void PreClusterSerializer::storePreClusters(PreClusterFinder& preClusterFinder, char* buffer)
{
  const uint32_t SSizeOfInt = sizeof(int);
  const uint32_t SSizeOfUShort = sizeof(unsigned short int);
  const uint32_t SSizeOfDigit = sizeof(Digit);

  // get the total size needed for the buffer
  auto sizeOut = getPreClustersBufferSize(preClusterFinder);

  int nUsedDigits(0);
  int nDEWithPreClusters = preClusterFinder.getNDEWithPreClusters(nUsedDigits);

  // store the number of DE with preclusters
  printf("[PreClusterSerializer::storePreClusters] nDEWithPreClusters=%d -> %p\n", nDEWithPreClusters, buffer);
  memcpy(buffer, &nDEWithPreClusters, SSizeOfInt);
  buffer += SSizeOfInt;
  sizeOut -= SSizeOfInt;

  /// store the preclusters in the given buffer
  const PreClusterFinder::PreCluster* cluster(nullptr);
  const o2::mch::Digit* digit(nullptr);
  uint32_t* bytesUsed(nullptr);
  uint32_t totalBytesUsed(0);

  for (int iDE = 0, nDEs = preClusterFinder.getNDEs(); iDE < nDEs; ++iDE) {

    if (!preClusterFinder.hasPreClusters(iDE)) {
      continue;
    }

    // store the DE ID
    if (sizeOut - totalBytesUsed >= SSizeOfInt) {
      auto deId(reinterpret_cast<int*>(buffer + totalBytesUsed));
      *deId = preClusterFinder.getDEId(iDE);
      totalBytesUsed += SSizeOfInt;
    } else {
      throw length_error("cannot store DE ID");
    }

    // prepare to store the size of the PreClusterSerializer
    if (sizeOut - totalBytesUsed >= SSizeOfInt) {
      bytesUsed = reinterpret_cast<uint32_t*>(buffer + totalBytesUsed);
      totalBytesUsed += SSizeOfInt;
    } else {
      throw length_error("cannot store size of the PreClusterSerializer");
    }

    // prepare to store the preclusters of this DE
    printf("[PreClusterSerializer::storePreClusters] totalBytesUsed: %d\n", (int)(totalBytesUsed));
    if (reset(buffer + totalBytesUsed, sizeOut - totalBytesUsed, true) < 0) {
      throw length_error("cannot reset the cluster block");
    }

    for (int iPlane = 0; iPlane < 2; ++iPlane) {
      for (int iCluster = 0, nClusters = preClusterFinder.getNPreClusters(iDE, iPlane);
           iCluster < nClusters; ++iCluster) {

        cluster = preClusterFinder.getPreCluster(iDE, iPlane, iCluster);
        if (!cluster->storeMe) {
          continue;
        }

        // add the precluster with its first digit
        printf("[PreClusterSerializer::storePreClusters] adding pre-cluster of size %d\n", (int)(cluster->lastPad+1-cluster->firstPad));
        // Fix code to use O2 digits
        digit = preClusterFinder.getDigit(iDE, cluster->firstPad);
        if (startPreCluster(*digit) < 0) {
          throw length_error("cannot store a new precluster");
        }

        // loop over other pads and add corresponding digits
        for (uint16_t iOrderedPad = cluster->firstPad + 1; iOrderedPad <= cluster->lastPad; ++iOrderedPad) {
          digit = preClusterFinder.getDigit(iDE, iOrderedPad);
          if (addDigit(*digit) < 0) {
            throw length_error("cannot store a new digit");
          }
        }
      }
    }

    // store the size of the PreClusterSerializer
    *bytesUsed = getCurrentSize();
    printf("[PreClusterSerializer::storePreClusters] bytesUsed: %d -> %p\n", (int)(*bytesUsed), bytesUsed);
    totalBytesUsed += *bytesUsed;

    //if (mPrint) {
    //  LOG(INFO) << "block: " << mPreClusterSerializer;
    //}
  }

  if (totalBytesUsed != sizeOut) {
    throw length_error("incorrect payload");
  }
}

//_________________________________________________________________________________________________
void PreClusterSerializer::readPreClusters(std::vector<PreClusterStruct>& preClusters, char* buffer, uint32_t size)
{
  const uint32_t SSizeOfInt = sizeof(int);
  const uint32_t SSizeOfUShort = sizeof(unsigned short int);
  const uint32_t SSizeOfDigit = sizeof(Digit);

  uint32_t sizeIn = size;
  uint32_t totalBytesUsed(0);

  int nDEWithPreClusters;
  // get the number of DE with preclusters
  memcpy(&nDEWithPreClusters, buffer, SSizeOfInt);
  printf("[PreClusterSerializer::readPreClusters] nDEWithPreClusters: %d -> %p\n", (int)(nDEWithPreClusters), buffer);
  buffer += SSizeOfInt;
  sizeIn -= SSizeOfInt;
  //totalBytesUsed += SSizeOfInt;

  /// store the preclusters in the given buffer
  const PreClusterFinder::PreCluster* cluster(nullptr);
  const o2::mch::Digit* digit(nullptr);
  int deId;
  uint32_t bytesUsed;


  for (int iDE = 0; iDE < nDEWithPreClusters; ++iDE) {

    // retrieve the DE ID
    if (sizeIn - totalBytesUsed >= SSizeOfInt) {
      auto iptr(reinterpret_cast<int*>(buffer + totalBytesUsed));
      deId = *iptr;
      totalBytesUsed += SSizeOfInt;
    } else {
      throw length_error("cannot retrieve DE ID");
    }

    // retrieve the size of the PreClusterSerializer
    if (sizeIn - totalBytesUsed >= SSizeOfInt) {
      auto iptr = reinterpret_cast<uint32_t*>(buffer + totalBytesUsed);
      bytesUsed = *iptr;
      printf("[PreClusterSerializer::readPreClusters] bytesUsed: %d -> %p, totalBytesUsed: %d\n", (int)(bytesUsed), iptr, (int)(totalBytesUsed));
     totalBytesUsed += SSizeOfInt;
    } else {
      throw length_error("cannot store size of the PreClusterSerializer");
    }

    // prepare to read the preclusters of this DE
    printf("[PreClusterSerializer::readPreClusters] totalBytesUsed: %d\n", (int)(totalBytesUsed));
    if (reset(buffer + totalBytesUsed, bytesUsed/*sizeIn - totalBytesUsed*/, false) < 0) {
      throw length_error("cannot reset the cluster block");
    }

    if (readBuffer(preClusters) < 0) {
      throw length_error("cannot read the cluster block");
    }

    // increment the total size
    totalBytesUsed += bytesUsed;

    //if (mPrint) {
    //  LOG(INFO) << "block: " << mPreClusterSerializer;
    //}
  }

  if (totalBytesUsed != sizeIn) {
    printf("[PreClusterSerializer::readPreClusters] %d != %d\n", (int)totalBytesUsed, (int)sizeIn);
    throw length_error("incorrect payload");
  }
}

} // namespace mch
} // namespace o2
