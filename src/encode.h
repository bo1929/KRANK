#ifndef _ENCODE_H
#define _ENCODE_H

#include "common.h"
#include "io.h"

void
kmerEncodingBPCompute(const char* seq, uint64_t& enc_bp);

void
kmerEncodingCompute(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingComputeC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingUpdate(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingUpdateC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
retrieveEncodings(char* fpath, uint64_t*& enc_arr, uint32_t num_kmers, unsigned int batch_size, unsigned int num_threads);

uint8_t
computeHammingDistanceEnc64(uint64_t x, uint64_t y); // TODO: CHECK

#endif
