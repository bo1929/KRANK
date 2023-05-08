#include "index.h"

uint32_t
fillIndexTable(char* fpath, std::vector<std::vector<uint32_t>>& index_table, uint8_t k, uint8_t h, unsigned int batch_size, unsigned int num_threads)
{
  vectorLSH mask_lsh = generateMaskLSH(k, h);

  kseq_t* reader = getReader(fpath);
  batch_size = adjustBatchSize(batch_size, num_threads);
  std::vector<sseq_t> seqBatch = readBatch(reader, batch_size);
  uint32_t total_ix = 0;

  while (!(seqBatch.empty())) {
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (uint32_t ix = 0; ix < (uint32_t)seqBatch.size(); ++ix) {
      const char* kmer_seq;
      uint64_t enc_bp;
      uint64_t row_val;

      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingBPCompute(kmer_seq, enc_bp);
      row_val = computeValueLSH(enc_bp, mask_lsh);
      assert(row_val < index_table.size());
#pragma omp critical
      {
        index_table[row_val].push_back(total_ix + ix);
      }
    }
    total_ix += (uint32_t)seqBatch.size();
    if (seqBatch.size() == batch_size)
      seqBatch = readBatch(reader, batch_size);
    else
      seqBatch.clear();
  }
  return total_ix;
}

std::vector<uint64_t>
computeHistNumCols(std::vector<std::vector<uint32_t>> index_table)
{
  std::vector<uint64_t> vecHist;
  for (auto i = 0; i < index_table.size(); ++i) {
    uint64_t c = index_table[i].size();
    if (vecHist.size() < c)
      vecHist.resize(c + 1);
    ++vecHist[c];
  }
  return vecHist;
}
