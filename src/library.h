#ifndef _LIBRARY_H
#define _LIBRARY_H

#include "common.h"
#include "io.h"
#include "lsh.h"
#include "table.h"
#include "taxonomy.h"
#include "filesystem.hpp"

class Library
{
public:
  Library(const char *library_dirpath,
          const char *taxonomy_dirpath,
          const char *input_filepath,
          uint8_t k,
          uint8_t w,
          uint8_t h,
          uint8_t b,
          RankingMethod ranking_method,
          bool adaptive_size,
          uint64_t capacitiy_size,
          uint32_t num_batch_rows,
          bool from_library = false,
          bool from_kmers = false,
          uint16_t target_batch = 0,
          bool only_init = false,
          bool update_annotations = false,
          bool fast_mode = false,
          bool verbose = true,
          bool log = true);

private:
  TaxonomyInput _taxonomy_input;
  TaxonomyRecord<tT> _taxonomy_record;
  const char *_library_dirpath;
  const char *_taxonomy_dirpath;
  const char *_input_filepath;
  uint8_t _k;
  uint8_t _w;
  uint8_t _h;
  uint8_t _b;
  RankingMethod _ranking_method;
  bool _adaptive_size;
  uint64_t _num_rows; // Total number of rows: pow(2, 2h).
  uint64_t _capacity_size;
  uint32_t _tbatch_size;
  uint16_t _total_batches;
  uint64_t _root_size;
  maskLSH _lsh_vg;
  std::vector<uint8_t> _positions;
  std::vector<uint8_t> _npositions;
  std::unordered_map<tT, inputStream<encT>> _inputStream_map;
  std::unordered_map<tT, uint64_t> _basis_to_size;
  std::unordered_map<tT, size_t> _basis_to_ninput;
  std::unordered_map<tT, uint64_t> _tID_to_size;
  std::unordered_map<tT, float> _tID_to_length;
  std::unordered_map<tT, uint32_t> _tID_to_ngenomes;
  std::vector<tT> _tID_vec;
  uint32_t _num_species;
  uint32_t _num_nodes;
  bool _from_library;
  bool _input_kmers;
  uint16_t _target_batch;
  bool _only_init;
  bool _update_annotations;
  bool _fast_mode;
  bool _log;
  bool _verbose;
  const tT _rootID = 1;

public:
  uint64_t getConstrainedSizeKC(tT curr_tID);
  uint64_t getConstrainedSizeSC(uint64_t num_basis);
  void getBatchHTd(HTd<encT> *td, unsigned int curr_batch);
  void getBatchHTs(HTs<encT> *ts, unsigned int curr_batch);
  bool saveBatchHTs(HTs<encT> &ts, unsigned int curr_batch);
  bool loadBatchHTs(HTs<encT> &ts, unsigned int curr_batch);
  bool saveMetadata();
  bool loadMetadata();
  void getRandomPositions();
  void skipBatch();
  void buildTables();
  void annotateInfo();
  void processLeaf(tT tID_key);
  void resetInfo(HTs<encT> &ts, bool reset_scount, bool reset_tlca);
  void softLCA(HTs<encT> &ts, unsigned int curr_batch);
  void countBasis(HTs<encT> &ts, unsigned int curr_batch);
  decltype(_npositions) &npositions() { return _npositions; }
  decltype(_positions) &positions() { return _positions; }
  decltype(_lsh_vg) &lsh_vg() { return _lsh_vg; }
};

#endif
