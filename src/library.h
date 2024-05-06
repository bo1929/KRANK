#ifndef _LIBRARY_H
#define _LIBRARY_H

#include "common.h"
#include "io.h"
#include "lsh.h"
#include "table.h"
#include "taxonomy.h"
#include "filesystem.hpp"

enum LabelsLCA
{
  hard_lca,
  soft_lca,
};

class Library
{
public:
  Library(const char *library_dirpath,
          const char *tax_dirpath,
          const char *input_filepath,
          uint8_t k,
          uint8_t w,
          uint8_t h,
          uint8_t b,
          RankingMethod ranking_method,
          LabelsLCA labels_lca,
          bool adaptive_size,
          uint64_t capacitiy_size,
          uint32_t num_batch_rows,
          bool from_library = false,
          bool from_kmers = false,
          uint16_t target_batch = 0,
          bool only_init = false,
          bool update_annotations = false,
          bool fast_mode = false,
          bool remove_intermediate = true,
          bool verbose = true,
          bool log = true);

private:
  TaxonomyInput _tax_input;
  TaxonomyRecord<tT> _tax_record;
  const char *_library_dirpath;
  const char *_tax_dirpath;
  const char *_input_filepath;
  uint8_t _k;
  uint8_t _w;
  uint8_t _h;
  uint8_t _b;
  RankingMethod _ranking_method;
  LabelsLCA _labels_lca;
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
  std::unordered_map<tT, uint64_t> _trID_to_size;
  std::unordered_map<tT, float> _trID_to_length;
  std::unordered_map<tT, uint32_t> _trID_to_ngenomes;
  std::vector<tT> _trID_vec;
  uint32_t _num_species;
  uint32_t _num_nodes;
  bool _from_library;
  bool _input_kmers;
  uint16_t _target_batch;
  bool _only_init;
  bool _update_annotations;
  bool _fast_mode;
  bool _remove_intermediate;
  bool _log;
  bool _verbose;
  const tT _rootrID = 1;

public:
  uint64_t getConstrainedSizeKC(tT curr_trID);
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
  void processLeaf(tT trID_key);
  void resetInfo(HTs<encT> &ts, bool reset_scount, bool reset_tlca);
  void labelLCAs(HTs<encT> &ts, unsigned int curr_batch);
  void countBasis(HTs<encT> &ts, unsigned int curr_batch);
  decltype(_npositions) &npositions() { return _npositions; }
  decltype(_positions) &positions() { return _positions; }
  decltype(_lsh_vg) &lsh_vg() { return _lsh_vg; }
};

#endif
