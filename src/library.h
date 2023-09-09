#ifndef _LIBRARY_H
#define _LIBRARY_H

#include "common.h"
#include "io.h"
#include "lsh.h"
#include "table.h"
#include "taxonomy.h"

class Library
{
public:
  Library(const char* library_dirpath,
          const char* nodes_filepath,
          const char* input_filepath,
          uint8_t k,
          uint8_t h,
          uint8_t b,
          uint64_t capacitiy_size,
          uint32_t num_batch_rows,
          bool switch_ranking = true,
          bool in_library = false,
          bool on_disk = true,
          uint8_t specified_batch = 0,
          bool log = true);

private:
  TaxonomyNCBI _taxonomy_ncbi;
  TaxonomyRecord<tT> _taxonomy_record;
  const char* _library_dirpath;
  const char* _nodes_filepath;
  const char* _input_filepath;
  uint8_t _k;
  uint8_t _h;
  uint8_t _b;
  uint64_t _num_rows; // Total number of rows: pow(2, 2h).
  uint64_t _capacity_size;
  uint32_t _tbatch_size;
  uint16_t _total_batches;
  uint64_t _root_size;
  maskLSH _lsh_vg;
  std::vector<uint8_t> _positions;
  std::vector<uint8_t> _npositions;
  std::unordered_map<tT, StreamIM<encT>> _streamIM_map;
  std::unordered_map<tT, StreamOD<encT>> _streamOD_map;
  std::unordered_map<tT, uint64_t> _basis_to_size;
  std::vector<tT> _tID_vec;
  uint64_t _num_species;
  bool _switch_ranking;
  bool _on_disk;
  bool _in_library;
  uint8_t _specified_batch;
  bool _log;
  const uint16_t _rootID = 1;

public:
  uint64_t getConstrainedSize(std::set<tT> tIDsBasis);
  void getBatchHTd(HTd<encT>* td);
  void getBatchHTs(HTs<encT>* ts, uint8_t curr_depth = 0, uint8_t last_depth = 1);
  bool saveBatchHTs(HTs<encT>& ts, uint16_t curr_batch);
  bool saveMetadata();
  bool loadMetadata();
  void getRandomPositions();
  void skipBatch();
  void run(uint8_t sdepth = 3, RankingMethod init_ranking = information_score);
  decltype(_npositions)& npositions() { return _npositions; }
  decltype(_positions)& positions() { return _positions; }
  decltype(_lsh_vg)& lsh_vg() { return _lsh_vg; }
};

#endif
