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
  Library(char* library_dirpath,
          char* nodes_filepath,
          char* input_filepath,
          uint8_t k,
          uint8_t h,
          uint8_t b,
          uint64_t capacitiy_size,
          uint32_t num_batch_rows,
          bool onDisk = true);
  uint64_t getConstrainedSize(std::set<tT> tIDsBasis);
  void getBatchHTd(HTd<encT>& td);
  void getBatchHTs(HTs<encT>& ts, uint8_t curr_depth = 0, uint8_t last_depth = 1);

private:
  TaxonomyNCBI _taxonomy_ncbi;
  TaxonomyRecord<tT> _taxonomy_record;
  char* _library_dirpath;
  char* _nodes_filepath;
  char* _input_filepath;
  uint8_t _k;
  uint8_t _h;
  uint8_t _b;
  uint64_t _capacity_size;
  uint32_t _tbatch_size;
  maskLSH _lsh_vg;
  std::unordered_map<tT, StreamIM<encT>> _streamIM_map;
  std::unordered_map<tT, StreamOD<encT>> _streamOD_map;
  std::unordered_map<tT, uint64_t> _basis_to_size;
  bool _onDisk;
  uint64_t _root_size;
  uint64_t _num_rows;
  const uint16_t _rootID = 1;
};

#endif
