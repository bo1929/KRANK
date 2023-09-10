#ifndef _QUERY_H
#define _QUERY_H

#include "common.h"
#include "encode.h"
#include "library.h"
#include "lsh.h"
#include "taxonomy.h"

class Query
{
public:
  Query(std::vector<std::string> library_dirpaths,
        const char* output_dirpath,
        const char* query_filepath,
        uint8_t max_match_hdist = 5,
        bool save_match_info = true,
        bool log = false);
  void run(uint64_t rbatch_size = DEFAULT_BATCH_SIZE);

private:
  struct SearchL
  {
    SearchL(const char* library_dirpath, bool log = false);
    bool _log;
    const char* _library_dirpath;
    uint8_t _k;
    uint8_t _h;
    uint8_t _b;
    uint64_t _num_rows; // Total number of rows: pow(2, 2h).
    uint64_t _capacity_size;
    uint32_t _tbatch_size;
    uint16_t _total_batches;
    uint16_t _num_species;
    uint64_t _root_size;
    maskLSH _lsh_vg;
    std::vector<uint8_t> _positions;
    std::vector<uint8_t> _npositions;
    std::unordered_map<tT, uint64_t> _tID_to_taxID;
    tT _tax_num_input;
    tT _tax_num_nodes;
    std::vector<tT> _tax_parent_vec;
    std::vector<uint8_t> _tax_depth_vec;
    const uint16_t _rootID = 1;
    std::vector<std::pair<tT, uint16_t>> _bases_sizes;
    encT* _enc_arr;
    tT* _tlca_arr;
    // scT* _scount_ar;
    uint8_t* _ind_arr;
    bool loadMetadata();
    bool loadTaxonomy();
    void run(uint64_t rbatch_size = DEFAULT_BATCH_SIZE);
    decltype(_npositions)& npositions() { return _npositions; }
    decltype(_positions)& positions() { return _positions; }
    decltype(_lsh_vg)& lsh_vg() { return _lsh_vg; }
  };
  bool _log;
  std::vector<std::string> _library_dirpaths;
  const char* _query_filepath;
  const char* _output_dirpath;
  uint8_t _num_libraries;
  uint8_t _max_match_hdist;
  uint8_t _save_match_info;
  std::vector<std::unique_ptr<SearchL>> _slib_ptr_v;
  std::unordered_map<std::string, std::string> _queryID_to_path;
  std::unordered_map<tT, uint64_t> _tID_to_taxID;
};

#endif
