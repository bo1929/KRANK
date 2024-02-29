#ifndef _QUERY_H
#define _QUERY_H

#include "common.h"
#include "encode.h"
#include "io.h"
#include "library.h"
#include "lsh.h"
#include "table.h"
#include "taxonomy.h"

class Query
{
public:
  Query(std::vector<std::string> library_dirpaths,
        const char *output_dirpath,
        const char *query_filepath,
        uint8_t max_match_hdist = 5,
        bool save_match_info = true,
        bool verbose = false,
        bool log = false);
  void perform(uint64_t rbatch_size = DEFAULT_BATCH_SIZE);
  void processBatch(std::vector<sseq_t> seqBatch,
                    vvec_string names_vec,
                    vvec_uint64 tlca_vec_or,
                    vvec_uint64 tlca_vec_rc,
                    vvec_uint8 hdist_vec_or,
                    vvec_uint8 hdist_vec_rc);
  uint64_t getLowestCommonAncestor(uint64_t a, uint64_t b)
  {
    if (!a || !b) // LCA(x,0) = LCA(0,x) = x
      return a ? a : b;
    while (a != b) {
      if (_depth_inmap[a] < _depth_inmap[b])
        b = _parent_inmap[b];
      else
        a = _parent_inmap[a];
    }
    return a;
  }

private:
  struct QLibrary
  {
    QLibrary(const char *library_dirpath, bool log = false);
    const char *_library_dirpath;
    const uint16_t _rootID = 1;
    uint8_t _k;
    uint8_t _h;
    uint8_t _b;
    uint64_t _num_rows; // Total number of rows: pow(2, 2h).
    uint64_t _capacity_size;
    uint32_t _tbatch_size;
    uint16_t _total_batches;
    uint64_t _num_species;
    uint64_t _num_nodes;
    uint64_t _root_size;
    maskLSH _lsh_vg;
    std::vector<uint8_t> _positions;
    std::vector<uint8_t> _npositions;
    std::unordered_map<tT, uint64_t> _tID_to_taxID;
    std::unordered_map<uint64_t, uint64_t> _parent_inmap;
    tT _tax_num_input;
    tT _tax_num_nodes;
    uint64_t _tax_full_size;
    std::vector<tT> _tax_parent_vec;
    std::vector<uint8_t> _tax_depth_vec;
    std::vector<std::pair<tT, uint64_t>> _bases_sizes;
    std::vector<std::pair<tT, uint64_t>> _tIDs_sizes;
    encT *_enc_arr;
    tT *_tlca_arr;
    // scT* _scount_ar;
    uint8_t *_ind_arr;
    bool _log;
    bool loadMetadata();
    bool loadTaxonomy();
    decltype(_npositions) &npositions() { return _npositions; }
    decltype(_positions) &positions() { return _positions; }
    decltype(_lsh_vg) &lsh_vg() { return _lsh_vg; }
  };
  bool _log;
  std::vector<std::string> _library_dirpaths;
  const char *_query_filepath;
  const char *_output_dirpath;
  uint8_t _k;
  uint8_t _num_libraries;
  uint8_t _max_match_hdist;
  uint8_t _save_match_info;
  uint64_t _mask_bp;
  uint64_t _mask_lr;
  std::vector<std::unique_ptr<QLibrary>> _slib_ptr_v;
  std::unordered_map<std::string, std::string> _queryID_to_path;
  std::unordered_map<uint64_t, uint64_t> _parent_inmap;
  std::unordered_map<uint64_t, uint8_t> _depth_inmap;
};

#endif
