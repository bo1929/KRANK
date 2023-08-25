#include "library.h"
#include <utility>

Library::Library(char* library_dirpath,
                 char* nodes_filepath,
                 char* input_filepath,
                 uint8_t k,
                 uint8_t h,
                 uint8_t b,
                 uint64_t capacity_size,
                 uint32_t tbatch_size,
                 bool onDisk)
  : _library_dirpath(library_dirpath)
  , _taxonomy_ncbi(nodes_filepath)
  , _taxonomy_record(input_filepath, _taxonomy_ncbi)
  , _nodes_filepath(nodes_filepath)
  , _input_filepath(input_filepath)
  , _k(k)
  , _h(h)
  , _b(b)
  , _capacity_size(capacity_size)
  , _tbatch_size(tbatch_size)
  , _onDisk(onDisk)
{
  bool isDir = IO::ensureDirectory(_library_dirpath);
  if (isDir) {
    std::cout << "Library will be created at " << _library_dirpath << std::endl;
  } else {
    std::cerr << "Library can not be created at " << _library_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }

  _lsh_vg = generateMaskLSH(_k, _h);
  _num_rows = pow(2, 2 * _h);

  uint64_t root_size = 0;

  for (auto kv : _taxonomy_record.tID_to_input()) {
    StreamIM<encT> sIM(kv.second.c_str(), _k, h, &_lsh_vg);
    uint64_t size_basis = sIM.processInput(static_cast<uint32_t>(DEFAULT_BATCH_SIZE));
    _basis_to_size[kv.first] = size_basis;
    root_size += size_basis;
    if (_onDisk) {
      std::string save_filepath(_library_dirpath);
      save_filepath = save_filepath + "/" + std::to_string(kv.first);
      bool is_ok = sIM.save(save_filepath.c_str());
      if (!is_ok) {
        std::cerr << "Error saving to " << library_dirpath << std::endl;
        exit(EXIT_FAILURE);
      }
      _streamOD_map.emplace(std::make_pair(kv.first, StreamOD<encT>(save_filepath.c_str())));
      _streamOD_map.at(kv.first).openStream();
    } else {
      _streamIM_map.insert(std::make_pair(kv.first, sIM));
    }
  }
  _root_size = root_size;
}

uint64_t
Library::getConstrainedSize(std::set<tT> tIDsBasis)
{
  uint64_t constrained_size;
  uint64_t sum_size;
  for (tT tID : tIDsBasis) {
    sum_size = sum_size + _basis_to_size[tID];
  }
  float sq_ratio = sqrt(static_cast<float>(sum_size) / static_cast<float>(_root_size));
  float batch_ratio = static_cast<float>(_tbatch_size) / static_cast<float>(_num_rows);
  constrained_size = static_cast<uint64_t>(_capacity_size * sq_ratio * batch_ratio);
  return constrained_size;
}

void
Library::getBatchHTs(HTs<encT>& ts, uint8_t curr_depth, uint8_t last_depth)
{
  if (_taxonomy_record.isBasis(ts.tID)) {
    HTd<encT> td(ts.tID, ts.k, ts.h, ts.num_rows, ts.ptr_lsh_vg, ts.kmer_priority);
    if (_onDisk) {
      _streamOD_map.at(td.tID).getBatch(td.enc_vvec, _tbatch_size);
    } else {
      _streamIM_map.at(td.tID).getBatch(td.enc_vvec, _tbatch_size);
    }
    td.initBasis(td.tID);
    td.convertHTs(ts);
  } else if (curr_depth <= last_depth) {
    HTd<encT> td(ts.tID, ts.k, ts.h, ts.num_rows, ts.ptr_lsh_vg, ts.kmer_priority);
    getBatchHTd(td);
    td.convertHTs(ts);
  } else {
    for (tT tID_c : _taxonomy_record.child_map()[ts.tID]) {
      HTs<encT> ts_c(ts.tID, ts.k, ts.b, ts.h, ts.num_rows, ts.ptr_lsh_vg, ts.kmer_priority);
      ts.childrenHT.push_back(ts_c);
      getBatchHTs(ts.childrenHT.back(), curr_depth - 1, last_depth);
    }
    for (auto& ts_c : ts.childrenHT) {
      ts.unionRows(ts_c);
    }
    if (ts.tID != _rootID) {
      int64_t num_rm;
      uint64_t constrained_size = getConstrainedSize(ts.tIDsBasis);
      num_rm = static_cast<int64_t>(constrained_size) - static_cast<int64_t>(ts.num_kmers);
      if (num_rm > 0) {
        ts.shrinkHT(static_cast<uint64_t>(num_rm));
      }
      ts.childrenHT.clear();
    }
  }
}

void
Library::getBatchHTd(HTd<encT>& td)
{
  if (_taxonomy_record.isBasis(td.tID)) {
    if (_onDisk) {
      _streamOD_map.at(td.tID).getBatch(td.enc_vvec, _tbatch_size);
    } else {
      _streamIM_map.at(td.tID).getBatch(td.enc_vvec, _tbatch_size);
    }
    td.initBasis(td.tID);
  } else {
    for (tT tID_c : _taxonomy_record.child_map()[td.tID]) {
      HTd<encT> td_c(td.tID, td.k, td.h, td.num_rows, td.ptr_lsh_vg, td.kmer_priority);
      td.childrenHT.push_back(td_c);
      getBatchHTd(td.childrenHT.back());
    }
    for (auto& td_c : td.childrenHT) {
      td.unionRows(td_c);
    }
    int64_t num_rm;
    uint64_t constrained_size = getConstrainedSize(td.tIDsBasis);
    num_rm = static_cast<int64_t>(constrained_size) - static_cast<int64_t>(td.num_kmers);
    if (num_rm > 0) {
      td.shrinkHT(static_cast<uint64_t>(num_rm), _b);
    }
    td.childrenHT.clear();
  }
}
