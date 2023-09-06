#include "library.h"
#include "common.h"

Library::Library(const char* library_dirpath,
                 const char* nodes_filepath,
                 const char* input_filepath,
                 uint8_t k,
                 uint8_t h,
                 uint8_t b,
                 uint64_t capacity_size,
                 uint32_t tbatch_size,
                 bool swicth_selection,
                 bool in_library,
                 bool on_disk,
                 uint8_t specified_batch)
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
  , _switch_selection(swicth_selection)
  , _in_library(in_library)
  , _on_disk(on_disk)
  , _specified_batch(specified_batch)
{
  bool isDir = IO::ensureDirectory(_library_dirpath);
  if (!in_library) {
    if (isDir) {
      std::cout << "Library will be created at " << _library_dirpath << std::endl;
    } else {
      std::cerr << "Library can not be created at " << _library_dirpath << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "Library is at " << _library_dirpath << std::endl;
    std::cout << "Sets of k-mers are already on disk." << std::endl;
  }

  if (!_in_library) {
    getRandomPositions();
    if (_specified_batch != 0) {
      std::cerr << "A specific batch can not be built before saving processed "
                   "k-mer sets into library."
                << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    Library::loadMetadata();
    std::puts("Given parameters will be ignored. Metadata of the given library "
              "in the disk is loaded.\n");
  }
  _lsh_vg = generateMaskLSH(_positions);
  _num_rows = pow(2, 2 * _h);
  _total_batches = _num_rows / _tbatch_size;

  uint64_t root_size = 0;
  _num_species = _taxonomy_record.tID_to_input().size();
  std::vector<tT> tID_vec;
  tID_vec.reserve(_num_species);
  for (auto const& kv : _taxonomy_record.tID_to_input())
    tID_vec.push_back(kv.first);

#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < tID_vec.size(); ++i) {
    tT tID_key = tID_vec[i];
    std::string filepath = _taxonomy_record.tID_to_input()[tID_key];
#pragma omp critical
    {
      std::cout << "Processing k-mer set for taxon " << _taxonomy_record.changeIDtax(tID_key)
                << std::endl;
    }
    uint64_t size_basis = 0;
    StreamIM<encT> sIM(filepath.c_str(), _k, h, &_lsh_vg, &_npositions);
    if (!in_library) {
      size_basis = sIM.processInput(static_cast<uint64_t>(DEFAULT_BATCH_SIZE));
    } else {
      std::ifstream input(filepath);
      std::string tmp;
      while (std::getline(input, tmp))
        ++size_basis;
      size_basis /= 2;
    }
#pragma omp critical
    {
      _basis_to_size[tID_key] = size_basis;
      root_size += size_basis;
    }
    assert(on_disk || !in_library);
    if (_on_disk) {
      std::string disk_path(_library_dirpath);
      disk_path = disk_path + "/lsh_enc_vec-" + std::to_string(tID_key);
      if (!in_library) {
        bool is_ok = sIM.save(disk_path.c_str());
        if (!is_ok) {
          std::cerr << "Error saving to " << disk_path << std::endl;
          exit(EXIT_FAILURE);
        }
      }
#pragma omp critical
      {
        _streamOD_map.emplace(std::make_pair(tID_key, StreamOD<encT>(disk_path.c_str())));
        _streamOD_map.at(tID_key).openStream();
      }
      sIM.clear();
    } else {
#pragma omp critical
      {
        _streamIM_map.insert(std::make_pair(tID_key, sIM));
      }
    }
  }
  _root_size = root_size;
  std::cout << "Total number of (non-unique) k-mers under the root: " << _root_size << std::endl;
  bool taxonomy_saved = _taxonomy_record.saveTaxonomyRecord(_library_dirpath);

  std::cout << "Library is initialized and k-mers sets are ready to be streamed." << std::endl;
  if (on_disk)
    std::cout << "Stream will be on disk." << std::endl;
  else
    std::cout << "Stream will be in memory." << std::endl;

  std::cout << "Details:" << std::endl;
  std::cout << "\tLength of the k-mer, k: " << std::to_string(_k) << std::endl;
  std::cout << "\tNumber of positions of LSH, h: " << std::to_string(_h) << std::endl;
  std::cout << "\tNumber of columns in the table, b: " << std::to_string(_b) << std::endl;
  std::cout << "\tMaximum capacity size: " << _capacity_size << std::endl;
  std::cout << "\tTable row batch size: " << _tbatch_size << std::endl;
  std::cout << "\tTotal number of rows: " << _num_rows << std::endl;
  std::cout << "\tNumber of species: " << _num_species << std::endl;
  std::cout << "\tTotal number of (non-unique) kmers: " << _root_size << std::endl;
  std::cout << std::endl;

  if (_specified_batch == 0)
    Library::saveMetadata();
}

void
Library::run(uint8_t sdepth, SelectionMethod init_selection)
{
  if (_switch_selection && init_selection == random_kmer)
    std::puts("No counterpart to current selection method, no switching and "
              "only random.\n");
  std::cout << "Total number of batches is  " << _total_batches << std::endl;
  unsigned int curr_batch = 0;
  for (unsigned int i = 0; i < _total_batches; ++i) {
    curr_batch++;
    if ((_specified_batch == 0) || (_specified_batch == curr_batch)) {
      std::cout << "Building the library for batch " << curr_batch << "..." << std::endl;
      HTs<encT> ts_root(1, _k, _h, _b, _tbatch_size, &_lsh_vg, init_selection);
      Library::getBatchHTs(&ts_root, 0, sdepth);
      Library::saveBatchHTs(ts_root, curr_batch);
    } else {
      Library::continueBatch();
    }
  }
}

void
Library::continueBatch()
{
  vvec<encT> tmp;
  tmp.resize(_tbatch_size);
  if (_on_disk) {
    for (auto& kv : _streamOD_map) {
      _streamOD_map.at(kv.first).getBatch(tmp, _tbatch_size);
    }
  } else {
    for (auto& kv : _streamIM_map) {
      _streamIM_map.at(kv.first).getBatch(tmp, _tbatch_size);
    }
  }
}

void
Library::getRandomPositions()
{
  uint8_t n;
  assert(_h <= 16);
  assert(_h < _k);
  for (uint8_t m = 0; m < _h; m++) {
    n = rand() % _k;
    if (count(_positions.begin(), _positions.end(), n)) {
      m -= 1;
    } else {
      _positions.push_back(n);
    }
  }
  sort(_positions.begin(), _positions.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < _k; ++i) {
    if (i != _positions[ix_pos])
      _npositions.push_back(i);
    else
      ix_pos++;
  }
}

uint64_t
Library::getConstrainedSize(std::set<tT> tIDsBasis)
{
  uint64_t constrained_size;
  uint64_t sum_size = 0;
  for (tT tID : tIDsBasis) {
    sum_size = sum_size + _basis_to_size[tID];
  }
  float sq_ratio = sqrt(static_cast<float>(sum_size) / static_cast<float>(_root_size));
  float batch_ratio = static_cast<float>(_tbatch_size) / static_cast<float>(_num_rows);
  constrained_size = static_cast<uint64_t>(_capacity_size * sq_ratio * batch_ratio);
  return constrained_size;
}

void
Library::getBatchHTs(HTs<encT>* ts, uint8_t curr_depth, uint8_t last_depth)
{
  if ((curr_depth >= last_depth) || _taxonomy_record.isBasis(ts->tID)) {
    HTd<encT> td(ts->tID, ts->k, ts->h, ts->num_rows, ts->ptr_lsh_vg, ts->kmer_selection);
    getBatchHTd(&td);
    td.convertHTs(ts);
  } else {
    SelectionMethod next_selection;
    if (_switch_selection && (curr_depth >= (last_depth - 1))) {
      if (ts->kmer_selection == information_score) {
        next_selection = large_scount;
      } else if (ts->kmer_selection == large_scount) {
        next_selection = information_score;
      } else {
        next_selection = random_kmer;
      }
    } else {
      next_selection = ts->kmer_selection;
    }
    std::set<tT>& children_set = _taxonomy_record.child_map()[ts->tID];
    size_t num_child = children_set.size();
    std::vector<double> children(num_child);
    std::copy(children_set.begin(), children_set.end(), children.begin());
    ts->childrenHT.assign(
      num_child, HTs<encT>(0, ts->k, ts->h, ts->b, ts->num_rows, ts->ptr_lsh_vg, next_selection));
    for (unsigned int ti = 0; ti < num_child; ++ti) {
      ts->childrenHT[ti].tID = children[ti];
      getBatchHTs(&(ts->childrenHT[ti]), curr_depth + 1, last_depth);
    }
    for (auto& ts_c : ts->childrenHT) {
      ts->unionRows(ts_c);
    }
    if (ts->childrenHT.size() > 1) {
      int64_t num_rm;
      uint64_t constrained_size = getConstrainedSize(ts->tIDsBasis);
      num_rm = static_cast<int64_t>(ts->num_kmers) - static_cast<int64_t>(constrained_size);
      ts->updateLCA();
      if (num_rm > 0) {
        ts->shrinkHT(static_cast<uint64_t>(num_rm));
      }
    }
    ts->childrenHT.clear();
    std::cout << "HTs has been constructed for " << _taxonomy_record.changeIDtax(ts->tID)
              << std::endl;
  }
}

void
Library::getBatchHTd(HTd<encT>* td)
{
  if (_taxonomy_record.isBasis(td->tID)) {
    if (_on_disk) {
      _streamOD_map.at(td->tID).getBatch(td->enc_vvec, _tbatch_size);
    } else {
      _streamIM_map.at(td->tID).getBatch(td->enc_vvec, _tbatch_size);
    }
    td->initBasis(td->tID);
    td->makeUnique();
    std::cout << "Basis HTd has been constructed for " << _taxonomy_record.changeIDtax(td->tID)
              << std::endl;
  } else {
    std::set<tT>& children_set = _taxonomy_record.child_map()[td->tID];
    size_t num_child = children_set.size();
    std::vector<double> children(num_child);
    std::copy(children_set.begin(), children_set.end(), children.begin());
    td->childrenHT.assign(
      num_child, HTd<encT>(0, td->k, td->h, td->num_rows, td->ptr_lsh_vg, td->kmer_selection));
    for (unsigned int ti = 0; ti < num_child; ++ti) {
      td->childrenHT[ti].tID = children[ti];
      getBatchHTd(&(td->childrenHT[ti]));
    }
    for (auto& td_c : td->childrenHT) {
      td->unionRows(td_c);
    }
    if (td->childrenHT.size() > 1) {
      int64_t num_rm;
      uint64_t constrained_size = getConstrainedSize(td->tIDsBasis);
      num_rm = static_cast<int64_t>(td->num_kmers) - static_cast<int64_t>(constrained_size);
      td->updateLCA();
      if (num_rm > 0) {
        td->shrinkHT(static_cast<uint64_t>(num_rm), _b);
      }
    }
    td->childrenHT.clear();
    std::cout << "HTd has been constructed for " << _taxonomy_record.changeIDtax(td->tID)
              << std::endl;
  }
}

bool
Library::saveBatchHTs(HTs<encT>& ts, uint16_t curr_batch)
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);

  {
    FILE* encf =
      IO::open_file((save_dirpath + "/enc_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.enc_arr, sizeof(encT), ts.num_rows * ts.b, encf);
    if (std::ferror(encf)) {
      std::puts("I/O error when writing encoding array to the library.\n");
      is_ok = false;
    }
    std::fclose(encf);
  }

  {
    FILE* indf =
      IO::open_file((save_dirpath + "/ind_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.ind_arr, sizeof(uint8_t), ts.num_rows, indf);
    if (std::ferror(indf)) {
      std::puts("I/O error when writing indicator array to the library.\n");
      is_ok = false;
    }
    std::fclose(indf);
  }

  {
    FILE* tlcaf = IO::open_file(
      (save_dirpath + "/tlca_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.tlca_arr, sizeof(tT), ts.num_rows * ts.b, tlcaf);
    if (std::ferror(tlcaf)) {
      std::puts("I/O error when writing taxon-LCA array to the library.\n");
      is_ok = false;
    }
    std::fclose(tlcaf);
  }

  {
    FILE* scountf = IO::open_file(
      (save_dirpath + "/scount_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.scount_arr, sizeof(scT), ts.num_rows * ts.b, scountf);
    if (std::ferror(scountf)) {
      std::puts("I/O error when writing species-count array to the library.\n");
      is_ok = false;
    }
    std::fclose(scountf);
  }

  return is_ok;
}

bool
Library::loadMetadata()
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint16_t>> bases_sizes;

  FILE* metadataf = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "rb");
  std::fread(&_k, sizeof(uint16_t), 1, metadataf);
  std::fread(&_h, sizeof(uint16_t), 1, metadataf);
  std::fread(&_b, sizeof(uint16_t), 1, metadataf);
  std::fread(&_capacity_size, sizeof(uint64_t), 1, metadataf);
  std::fread(&_total_batches, sizeof(uint16_t), 1, metadataf);
  std::fread(&_tbatch_size, sizeof(uint32_t), 1, metadataf);
  std::fread(&_num_rows, sizeof(uint64_t), 1, metadataf);
  std::fread(&_num_species, sizeof(uint64_t), 1, metadataf);
  std::fread(&_root_size, sizeof(uint64_t), 1, metadataf);
  bases_sizes.resize(_num_species);
  std::fread(bases_sizes.data(), sizeof(std::pair<tT, uint16_t>), _num_species, metadataf);
  _positions.resize(_h);
  std::fread(_positions.data(), sizeof(uint8_t), _positions.size(), metadataf);
  _npositions.resize(_k - _h);
  std::fread(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadataf);

  if (std::ferror(metadataf)) {
    std::puts("I/O error when reading metadata file from the library.\n");
    is_ok = false;
  }
  std::fclose(metadataf);

  return is_ok;
}

bool
Library::saveMetadata()
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint16_t>> bases_sizes(_basis_to_size.begin(), _basis_to_size.end());

  FILE* metadataf = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "wb");
  std::fwrite(&_k, sizeof(uint16_t), 1, metadataf);
  std::fwrite(&_h, sizeof(uint16_t), 1, metadataf);
  std::fwrite(&_b, sizeof(uint16_t), 1, metadataf);
  std::fwrite(&_capacity_size, sizeof(uint64_t), 1, metadataf);
  std::fwrite(&_total_batches, sizeof(uint16_t), 1, metadataf);
  std::fwrite(&_tbatch_size, sizeof(uint32_t), 1, metadataf);
  std::fwrite(&_num_rows, sizeof(uint64_t), 1, metadataf);
  std::fwrite(&_num_species, sizeof(uint64_t), 1, metadataf);
  std::fwrite(&_root_size, sizeof(uint64_t), 1, metadataf);
  std::fwrite(bases_sizes.data(), sizeof(std::pair<tT, uint16_t>), bases_sizes.size(), metadataf);
  std::fwrite(_positions.data(), sizeof(uint8_t), _positions.size(), metadataf);
  std::fwrite(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadataf);

  if (std::ferror(metadataf)) {
    std::puts("I/O error when writing metadata file to the library.\n");
    is_ok = false;
  }
  std::fclose(metadataf);

  return is_ok;
}
