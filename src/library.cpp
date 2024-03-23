#include "library.h"
#include "common.h"

Library::Library(const char *library_dirpath,
                 const char *nodes_filepath,
                 const char *input_filepath,
                 uint8_t k,
                 uint8_t w,
                 uint8_t h,
                 uint8_t b,
                 RankingMethod ranking_method,
                 bool adaptive_size,
                 uint64_t capacity_size,
                 uint32_t tbatch_size,
                 bool from_library,
                 bool input_kmers,
                 uint16_t target_batch,
                 bool only_init,
                 bool verbose,
                 bool log)
  : _library_dirpath(library_dirpath)
  , _taxonomy_input(nodes_filepath)
  , _taxonomy_record(input_filepath, _taxonomy_input)
  , _nodes_filepath(nodes_filepath)
  , _input_filepath(input_filepath)
  , _k(k)
  , _w(w)
  , _h(h)
  , _b(b)
  , _ranking_method(ranking_method)
  , _adaptive_size(adaptive_size)
  , _capacity_size(capacity_size)
  , _tbatch_size(tbatch_size)
  , _from_library(from_library)
  , _input_kmers(input_kmers)
  , _target_batch(target_batch)
  , _only_init(only_init)
  , _verbose(verbose)
  , _log(log)
{
  _root_size = 0;
  _num_species = _taxonomy_record.tID_to_input().size();
  _num_nodes = _taxonomy_record.get_num_nodes();
  if (!_from_library) {
    if (IO::ensureDirectory(_library_dirpath)) {
      std::cout << "Library directory exists, current files will be overwritten " << _library_dirpath << std::endl;
    } else if (ghc::filesystem::create_directories(_library_dirpath)) {
      std::cout << "Library directory has been created at the given path " << _library_dirpath << std::endl;
    } else {
      std::cerr << "Library directory can not be created at the given path " << _library_dirpath << std::endl;
      exit(EXIT_FAILURE);
    }
    getRandomPositions();
  } else {
    std::cout << "Library has been found at the given path " << _library_dirpath << std::endl;
    Library::loadMetadata();
    if (_verbose)
      std::puts("Given parameters will be ignored, parameters found in library metadata will be used.\n");
  }

  _lsh_vg = generateMaskLSH(_positions);
  _num_rows = pow(2, 2 * _h);
  _total_batches = _num_rows / _tbatch_size;

  if (_target_batch > _total_batches) {
    std::cerr << "Given target " << _target_batch << " can not be greater than # of batches " << _total_batches << std::endl;
    exit(EXIT_FAILURE);
  }

  _tID_vec.reserve(_num_species);
  for (auto const &kv : _taxonomy_record.tID_to_input()) {
    _tID_vec.push_back(kv.first);
    _basis_to_ninput[kv.first] = kv.second.size();
  }
  std::sort(_tID_vec.begin(), _tID_vec.end(), [&](const tT &x, const tT &y) {
    return _taxonomy_record.tID_to_input()[x].size() > _taxonomy_record.tID_to_input()[y].size();
  });

  if (!from_library) {
    for (int i = 1; i <= _total_batches; ++i) {
      std::string batch_dirpath = _library_dirpath;
      batch_dirpath += +"/batch" + std::to_string(i);
      ghc::filesystem::create_directories(batch_dirpath);
    }
    std::string rcounts_dirpath = _library_dirpath;
    rcounts_dirpath += "/rcounts/";
    ghc::filesystem::create_directories(rcounts_dirpath);
  }

  if (_verbose)
    std::cout << "Reading and processing the input sequences..." << std::endl;

#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < _tID_vec.size(); ++i) {
    tT tID_key = _tID_vec[i];
    processLeaf(tID_key);
  }

  if (_log || _verbose) {
    std::cout << "LSH positions:";
    for (auto p : _positions)
      std::cout << " " << std::to_string(p);
    std::cout << std::endl;
  }

  if (_log || _verbose)
    LOG(INFO) << "Total number of (non-distinct) k-mers under the root: " << _root_size << std::endl;

  if (_verbose) {
    std::cout << "Library has been set and k-mers sets are ready to be streamed" << std::endl;
  }

  if (_verbose) {
    std::cout << "Details:" << std::endl;
    std::cout << "\tSequences as input directly: " << std::noboolalpha << !_input_kmers << std::endl;
    std::cout << "\tLength of the k-mer, k: " << std::to_string(_k) << std::endl;
    std::cout << "\tLength of the minimizer window, w: " << std::to_string(_w) << std::endl;
    std::cout << "\tNumber of positions of LSH, h: " << std::to_string(_h) << std::endl;
    std::cout << "\tNumber of columns in the table, b: " << std::to_string(_b) << std::endl;
    std::cout << "\tMaximum capacity size: " << _capacity_size << std::endl;
    std::cout << "\tTable row batch size: " << _tbatch_size << std::endl;
    std::cout << "\tTotal number of rows: " << _num_rows << std::endl;
    std::cout << "\tNumber of species: " << _num_species << std::endl;
    std::cout << "\tNumber of nodes: " << _num_nodes << std::endl;
    std::cout << "\tTotal number of (non-distinct) kmers: " << _root_size << std::endl;
    std::cout << std::endl;
  }

  if (_log)
    LOG(INFO) << "Specified batch is " << static_cast<uint16_t>(_target_batch) << std::endl;

  if (_target_batch == 0 || !_from_library) {
    Library::saveMetadata();
    bool taxonomy_saved = _taxonomy_record.saveTaxonomyRecord(_library_dirpath);
    if (_log)
      LOG(INFO) << "Metadata has been saved to the library" << std::endl;
  } else {
    if (_log)
      LOG(INFO) << "Metadata has been read from the library, won't be saved again" << std::endl;
  }

  if (!_only_init)
    Library::build();
  /* Library::annotateInfo(); */
}

void Library::processLeaf(tT tID_key)
{
  std::vector<std::string> &filepath_v = _taxonomy_record.tID_to_input()[tID_key];
  if (_log) {
#pragma omp critical
    LOG(INFO) << "Processing k-mer set for taxon " << _taxonomy_record.changeIDtax(tID_key) << std::endl;
  }
  uint64_t size_basis = 0;
  inputHandler<encT> pI(filepath_v, _k, _w, _h, &_lsh_vg, &_npositions);
  if (!_from_library) {
    if (!pI.checkInput(_library_dirpath, tID_key, _total_batches)) {
      if (_input_kmers)
        size_basis = pI.readInput(static_cast<uint64_t>(DEFAULT_BATCH_SIZE));
      else
        size_basis = pI.extractInput(1);
    } else {
      bool is_ok = pI.loadInput(_library_dirpath, tID_key, _total_batches);
      size_basis = pI.lsh_enc_vec.size();
      if (_log) {
#pragma omp critical
        LOG(NOTICE) << "Encodings and hash values already exist for " << _taxonomy_record.changeIDtax(tID_key) << std::endl;
      }
    }
#pragma omp critical
    {
      _basis_to_size[tID_key] = size_basis;
      tT tmp_tID = tID_key;
      while (_taxonomy_record.parent_vec()[tmp_tID] != 0) {
        _tID_to_size[tmp_tID] += size_basis;
        tmp_tID = _taxonomy_record.parent_vec()[tmp_tID];
      }
      _root_size += size_basis;
      _tID_to_size[1] += size_basis;
    }
    bool is_ok = pI.saveInput(_library_dirpath, tID_key, _total_batches, _tbatch_size);
    if (!is_ok) {
      std::cerr << "Error saving LSH-value and encoding pairs for " << _taxonomy_record.changeIDtax(tID_key) << std::endl;
      exit(EXIT_FAILURE);
    }
    if (_log) {
#pragma omp critical
      LOG(INFO) << "LSH-value and encoding pairs has been saved for " << _taxonomy_record.changeIDtax(tID_key) << std::endl;
    }
  } else {
#pragma omp critical
    {
      _inputStream_map.emplace(std::make_pair(tID_key, inputStream<encT>(_library_dirpath, tID_key)));
    }
  }
  pI.clearInput();
}

void Library::countBasis(HTs<encT> &ts, unsigned int curr_batch)
{
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < _tID_vec.size(); ++i) {
    tT tID_key = _tID_vec[i];
    std::unordered_map<encT, uint64_t> rcounts;
    _inputStream_map.at(tID_key).loadCounts(rcounts);
    std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;
    _inputStream_map.at(tID_key).loadBatch(lsh_enc_vec, curr_batch);
    for (unsigned int i = 0; i < lsh_enc_vec.size(); ++i) {
      uint32_t rix = lsh_enc_vec[i].first - (curr_batch - 1) * _tbatch_size;
      for (unsigned int j = 0; j < ts.ind_arr[rix]; ++j) {
        if (ts.enc_arr[rix * _b + j] == lsh_enc_vec[i].second) {
#pragma omp atomic update
          ts.scount_arr[rix * _b + j] += (1 + rcounts[lsh_enc_vec[i].second]);
        }
      }
    }
  }
}

void Library::resetInfo(HTs<encT> &ts, bool reset_scount, bool reset_tlca)
{
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (uint32_t rix = 0; rix < ts.num_rows; ++rix) {
    for (unsigned int j = 0; j < ts.ind_arr[rix]; ++j) {
      if (reset_scount)
        ts.scount_arr[rix * _b + j] = 0;
      if (reset_tlca)
        ts.tlca_arr[rix * _b + j] = 0;
    }
  }
}

void Library::softLCA(HTs<encT> &ts, unsigned int curr_batch)
{
  double w = 6.0;
  unsigned int r = 5;
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < _tID_vec.size(); ++i) {
    tT tID_key = _tID_vec[i];
    std::unordered_map<encT, uint64_t> rcounts;
    _inputStream_map.at(tID_key).loadCounts(rcounts);
    std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;
    _inputStream_map.at(tID_key).loadBatch(lsh_enc_vec, curr_batch);
    double p_update, sc;
    for (unsigned int i = 0; i < lsh_enc_vec.size(); ++i) {
      uint32_t rix = lsh_enc_vec[i].first - (curr_batch - 1) * _tbatch_size;
      for (unsigned int j = 0; j < ts.ind_arr[rix]; ++j) {
        if ((ts.enc_arr[rix * _b + j] == lsh_enc_vec[i].second)) { // && nseen[rix * _b + j]
          sc = static_cast<double>(rcounts[lsh_enc_vec[i].second]) + 1.0;
          if (ts.scount_arr[rix * _b + j] <= r)
            p_update = 1.0;
          else
            p_update = 1 - pow((1 - 1 / log2(pow((ts.scount_arr[rix * _b + j] - 1) / w, 2) + 2.0)), sc);
          std::bernoulli_distribution bt(p_update);
#pragma omp critical
          {
            if ((p_update == 1.0) || bt(gen)) {
              ts.tlca_arr[rix * _b + j] = _taxonomy_record.getLowestCommonAncestor(ts.tlca_arr[rix * _b + j], tID_key);
            }
          }
        }
      }
    }
  }
}

void Library::annotateInfo()
{
  if (_verbose)
    std::cout << "Total number of batches is " << _total_batches << std::endl;
  unsigned int curr_batch = 0;
  for (unsigned int i = 0; i < _total_batches; ++i) {
    curr_batch++;
    if ((_target_batch == 0) || (_target_batch == curr_batch)) {
      std::cout << "Annotating k-mners in the library for batch " << curr_batch << "..." << std::endl;
      HTs<encT> ts_root(1, _k, _h, _b, _tbatch_size, &_lsh_vg, _ranking_method, &_taxonomy_record);
      Library::loadBatchHTs(ts_root, curr_batch);
      if (_log)
        LOG(INFO) << "The table has been loaded and contains " << ts_root.num_kmers << "/" << ts_root.num_rows * ts_root.b
                  << " k-mers" << std::endl;
      Library::resetInfo(ts_root, true, true);
      Library::countBasis(ts_root, curr_batch);
      Library::softLCA(ts_root, curr_batch);
      if (_log)
        LOG(INFO) << "Leaves are counted and k-mers are labeled with soft-LCA" << std::endl;
      Library::saveBatchHTs(ts_root, curr_batch);
    }
  }
}

void Library::build()
{
  if (_verbose)
    std::cout << "Total number of batches is " << _total_batches << std::endl;
  unsigned int curr_batch = 0;
  for (unsigned int i = 0; i < _total_batches; ++i) {
    curr_batch++;
    if ((_target_batch == 0) || (_target_batch == curr_batch)) {
      std::cout << "Building the library for batch " << curr_batch << "..." << std::endl;
      HTs<encT> ts_root(1, _k, _h, _b, _tbatch_size, &_lsh_vg, _ranking_method, &_taxonomy_record);
      Library::getBatchHTs(&ts_root, curr_batch);
      if (_log)
        LOG(INFO) << "The table has been built and contains " << ts_root.num_kmers << "/" << ts_root.num_rows * ts_root.b
                  << " k-mers" << std::endl;
      Library::resetInfo(ts_root, true, true);
      Library::countBasis(ts_root, curr_batch);
      Library::softLCA(ts_root, curr_batch);
      if (_log)
        LOG(INFO) << "Leaves are counted and k-mers are labeled with soft-LCA" << std::endl;
      Library::saveBatchHTs(ts_root, curr_batch);
    }
  }
}

void Library::getRandomPositions()
{
  uint8_t n;
  assert(_h <= 16);
  assert(_h < _k);
  std::uniform_int_distribution<uint8_t> distrib(0, _k - 1);
  for (uint8_t m = 0; m < _h; m++) {
    n = distrib(gen);
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

uint64_t Library::getConstrainedSizeKC(tT curr_tID)
{
  uint64_t constrained_size;
  uint64_t sum_size = _tID_to_size[curr_tID];
  /* uint64_t kingdom_size = _tID_to_size[_taxonomy_record.tID_to_lsroot()[curr_tID]]; */
  /* float sq_ratio = sqrt(static_cast<float>(sum_size) / static_cast<float>(kingdom_size)); */
  float sq_ratio = sqrt(static_cast<float>(sum_size) / static_cast<float>(_root_size));
  float batch_ratio = static_cast<float>(_tbatch_size) / static_cast<float>(_num_rows);
  constrained_size = static_cast<uint64_t>(_capacity_size * sq_ratio * batch_ratio);
  return constrained_size;
}

uint64_t Library::getConstrainedSizeSC(uint64_t num_basis)
{
  float sq_ratio = sqrt(static_cast<float>(num_basis) / static_cast<float>(_num_species));
  float batch_ratio = static_cast<float>(_tbatch_size) / static_cast<float>(_num_rows);
  uint64_t constrained_size = static_cast<uint64_t>(_capacity_size * sq_ratio * batch_ratio);
  return constrained_size;
}

void Library::getBatchHTs(HTs<encT> *ts, unsigned int curr_batch)
{
  if (_verbose)
    std::cout << "Constructing the table for " << _taxonomy_record.changeIDtax(ts->tID) << std::endl;
  HTd<encT> td(ts->tID, ts->k, ts->h, ts->num_rows, ts->ptr_lsh_vg, _ranking_method, &_taxonomy_record);
  omp_set_nested(1);
  num_tasks = static_cast<unsigned int>(sqrt(2 * num_threads));
  if (num_tasks < LNUM_TASKS)
    num_tasks = std::min(static_cast<unsigned int>(LNUM_TASKS), num_threads);
  omp_set_num_threads(num_tasks);
#pragma omp parallel
  {
#pragma omp single
    {
      getBatchHTd(&td, curr_batch);
    }
  }
#pragma omp taskwait
  if (_log)
    LOG(INFO) << "Converting from HTd to HTs for the taxon " << _taxonomy_record.changeIDtax(ts->tID) << std::endl;
  td.convertHTs(ts);
  if (_log) {
    LOG(INFO) << "The static hash table has been constructed for " << _taxonomy_record.changeIDtax(ts->tID) << std::endl;
    std::cout << "Histogram for # of table columns:" << std::endl;
    for (auto kv : ts->histRowSizes())
      std::cout << "\t" << static_cast<unsigned int>(kv.first) << " : " << kv.second << std::endl;
  }
}

void Library::getBatchHTd(HTd<encT> *td, unsigned int curr_batch)
{
  uint64_t curr_taxID = _taxonomy_record.changeIDtax(td->tID);
  uint64_t num_batch_kmers;
  if (_verbose)
    std::cout << "Constructing the table for " << curr_taxID << std::endl;
  if (_taxonomy_record.isBasis(td->tID)) {
    num_batch_kmers = _inputStream_map.at(td->tID).retrieveBatch(td->enc_vvec, _tbatch_size, curr_batch);
    td->initBasis(td->tID);
    /* td->makeUnique(true); */
    if (_log) {
#pragma omp critical
      {
        LOG(INFO) << "The dynamic hash table has been constructed for the leaf " << curr_taxID << std::endl;
        LOG(INFO) << "The number of k-mers read for this batch is " << num_batch_kmers << "/" << _basis_to_size[td->tID]
                  << std::endl;
      }
    }
  } else {
    std::set<tT> &children_set = _taxonomy_record.child_map()[td->tID];
    size_t num_child = children_set.size();
    std::vector<tT> children(num_child);
    std::copy(children_set.begin(), children_set.end(), children.begin());
    td->childrenHT.assign(num_child,
                          HTd<encT>(0, td->k, td->h, td->num_rows, td->ptr_lsh_vg, td->ranking_method, td->taxonomy_record));
    for (unsigned int ti = 0; ti < num_child; ++ti) {
      if (_log)
#pragma omp critical
        LOG(INFO) << "Building for the child " << _taxonomy_record.changeIDtax(children[ti]) << " of " << curr_taxID
                  << std::endl;
      td->childrenHT[ti].tID = children[ti];
#pragma omp task untied
      {
        getBatchHTd(&(td->childrenHT[ti]), curr_batch);
      }
    }
#pragma omp taskwait

    if (_log)
      LOG(INFO) << "Taking the union of tables of children of " << curr_taxID << std::endl;

    if (_taxonomy_record.depth_vec()[td->tID] <= LSR)
      omp_set_num_threads(num_threads);
    else
      omp_set_num_threads((num_threads + num_tasks - 1) / num_tasks + 1);

    for (auto &td_c : td->childrenHT) {
      td->unionRows(td_c, false);
    }
    td->updateSize();

    if (td->childrenHT.size() > 1) {
      /* if (td->tID != 0 && td->tID != 1) { */
      /*   if (_log) */
      /*     LOG(INFO) << curr_taxID << " has more than one child, updating LCA labels" << std::endl; */
      /*   td->updateLCA(); // computes hard-LCA labels along the way */
      /* } */

      /* if (_taxonomy_record.depth_vec()[td->tID] <= LSR) { */
      /*   td->filterLSR(_taxonomy_record.depth_vec(), LSR); // removes k-mers with hard-LCA above LSR rank */
      /* } */

      uint64_t constraint_size = getConstrainedSizeKC(td->tID);
      /* uint64_t constraint_size = getConstrainedSizeSC(td->num_basis); */
      int64_t num_rm = static_cast<int64_t>(td->num_kmers) - static_cast<int64_t>(constraint_size);
      if ((num_rm > 0) && _adaptive_size) {
        if (_log)
#pragma omp critical
          LOG(INFO) << curr_taxID << " has " << num_rm << " k-mers above the constraint" << std::endl;
        td->shrinkHT(static_cast<uint64_t>(num_rm), _b);
      }
    }

    td->childrenHT.clear();
    if (_log)
#pragma omp critical
      LOG(INFO) << "The dynamic hash table has been constructed for " << curr_taxID << std::endl;
  }
  if (_log) {
    std::cout << "Histogram for # of table columns:" << std::endl;
    for (auto kv : td->histRowSizes())
      std::cout << "\t" << kv.first << " : " << kv.second << std::endl;
  }
}

bool Library::loadBatchHTs(HTs<encT> &ts, unsigned int curr_batch)
{
  bool is_ok = true;
  std::string load_dirpath(_library_dirpath);

  if (_log)
    LOG(INFO) << "Loading the library for batch " << curr_batch << std::endl;

  FILE *encf = IO::open_file((load_dirpath + "/enc_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(encf)) {
    std::puts("I/O error when reading encoding array from the library\n");
    is_ok = false;
  } else
    std::fread(ts.enc_arr, sizeof(encT), _tbatch_size * _b, encf);
  std::fclose(encf);

  FILE *indf = IO::open_file((load_dirpath + "/ind_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(indf)) {
    std::puts("I/O error when reading indicator array from the library\n");
    is_ok = false;
  } else
    std::fread(ts.ind_arr, sizeof(uint8_t), _tbatch_size, indf);
  std::fclose(indf);

  FILE *tlcaf = IO::open_file((load_dirpath + "/tlca_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(tlcaf)) {
    std::puts("I/O error when reading taxon-LCA array from the library\n");
    is_ok = false;
  } else
    std::fread(ts.tlca_arr, sizeof(tT), _tbatch_size * _b, tlcaf);
  std::fclose(tlcaf);

  FILE *scountf = IO::open_file((load_dirpath + "/scount_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(scountf)) {
    std::puts("I/O error when reading species-count array from the library\n");
    is_ok = false;
  } else
    std::fread(ts.scount_arr, sizeof(tT), _tbatch_size * _b, scountf);
  std::fclose(scountf);

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded the library for batch " << curr_batch << std::endl;
    else
      LOG(ERROR) << "Failed loading the library for batch " << curr_batch << std::endl;
  }

  return is_ok;
}

bool Library::saveBatchHTs(HTs<encT> &ts, unsigned int curr_batch)
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);

  if (_verbose)
    std::cout << "Saving the library batch " << curr_batch << std::endl;

  {
    FILE *encf = IO::open_file((save_dirpath + "/enc_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.enc_arr, sizeof(encT), ts.num_rows * ts.b, encf);
    if (std::ferror(encf)) {
      std::puts("I/O error when writing encoding array to the library.\n");
      is_ok = false;
    }
    std::fclose(encf);
  }

  {
    FILE *indf = IO::open_file((save_dirpath + "/ind_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.ind_arr, sizeof(uint8_t), ts.num_rows, indf);
    if (std::ferror(indf)) {
      std::puts("I/O error when writing indicator array to the library.\n");
      is_ok = false;
    }
    std::fclose(indf);
  }

  {
    FILE *tlcaf = IO::open_file((save_dirpath + "/tlca_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.tlca_arr, sizeof(tT), ts.num_rows * ts.b, tlcaf);
    if (std::ferror(tlcaf)) {
      std::puts("I/O error when writing taxon-LCA array to the library.\n");
      is_ok = false;
    }
    std::fclose(tlcaf);
  }

  {
    FILE *scountf = IO::open_file((save_dirpath + "/scount_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.scount_arr, sizeof(scT), ts.num_rows * ts.b, scountf);
    if (std::ferror(scountf)) {
      std::puts("I/O error when writing species-count array to the library.\n");
      is_ok = false;
    }
    std::fclose(scountf);
  }

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully saved the library for batch " << curr_batch << std::endl;
    else
      LOG(ERROR) << "Failed saving the library for batch " << curr_batch << std::endl;
  }

  return is_ok;
}

bool Library::loadMetadata()
{
  bool is_ok = true;

  if (_verbose)
    std::cout << "Loading metadata of the library" << std::endl;
  std::string save_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint64_t>> bases_sizes;
  std::vector<std::pair<tT, uint64_t>> tIDs_sizes;

  FILE *metadata_f = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "rb");
  std::fread(&_k, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_h, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_b, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_capacity_size, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_total_batches, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_tbatch_size, sizeof(uint32_t), 1, metadata_f);
  std::fread(&_num_rows, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_num_species, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_num_nodes, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_root_size, sizeof(uint64_t), 1, metadata_f);

  bases_sizes.resize(_num_species);
  tIDs_sizes.resize(_num_nodes);
  _positions.resize(_h);
  _npositions.resize(_k - _h);

  std::fread(bases_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_species, metadata_f);
  std::fread(tIDs_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_nodes, metadata_f);
  std::fread(_positions.data(), sizeof(uint8_t), _positions.size(), metadata_f);
  std::fread(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadata_f);

  if (std::ferror(metadata_f)) {
    std::puts("I/O error when reading metadata file from the library.\n");
    is_ok = false;
  }
  std::fclose(metadata_f);

  for (auto kv : bases_sizes) {
    _basis_to_size[kv.first] = kv.second;
  }
  for (auto kv : tIDs_sizes) {
    _tID_to_size[kv.first] = kv.second;
  }

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded library metadata" << std::endl;
    else
      LOG(ERROR) << "Failed loading library metadata" << std::endl;
  }

  return is_ok;
}

bool Library::saveMetadata()
{
  bool is_ok = true;
  if (_verbose)
    std::cout << "Saving metadata of the library" << std::endl;
  std::string save_dirpath(_library_dirpath);

  std::vector<std::pair<tT, uint64_t>> bases_sizes(_basis_to_size.begin(), _basis_to_size.end());
  std::vector<std::pair<tT, uint64_t>> tIDs_sizes(_tID_to_size.begin(), _tID_to_size.end());
  assert(_basis_to_size.size() == _num_species);
  assert(_tID_to_size.size() == _num_nodes);

  FILE *metadata_f = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "wb");
  std::fwrite(&_k, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_h, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_b, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_capacity_size, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_total_batches, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_tbatch_size, sizeof(uint32_t), 1, metadata_f);
  std::fwrite(&_num_rows, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_num_species, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_num_nodes, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_root_size, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(bases_sizes.data(), sizeof(std::pair<tT, uint64_t>), bases_sizes.size(), metadata_f);
  std::fwrite(tIDs_sizes.data(), sizeof(std::pair<tT, uint64_t>), tIDs_sizes.size(), metadata_f);
  std::fwrite(_positions.data(), sizeof(uint8_t), _positions.size(), metadata_f);
  std::fwrite(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadata_f);

  if (std::ferror(metadata_f)) {
    std::puts("I/O error when writing metadata file to the library.\n");
    is_ok = false;
  }
  std::fclose(metadata_f);

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully saved library metadata" << std::endl;
    else
      LOG(ERROR) << "Failed saving library metadata" << std::endl;
  }

  return is_ok;
}
