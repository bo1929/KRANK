#include "library.h"

Library::Library(const char *library_dirpath,
                 const char *tax_dirpath,
                 const char *input_filepath,
                 uint8_t k,
                 uint8_t w,
                 uint8_t h,
                 uint8_t b,
                 RankingMethod ranking_method,
                 LabelsLCA labels_lca,
                 bool adaptive_size,
                 uint64_t capacity_size,
                 uint32_t tbatch_size,
                 bool from_library,
                 bool input_kmers,
                 uint16_t target_batch,
                 bool only_init,
                 bool update_annotations,
                 bool fast_mode,
                 bool remove_intermediate,
                 bool verbose,
                 bool log)
  : _library_dirpath(library_dirpath)
  , _tax_dirpath(tax_dirpath)
  , _tax_input(tax_dirpath)
  , _tax_record(input_filepath, _tax_input)
  , _input_filepath(input_filepath)
  , _k(k)
  , _w(w)
  , _h(h)
  , _b(b)
  , _ranking_method(ranking_method)
  , _labels_lca(labels_lca)
  , _adaptive_size(adaptive_size)
  , _capacity_size(capacity_size)
  , _tbatch_size(tbatch_size)
  , _from_library(from_library)
  , _input_kmers(input_kmers)
  , _target_batch(target_batch)
  , _only_init(only_init)
  , _update_annotations(update_annotations)
  , _fast_mode(fast_mode)
  , _remove_intermediate(remove_intermediate)
  , _verbose(verbose)
  , _log(log)
{
  if (_only_init)
    std::cout << "Initializing the library..." << std::endl;
  else
    std::cout << "Building the library..." << std::endl;

  _root_size = 0;
  _num_species = _tax_record.trID_to_input().size();
  _num_nodes = _tax_record.get_num_nodes();
  _num_rows = pow(2, 2 * _h);
  _total_batches = _num_rows / _tbatch_size;

  if (!_from_library) {
    getRandomPositions();
    if (ghc::filesystem::exists(_library_dirpath)) {
      std::cout << "Library directory exists, current files will be overwritten " << _library_dirpath << std::endl;
    } else if (ghc::filesystem::create_directories(_library_dirpath)) {
      std::cout << "Library directory has been created at the given path " << _library_dirpath << std::endl;
    } else {
      std::cerr << "Library directory cannot be created at the given path " << _library_dirpath << std::endl;
      exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= _total_batches; ++i) {
      // TODO: An alternative to this would be creating nested directories based on
      // hash or modulo to prevent having too many files in a single directory.
      std::string batch_dirpath = _library_dirpath;
      batch_dirpath += +"/batch" + std::to_string(i);
      ghc::filesystem::create_directories(batch_dirpath);
    }
    std::string rcounts_dirpath = _library_dirpath;
    rcounts_dirpath += "/rcounts/";
    ghc::filesystem::create_directories(rcounts_dirpath);
  } else {
    std::cout << "Library has been found at the given path " << _library_dirpath << std::endl;
    Library::loadMetadata();
    if (_verbose)
      std::puts("Given parameters will be ignored, parameters found in library metadata will be used");
  }

  _lsh_vg = generateMaskLSH(_positions);

  _trID_vec.reserve(_num_species);
  tT tmp_trID;
  for (auto const &kv : _tax_record.trID_to_input()) {
    _trID_vec.push_back(kv.first);
    _basis_to_ninput[kv.first] = kv.second.size();
    tmp_trID = kv.first;
    while (tmp_trID != 0) {
      _trID_to_ngenomes[tmp_trID] += kv.second.size();
      tmp_trID = _tax_record.parent_vec()[tmp_trID];
    }
  }
  std::sort(_trID_vec.begin(), _trID_vec.end(), [&](const tT &x, const tT &y) {
    return _tax_record.trID_to_input()[x].size() > _tax_record.trID_to_input()[y].size();
  });

  if (_verbose)
    std::cout << "Reading and processing the input sequences..." << std::endl;
  omp_set_nested(1);
  omp_set_num_threads(num_threads);
#pragma omp parallel
  {
#pragma omp single
    {
      for (unsigned int i = 0; i < _trID_vec.size(); ++i) {
        tT trID_key = _trID_vec[i];
#pragma omp task untied
        {
          processLeaf(trID_key);
        }
      }
    }
#pragma omp taskwait
  }

  if (_verbose) {
    std::cout << "LSH positions:";
    for (auto p : _positions)
      std::cout << " " << std::to_string(p);
    std::cout << std::endl;
  }
  LOG(INFO) << COND(_log) << "Total number of (non-distinct) k-mers under the root: " << _root_size << std::endl;
  if (_verbose) {
    std::cout << "Input references has been successfully processed" << std::endl;
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
    /* std::cout << "\tTotal number of (non-distinct) kmers: " << _root_size << std::endl; */
  }
  if (_verbose && !_only_init)
    std::cout << "Specified batch is " << static_cast<uint16_t>(_target_batch) << std::endl;

  if (_target_batch == 0 || !_from_library) {
    Library::saveMetadata();
    bool taxonomy_saved = _tax_record.saveTaxonomyRecord(_library_dirpath);
    LOG(INFO) << COND(_log) << "Metadata has been saved to the library" << std::endl;
  } else {
    LOG(INFO) << COND(_log) << "Metadata has been read from the library, won't be saved again" << std::endl;
  }

  if (!_only_init) {
    if (_update_annotations) {
      Library::annotateInfo();
      std::cout << "LCA annotations of k-mers has been updated and saved" << std::endl;
    } else {
      Library::buildTables();
      std::cout << "Library has been built and saved" << std::endl;
    }
    if (_remove_intermediate) {
      for (int i = 1; i <= _total_batches; ++i) {
        if ((_target_batch == 0) || (i == _target_batch)) {
          std::string batch_dirpath = _library_dirpath;
          batch_dirpath += +"/batch" + std::to_string(i);
          ghc::filesystem::remove_all(batch_dirpath);
        }
      }
      std::string check_dirpath(_library_dirpath);
      bool all_built = true;
      for (int i = 1; i <= _total_batches; ++i) {
        if (!(ghc::filesystem::exists(check_dirpath + "/enc_arr-" + std::to_string(i)) &&
              ghc::filesystem::exists(check_dirpath + "/scount_arr-" + std::to_string(i)) &&
              ghc::filesystem::exists(check_dirpath + "/ind_arr-" + std::to_string(i)) &&
              ghc::filesystem::exists(check_dirpath + "/tlca_arr-" + std::to_string(i))))
          all_built = false;
      }
      if (all_built) {
        std::string rcounts_dirpath = _library_dirpath;
        rcounts_dirpath += "/rcounts/";
        ghc::filesystem::remove_all(rcounts_dirpath);
      }
    }
  } else {
    std::cout << "Library has been successfully initialized" << std::endl;
  }
}

void Library::processLeaf(tT trID_key)
{
  std::vector<std::string> &filepath_v = _tax_record.trID_to_input()[trID_key];
#pragma omp critical(log_processing)
  LOG(INFO) << COND(_log) << "Processing k-mer set for taxon " << _tax_record.tiID_from_trID(trID_key) << std::endl;
  inputHandler<encT> pI(filepath_v, _k, _w, _h, &_lsh_vg, &_npositions);
  if (!_from_library) {
    float total_genome_len;
    if (_input_kmers)
      total_genome_len = pI.readInput(static_cast<uint64_t>(DEFAULT_BATCH_SIZE));
    else
      total_genome_len = pI.extractInput(1);
    uint64_t size_basis = pI.lsh_enc_vec.size();
#pragma omp critical(update_metadata)
    {
      _basis_to_size[trID_key] = size_basis;
      tT tmp_trID = trID_key;
      while (tmp_trID != 0) {
        _trID_to_size[tmp_trID] += size_basis;
        _trID_to_length[tmp_trID] += total_genome_len / _trID_to_ngenomes[tmp_trID];
        tmp_trID = _tax_record.parent_vec()[tmp_trID];
      }
      _root_size += size_basis;
    }
    bool is_ok = pI.saveInput(_library_dirpath, trID_key, _total_batches, _tbatch_size);
    if (!is_ok) {
      std::cerr << "Error saving LSH-value and encoding pairs for " << _tax_record.tiID_from_trID(trID_key) << std::endl;
      exit(EXIT_FAILURE);
    }
#pragma omp critical(print_asave)
    LOG(INFO) << COND(_log) << "LSH & encodings are saved for " << _tax_record.tiID_from_trID(trID_key) << std::endl;
  }
#pragma omp critical(emplace_iS)
  _inputStream_map.emplace(std::make_pair(trID_key, inputStream<encT>(_library_dirpath, trID_key)));
  pI.clearInput();
}

void Library::countBasis(HTs<encT> &ts, unsigned int curr_batch)
{
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < _trID_vec.size(); ++i) {
    tT trID_key = _trID_vec[i];
    std::unordered_map<encT, uint32_t> rcounts;
    std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;
    _inputStream_map.at(trID_key).loadCounts(rcounts);
    _inputStream_map.at(trID_key).loadBatch(lsh_enc_vec, curr_batch);
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
#pragma omp parallel for schedule(static), num_threads(num_threads)
  for (uint32_t rix = 0; rix < ts.num_rows; ++rix) {
    for (unsigned int j = 0; j < ts.ind_arr[rix]; ++j) {
      if (reset_scount)
        ts.scount_arr[rix * _b + j] = 0;
      if (reset_tlca)
        ts.tlca_arr[rix * _b + j] = 0;
    }
  }
}

void Library::labelLCAs(HTs<encT> &ts, unsigned int curr_batch)
{
  static double w = 2.0;
  static unsigned int r = 4;
#pragma omp parallel for schedule(dynamic), num_threads(num_threads)
  for (unsigned int i = 0; i < _trID_vec.size(); ++i) {
    tT trID_key = _trID_vec[i];
    std::unordered_map<encT, uint32_t> rcounts;
    std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;
    _inputStream_map.at(trID_key).loadCounts(rcounts);
    _inputStream_map.at(trID_key).loadBatch(lsh_enc_vec, curr_batch);
    double p_update, sc;
    for (unsigned int i = 0; i < lsh_enc_vec.size(); ++i) {
      uint32_t rix = lsh_enc_vec[i].first - (curr_batch - 1) * _tbatch_size;
      bool update_label = false;
      for (unsigned int j = 0; j < ts.ind_arr[rix]; ++j) {
        if ((ts.enc_arr[rix * _b + j] == lsh_enc_vec[i].second)) {
          if (_labels_lca == soft_lca) {
            sc = static_cast<double>(rcounts[lsh_enc_vec[i].second]) + 1.0;
            if (ts.scount_arr[rix * _b + j] <= r) {
              update_label = true;
            } else {
              p_update = 1 - pow(1 - 1 / log2(pow((ts.scount_arr[rix * _b + j] - 1) / w, 2) + 2.0), sc);
              std::bernoulli_distribution bt(p_update);
              if (bt(gen))
                update_label = true;
            }
          } else {
            update_label = true;
          }
          if (update_label)
#pragma omp critical
            ts.tlca_arr[rix * _b + j] = _tax_record.getLowestCommonAncestor(ts.tlca_arr[rix * _b + j], trID_key);
          break;
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
  uint64_t tnb_mer = _tbatch_size * _b;
  for (unsigned int i = 0; i < _total_batches; ++i) {
    curr_batch++;
    if ((_target_batch == 0) || (_target_batch == curr_batch)) {
      std::cout << "Annotating k-mners in the library for batch " << curr_batch << "..." << std::endl;
      HTs<encT> ts_r(_rootrID, _k, _h, _b, _tbatch_size, &_lsh_vg, _ranking_method, &_tax_record);
      Library::loadBatchHTs(ts_r, curr_batch);
      LOG(INFO) << COND(_log) << "The table is loaded with " << ts_r.num_kmers << "/" << tnb_mer << " k-mers" << std::endl;
      Library::resetInfo(ts_r, true, true);
      Library::countBasis(ts_r, curr_batch);
      Library::labelLCAs(ts_r, curr_batch);
      LOG(INFO) << COND(_log) << "Leaves are counted and k-mers are annotated with soft-LCA" << std::endl;
      Library::saveBatchHTs(ts_r, curr_batch);
    }
  }
}

void Library::buildTables()
{
  if (_verbose)
    std::cout << "Total number of batches is " << _total_batches << std::endl;
  unsigned int curr_batch = 0;
  uint64_t tnb_mer = _tbatch_size * _b;
  for (unsigned int i = 0; i < _total_batches; ++i) {
    curr_batch++;
    if ((_target_batch == 0) || (_target_batch == curr_batch)) {
      std::cout << "Building the library for batch " << curr_batch << "..." << std::endl;
      HTs<encT> ts_r(_rootrID, _k, _h, _b, _tbatch_size, &_lsh_vg, _ranking_method, &_tax_record);
      Library::getBatchHTs(&ts_r, curr_batch);
      LOG(INFO) << COND(_log) << "The table is built with " << ts_r.num_kmers << "/" << tnb_mer << " k-mers" << std::endl;
      Library::resetInfo(ts_r, true, true);
      Library::countBasis(ts_r, curr_batch);
      Library::labelLCAs(ts_r, curr_batch);
      LOG(INFO) << COND(_log) << "Leaves are counted and k-mers are annotated with soft-LCA" << std::endl;
      Library::saveBatchHTs(ts_r, curr_batch);
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

uint64_t Library::getConstrainedSizeKC(tT curr_trID)
{
  uint64_t constrained_size;
  uint64_t sum_size = _trID_to_size[curr_trID];
  /* uint64_t kingdom_size = _trID_to_size[_tax_record.trID_to_lsroot()[curr_trID]]; */
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
  HTd<encT> td(ts->trID, ts->k, ts->h, ts->num_rows, ts->ptr_lsh_vg, _ranking_method, &_tax_record);
  if (_fast_mode) {
    uint64_t nkmers_batch;
    uint64_t constraint_size = getConstrainedSizeKC(td.trID);
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic), shared(td)
    for (unsigned int i = 0; i < _trID_vec.size(); ++i) {
      bool is_below = false;
      tT tmp_trID = _trID_vec[i];
      while ((tmp_trID != _rootrID) || (td.trID != _rootrID)) {
        if (tmp_trID == ts->trID) {
          is_below = true;
          break;
        } else {
          tmp_trID = _tax_record.parent_vec()[tmp_trID];
        }
      }
      if (is_below || (td.trID == _rootrID)) {
        nkmers_batch = _inputStream_map.at(_trID_vec[i]).retrieveBatch(td.enc_vvec, _tbatch_size, curr_batch, true);
      }
    }
    td.initBasis(td.trID);
    td.makeUnique();
    int64_t num_rm = static_cast<int64_t>(td.num_kmers) - static_cast<int64_t>(constraint_size);
    if (num_rm > 0)
      td.shrinkHT(static_cast<uint64_t>(num_rm), _b);
  } else {
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
  }
  LOG(INFO) << COND(_log) << "Converting HTd to HTs for taxon " << _tax_record.tiID_from_trID(ts->trID) << std::endl;
  td.convertHTs(ts);
  LOG(INFO) << COND(_log) << "The hash table has been built for " << _tax_record.tiID_from_trID(ts->trID) << std::endl;
  if (_log) {
    std::cout << "Histogram for # of table columns:" << std::endl;
    for (auto kv : ts->histRowSizes())
      std::cout << "\t" << static_cast<unsigned int>(kv.first) << " : " << kv.second << std::endl;
  }
}

void Library::getBatchHTd(HTd<encT> *td, unsigned int curr_batch)
{
  uint32_t curr_tiID = _tax_record.tiID_from_trID(td->trID);
  uint64_t nkmers_batch;
  if (_verbose)
#pragma omp critical
    std::cout << "Constructing the table for " << curr_tiID << std::endl;
  if (_tax_record.isBasis(td->trID)) {
    nkmers_batch = _inputStream_map.at(td->trID).retrieveBatch(td->enc_vvec, _tbatch_size, curr_batch, false);
    td->initBasis(td->trID);
    /* td->makeUnique(true); */
#pragma omp critical
    {
      LOG(INFO) << COND(_log) << "The dynamic hash table has been constructed for the leaf " << curr_tiID << std::endl;
      LOG(INFO) << COND(_log) << "# of k-mers in batch: " << nkmers_batch << "/" << _basis_to_size[td->trID] << std::endl;
    }
  } else {
    std::set<tT> &children_set = _tax_record.child_map()[td->trID];
    size_t nchild = children_set.size();
    std::vector<tT> sibl(nchild);
    std::copy(children_set.begin(), children_set.end(), sibl.begin());
    td->childrenHT.assign(nchild,
                          HTd<encT>(0, td->k, td->h, td->num_rows, td->ptr_lsh_vg, td->ranking_method, td->tax_record));
    for (unsigned int ti = 0; ti < nchild; ++ti) {
#pragma omp critical
      LOG(INFO) << COND(_log) << "Building child " << _tax_record.tiID_from_trID(sibl[ti]) << "/" << curr_tiID << std::endl;
      td->childrenHT[ti].trID = sibl[ti];
#pragma omp task untied
      {
        getBatchHTd(&(td->childrenHT[ti]), curr_batch);
      }
    }
#pragma omp taskwait

    LOG(INFO) << COND(_log) << "Taking the union of children of " << curr_tiID << std::endl;

    if (_tax_record.depth_vec()[td->trID] <= LSR)
      omp_set_num_threads(num_threads);
    else
      omp_set_num_threads((num_threads + num_tasks - 1) / num_tasks + 1);

    for (auto &td_c : td->childrenHT) {
      td->unionRows(td_c, false);
    }
    td->updateSize();

    if (td->childrenHT.size() > 1) {
      /* if (td->trID != 0 && td->trID != 1) { */
      /*    */
      /*     LOG(INFO) << COND(_log) << curr_tiID << " has more than one child, updating LCA labels" << std::endl; */
      /*   td->updateLCA(); // computes hard-LCA labels along the way */
      /* } */

      /* if (_tax_record.depth_vec()[td->trID] <= LSR) { */
      /*   td->filterLSR(_tax_record.depth_vec(), LSR); // removes k-mers with hard-LCA above LSR rank */
      /* } */
      uint64_t constraint_size = getConstrainedSizeKC(td->trID);
      /* uint64_t constraint_size = getConstrainedSizeSC(td->num_basis); */
      int64_t num_rm = static_cast<int64_t>(td->num_kmers) - static_cast<int64_t>(constraint_size);
      if ((num_rm > 0) && _adaptive_size) {
#pragma omp critical
        LOG(INFO) << COND(_log) << curr_tiID << " has " << num_rm << " k-mers above the constraint" << std::endl;
        td->shrinkHT(static_cast<uint64_t>(num_rm), _b);
      }
    }
    td->childrenHT.clear();
#pragma omp critical
    LOG(INFO) << COND(_log) << "The dynamic hash table has been constructed for " << curr_tiID << std::endl;
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

  LOG(INFO) << COND(_log) << "Loading the library for batch " << curr_batch << std::endl;

  FILE *encf = IO::open_file((load_dirpath + "/enc_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(encf)) {
    std::puts("I/O error when reading encoding array from the library");
    is_ok = false;
  } else
    std::fread(ts.enc_arr, sizeof(encT), _tbatch_size * _b, encf);
  std::fclose(encf);

  FILE *indf = IO::open_file((load_dirpath + "/ind_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(indf)) {
    std::puts("I/O error when reading indicator array from the library");
    is_ok = false;
  } else
    std::fread(ts.ind_arr, sizeof(uint8_t), _tbatch_size, indf);
  std::fclose(indf);

  FILE *tlcaf = IO::open_file((load_dirpath + "/tlca_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(tlcaf)) {
    std::puts("I/O error when reading taxon-LCA array from the library");
    is_ok = false;
  } else
    std::fread(ts.tlca_arr, sizeof(tT), _tbatch_size * _b, tlcaf);
  std::fclose(tlcaf);

  FILE *scountf = IO::open_file((load_dirpath + "/scount_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "rb");
  if (std::ferror(scountf)) {
    std::puts("I/O error when reading species-count array from the library");
    is_ok = false;
  } else
    std::fread(ts.scount_arr, sizeof(tT), _tbatch_size * _b, scountf);
  std::fclose(scountf);

  if (is_ok)
    LOG(NOTICE) << COND(_log) << "Successfully loaded the library for batch " << curr_batch << std::endl;
  else
    LOG(ERROR) << COND(_log) << "Failed loading the library for batch " << curr_batch << std::endl;

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
      std::puts("I/O error when writing encoding array to the library");
      is_ok = false;
    }
    std::fclose(encf);
  }

  {
    FILE *indf = IO::open_file((save_dirpath + "/ind_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.ind_arr, sizeof(uint8_t), ts.num_rows, indf);
    if (std::ferror(indf)) {
      std::puts("I/O error when writing indicator array to the library");
      is_ok = false;
    }
    std::fclose(indf);
  }

  {
    FILE *tlcaf = IO::open_file((save_dirpath + "/tlca_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.tlca_arr, sizeof(tT), ts.num_rows * ts.b, tlcaf);
    if (std::ferror(tlcaf)) {
      std::puts("I/O error when writing taxon-LCA array to the library");
      is_ok = false;
    }
    std::fclose(tlcaf);
  }

  {
    FILE *scountf = IO::open_file((save_dirpath + "/scount_arr-" + std::to_string(curr_batch)).c_str(), is_ok, "wb");
    std::fwrite(ts.scount_arr, sizeof(scT), ts.num_rows * ts.b, scountf);
    if (std::ferror(scountf)) {
      std::puts("I/O error when writing species-count array to the library");
      is_ok = false;
    }
    std::fclose(scountf);
  }

  {
    if (is_ok)
      LOG(NOTICE) << COND(_log) << "Successfully saved the library for batch " << curr_batch << std::endl;
    else
      LOG(ERROR) << COND(_log) << "Failed saving the library for batch " << curr_batch << std::endl;
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
  std::vector<std::pair<tT, uint64_t>> trIDs_sizes;
  std::vector<std::pair<tT, uint32_t>> trIDs_ngenomes;
  std::vector<std::pair<tT, float>> trIDs_lengths;

  FILE *metadata_f = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "rb");
  std::fread(&_k, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_h, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_b, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_capacity_size, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_total_batches, sizeof(uint16_t), 1, metadata_f);
  std::fread(&_tbatch_size, sizeof(uint32_t), 1, metadata_f);
  std::fread(&_num_rows, sizeof(uint64_t), 1, metadata_f);
  std::fread(&_num_species, sizeof(uint32_t), 1, metadata_f);
  std::fread(&_num_nodes, sizeof(uint32_t), 1, metadata_f);
  std::fread(&_root_size, sizeof(uint64_t), 1, metadata_f);

  bases_sizes.resize(_num_species);
  trIDs_sizes.resize(_num_nodes);
  trIDs_ngenomes.resize(_num_nodes);
  trIDs_lengths.resize(_num_nodes);
  _positions.resize(_h);
  _npositions.resize(_k - _h);

  std::fread(bases_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_species, metadata_f);
  std::fread(trIDs_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_nodes, metadata_f);
  std::fread(trIDs_ngenomes.data(), sizeof(std::pair<tT, uint32_t>), _num_nodes, metadata_f);
  std::fread(trIDs_lengths.data(), sizeof(std::pair<tT, float>), _num_nodes, metadata_f);
  std::fread(_positions.data(), sizeof(uint8_t), _positions.size(), metadata_f);
  std::fread(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadata_f);

  if (std::ferror(metadata_f)) {
    std::puts("I/O error when reading metadata file from the library");
    is_ok = false;
  }
  std::fclose(metadata_f);

  for (auto kv : bases_sizes) {
    _basis_to_size[kv.first] = kv.second;
  }
  for (auto kv : trIDs_sizes) {
    _trID_to_size[kv.first] = kv.second;
  }
  for (auto kv : trIDs_ngenomes) {
    _trID_to_ngenomes[kv.first] = kv.second;
  }
  for (auto kv : trIDs_lengths) {
    _trID_to_length[kv.first] = kv.second;
  }

  if (is_ok)
    LOG(NOTICE) << COND(_log) << "Successfully loaded library metadata" << std::endl;
  else
    LOG(ERROR) << COND(_log) << "Failed loading library metadata" << std::endl;

  return is_ok;
}

bool Library::saveMetadata()
{
  bool is_ok = true;
  if (_verbose)
    std::cout << "Saving metadata of the library" << std::endl;
  std::string save_dirpath(_library_dirpath);

  std::vector<std::pair<tT, uint64_t>> bases_sizes(_basis_to_size.begin(), _basis_to_size.end());
  std::vector<std::pair<tT, uint64_t>> trIDs_sizes(_trID_to_size.begin(), _trID_to_size.end());
  std::vector<std::pair<tT, uint32_t>> trIDs_ngenomes(_trID_to_ngenomes.begin(), _trID_to_ngenomes.end());
  std::vector<std::pair<tT, float>> trIDs_lengths(_trID_to_length.begin(), _trID_to_length.end());

  assert(_basis_to_size.size() == _num_species);
  assert(_trID_to_size.size() == _num_nodes);
  assert(trIDs_ngenomes.size() == _num_nodes);
  assert(trIDs_lengths.size() == _num_nodes);

  FILE *metadata_f = IO::open_file((save_dirpath + "/metadata").c_str(), is_ok, "wb");
  std::fwrite(&_k, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_h, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_b, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_capacity_size, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_total_batches, sizeof(uint16_t), 1, metadata_f);
  std::fwrite(&_tbatch_size, sizeof(uint32_t), 1, metadata_f);
  std::fwrite(&_num_rows, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(&_num_species, sizeof(uint32_t), 1, metadata_f);
  std::fwrite(&_num_nodes, sizeof(uint32_t), 1, metadata_f);
  std::fwrite(&_root_size, sizeof(uint64_t), 1, metadata_f);
  std::fwrite(bases_sizes.data(), sizeof(std::pair<tT, uint64_t>), bases_sizes.size(), metadata_f);
  std::fwrite(trIDs_sizes.data(), sizeof(std::pair<tT, uint64_t>), trIDs_sizes.size(), metadata_f);
  std::fwrite(trIDs_ngenomes.data(), sizeof(std::pair<tT, uint32_t>), trIDs_ngenomes.size(), metadata_f);
  std::fwrite(trIDs_lengths.data(), sizeof(std::pair<tT, float>), trIDs_lengths.size(), metadata_f);
  std::fwrite(_positions.data(), sizeof(uint8_t), _positions.size(), metadata_f);
  std::fwrite(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadata_f);

  if (std::ferror(metadata_f)) {
    std::puts("I/O error when writing metadata file to the library");
    is_ok = false;
  }
  std::fclose(metadata_f);

  if (is_ok)
    LOG(NOTICE) << COND(_log) << "Successfully saved library metadata" << std::endl;
  else
    LOG(ERROR) << COND(_log) << "Failed saving library metadata" << std::endl;

  return is_ok;
}
