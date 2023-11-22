#include "query.h"

Query::Query(std::vector<std::string> library_dirpaths,
             const char *output_dirpath,
             const char *query_filepath,
             uint8_t max_match_hdist,
             bool save_match_info,
             bool verbose,
             bool log)
  : _library_dirpaths(library_dirpaths)
  , _output_dirpath(output_dirpath)
  , _query_filepath(query_filepath)
  , _max_match_hdist(max_match_hdist)
  , _save_match_info(save_match_info)
  , _verbose(verbose)
  , _log(log)
{
  _num_libraries = _library_dirpaths.size();
  _slib_ptr_v.resize(_num_libraries);

  for (unsigned int i = 0; i < _num_libraries; ++i) {
    _slib_ptr_v[i] = std::make_unique<SearchL>(_library_dirpaths[i].c_str(), _verbose, _log);
    if (i == 0) {
      _tID_to_taxID = _slib_ptr_v[i]->_tID_to_taxID;
      _k = _slib_ptr_v[i]->_k;
      _tax_parent_vec = _slib_ptr_v[i]->_tax_parent_vec;
      _tax_depth_vec = _slib_ptr_v[i]->_tax_depth_vec;
    } else {
      if (!map_compare(_slib_ptr_v[i]->_tID_to_taxID, _tID_to_taxID) || (_tax_depth_vec != _slib_ptr_v[i]->_tax_depth_vec) ||
          (_tax_parent_vec != _slib_ptr_v[i]->_tax_parent_vec)) {
        std::cerr << "All libraries to search in have to be built with the same taxonomy" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (!(_k == _slib_ptr_v[i]->_k)) {
        std::cerr << "All libraries to search in must have the same k-mer length" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  if (IO::ensureDirectory(_output_dirpath)) {
    if (_verbose)
      std::cout << "Results will be written to " << _output_dirpath << std::endl;
  } else {
    std::cerr << "Directory to output results can not be found " << _output_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream query_file(query_filepath);
  if (!query_file.good()) {
    std::cerr << "Error opening " << query_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(query_file, line)) {
    std::istringstream iss(line);
    std::string queryID, fpath;
    if (!(std::getline(iss, queryID, '\t') && std::getline(iss, fpath, '\t'))) {
      std::cerr << "Failed to read file for query ID to query path map" << std::endl;
      exit(EXIT_FAILURE);
    }
    _queryID_to_path[queryID] = fpath;
  }

  Query::run();
}

void Query::run(uint64_t rbatch_size)
{
  std::vector<std::string> keys;
  for (auto &kv : _queryID_to_path) {
    keys.push_back(kv.first);
  }
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  uint64_t mask_bp = u64m >> (32 - _k) * 2;
  uint64_t mask_lr = ((u64m >> (64 - _k)) << 32) + ((u64m << 32) >> (64 - _k));
  for (unsigned int i = 0; i < _queryID_to_path.size(); ++i) {
    std::string queryID = keys[i];
    std::cout << "Searching for k-mer matches for the query " << queryID << std::endl;
    std::fstream ofs_match_info;
    if (_save_match_info) {
      if (_log)
        LOG(INFO) << "Opening file to output match information" << std::endl;
      std::string output_file(_output_dirpath);
      output_file = output_file + "/match_info-" + queryID;
      ofs_match_info.open(output_file, std::fstream::out);
      if (!ofs_match_info.is_open())
        std::cerr << "Failed to open " << output_file << std::endl;
    }
    kseq_t *reader = IO::getReader(_queryID_to_path[queryID].c_str());
    rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
    std::vector<sseq_t> seqBatch;
    IO::readBatch(seqBatch, reader, rbatch_size);
    uint64_t tnum_reads = 0;
    if (_log)
      LOG(INFO) << "Batch size for the query reads in the current query file is " << seqBatch.size() << std::endl;
    while (!(seqBatch.empty())) {
      std::vector<std::string> names_vec(seqBatch.size());
      std::vector<std::vector<uint8_t>> hdist_vec_or(seqBatch.size());
      std::vector<std::vector<uint8_t>> hdist_vec_rc(seqBatch.size());
      std::vector<std::vector<tT>> tlca_vec_or(seqBatch.size());
      std::vector<std::vector<tT>> tlca_vec_rc(seqBatch.size());
#pragma omp parallel for num_threads(num_threads), schedule(static)
      for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
        std::vector<uint8_t> hdist_vor;
        std::vector<uint8_t> hdist_vrc;
        std::vector<tT> tlca_vor;
        std::vector<tT> tlca_vrc;
        names_vec[ix] = seqBatch[ix].name;
        std::string or_read(seqBatch[ix].nseq.begin(), seqBatch[ix].nseq.end());
        std::string rc_read(seqBatch[ix].nseq.rbegin(), seqBatch[ix].nseq.rend());
        for (bool r : {false, true}) {
          std::istringstream iss_read;
          if (r)
            iss_read.str(rc_read);
          else
            iss_read.str(or_read);
          std::string curr_cont_read;
          while (getline(iss_read, curr_cont_read, 'N')) {
            std::string kmer_seq;
            uint64_t enc64_bp;
            uint64_t enc64_lr;
            uint32_t enc32_lr;
            uint32_t enc32_bp;
            uint32_t rix;
            uint8_t min_dist;
            uint8_t closest_di;
            if (curr_cont_read.length() >= _k) {
              for (unsigned int kix = 0; kix < curr_cont_read.length(); kix++) {
                if (kix == 0) {
                  kmer_seq = curr_cont_read.substr(kix, _k);
                  if (r)
                    kmerEncodingComputeC(kmer_seq.c_str(), enc64_lr, enc64_bp);
                  else
                    kmerEncodingCompute(kmer_seq.c_str(), enc64_lr, enc64_bp);
                  kix = _k - 1;
                } else {
                  kmer_seq = curr_cont_read.substr(kix, 1);
                  if (r)
                    kmerEncodingUpdateC(kmer_seq.c_str(), enc64_lr, enc64_bp);
                  else
                    kmerEncodingUpdate(kmer_seq.c_str(), enc64_lr, enc64_bp);
                }
                enc64_bp = enc64_bp & mask_bp;
                enc64_lr = enc64_lr & mask_lr;
                bool pm = false;
                for (unsigned int lix = 0; lix < _num_libraries; ++lix) {
                  min_dist = _k;
                  closest_di = 0;
                  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
                  max_rix = max_rix >> (32 - 2 * _slib_ptr_v[lix]->_h);
                  rix = computeValueLSH(enc64_bp, _slib_ptr_v[lix]->_lsh_vg);
                  assert(rix <= max_rix);
                  for (uint8_t di = 0; di < _slib_ptr_v[lix]->_ind_arr[rix]; ++di) {
                    uint8_t dist;
                    if (std::is_same<encT, uint64_t>::value) {
                      dist = computeHammingDistance64(enc64_lr, _slib_ptr_v[lix]->_enc_arr[rix * _slib_ptr_v[lix]->_b + di]);
                    } else if (std::is_same<encT, uint32_t>::value) {
                      drop64Encoding32(_slib_ptr_v[lix]->_npositions, enc64_bp, enc64_lr, enc32_bp, enc32_lr);
                      dist = computeHammingDistance32(enc32_lr, _slib_ptr_v[lix]->_enc_arr[rix * _slib_ptr_v[lix]->_b + di]);
                    } else {
                      std::puts("Available encoding types are 'uint64_t' and 'uint32_t'\n!");
                      exit(EXIT_FAILURE);
                    }
                    if (dist <= min_dist) {
                      min_dist = dist;
                      closest_di = di;
                    }
                    if (dist == 0)
                      break;
                  }
                  if (min_dist <= _max_match_hdist) {
                    if (!pm) {
                      if (r) {
                        hdist_vrc.push_back(min_dist);
                        tlca_vrc.push_back(_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]);
                      } else {
                        hdist_vor.push_back(min_dist);
                        tlca_vor.push_back(_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]);
                      }
                      pm = true;
                    } else {
                      if (r) {
                        if (hdist_vrc.back() == min_dist) {
                          tlca_vrc.back() = getLowestCommonAncestor(
                            tlca_vrc.back(), _slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]);
                        } else if (hdist_vrc.back() > min_dist) {
                          hdist_vrc.back() = min_dist;
                          tlca_vrc.back() = _slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di];
                        }
                      } else {
                        if (hdist_vor.back() == min_dist) {
                          tlca_vor.back() = getLowestCommonAncestor(
                            tlca_vor.back(), _slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]);
                        } else if (hdist_vor.back() > min_dist) {
                          hdist_vor.back() = min_dist;
                          tlca_vor.back() = _slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
#pragma omp critical
        {
          hdist_vec_or[ix] = hdist_vor;
          hdist_vec_rc[ix] = hdist_vrc;
          tlca_vec_or[ix] = tlca_vor;
          tlca_vec_rc[ix] = tlca_vrc;
        }
      }
      if (_save_match_info) {
        for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
          ofs_match_info << names_vec[ix] << std::endl;
          ofs_match_info << "or";
          for (unsigned int mi = 0; mi < tlca_vec_or[ix].size(); ++mi) {
            ofs_match_info << " " << _tID_to_taxID[tlca_vec_or[ix][mi]] << ":" << std::to_string(hdist_vec_or[ix][mi]);
          }
          ofs_match_info << std::endl;
          ofs_match_info << "rc";
          for (unsigned int mi = 0; mi < tlca_vec_rc[ix].size(); ++mi) {
            ofs_match_info << " " << _tID_to_taxID[tlca_vec_rc[ix][mi]] << ":" << std::to_string(hdist_vec_rc[ix][mi]);
          }
          ofs_match_info << std::endl;
        }
      }
      tnum_reads += seqBatch.size();
      IO::readBatch(seqBatch, reader, rbatch_size);
    }
    kseq_destroy(reader);
    gzclose(reader->f->f);
    if (_save_match_info) {
      if (_log)
        LOG(INFO) << "Closing match information output file" << std::endl;
      ofs_match_info.close();
    }
  }
}

Query::SearchL::SearchL(const char *library_dirpath, bool verbose, bool log)
  : _library_dirpath(library_dirpath)
  , _verbose(verbose)
  , _log(log)
{
  if (IO::ensureDirectory(_library_dirpath)) {
    if (_verbose)
      std::cout << "Library will be read from " << _library_dirpath << std::endl;
  } else {
    std::cerr << "Library can not be found at " << _library_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }

  SearchL::loadMetadata();
  _lsh_vg = generateMaskLSH(_positions);
  SearchL::loadTaxonomy();

  try {
    _enc_arr = new encT[_num_rows * _b];
    std::fill(_enc_arr, _enc_arr + _num_rows * _b, 0);
    if (_log)
      LOG(INFO) << "Allocated memory for the encoding array" << std::endl;
    _tlca_arr = new tT[_num_rows * _b];
    std::fill(_tlca_arr, _tlca_arr + _num_rows * _b, 0);
    if (_log)
      LOG(INFO) << "Allocated memory for the taxon-LCA array" << std::endl;
    /* scount_arr = new scT[_num_rows * _b]; */
    /* std::fill(scount_arr, scount_arr + _num_rows * _b, 0); */
    _ind_arr = new uint8_t[_num_rows];
    std::fill(_ind_arr, _ind_arr + _num_rows, 0);
    if (_log)
      LOG(INFO) << "Allocated memory for the indicator array" << std::endl;
  } catch (std::bad_alloc &ba) {
    std::cerr << "Failed to allocate memory for the library" << ba.what() << std::endl;
  }
  std::cout << "Memory allocation for the library is completed successfully" << std::endl;

  std::string load_dirpath(_library_dirpath);

#pragma omp parallel for num_threads(num_threads), schedule(static) shared(_enc_arr, _tlca_arr, _ind_arr)
  for (unsigned int i = 1; i <= _total_batches; ++i) {
    bool is_ok = true;
    FILE *encf = IO::open_file((load_dirpath + "/enc_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(encf)) {
      std::puts("I/O error when reading encoding array from the library\n");
      is_ok = false;
    } else
      std::fread(_enc_arr + (i - 1) * _tbatch_size * _b, sizeof(encT), _tbatch_size * _b, encf);
    std::fclose(encf);
    FILE *tlcaf = IO::open_file((load_dirpath + "/tlca_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(tlcaf)) {
      std::puts("I/O error when reading taxon-LCA array from the library\n");
      is_ok = false;
    } else
      std::fread(_tlca_arr + (i - 1) * _tbatch_size * _b, sizeof(tT), _tbatch_size * _b, tlcaf);
    std::fclose(tlcaf);
    FILE *indf = IO::open_file((load_dirpath + "/ind_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(indf)) {
      std::puts("I/O error when reading indicator array from the library\n");
      is_ok = false;
    } else
      std::fread(_ind_arr + (i - 1) * _tbatch_size, sizeof(uint8_t), _tbatch_size, indf);
    std::fclose(indf);
    if (!is_ok) {
      if (_log)
        LOG(INFO) << "Failed to load the library into the memory" << std::endl;
      exit(EXIT_FAILURE);
    } else {
      if (_log)
        LOG(INFO) << "Library has been loaded into the memory" << std::endl;
    }
  }
  if (_verbose) {
    std::cout << "Library is loaded and ready for performing queries" << std::endl;
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
  }

  if (_log) {
    std::unordered_map<uint8_t, uint64_t> hist_map;
    for (unsigned int rix = 0; rix < _num_rows; ++rix) {
      hist_map[_ind_arr[rix]]++;
    }
    LOG(INFO) << "Percentages of # columns in the table across rows:" << std::endl;
    std::cout << "\tNumber of Columns"
              << " : "
              << "Percentage" << std::endl;
    for (auto kv : hist_map)
      std::cout << "\t" << static_cast<uint16_t>(kv.first) << " : " << static_cast<float>(kv.second) / _num_rows
                << std::endl;
  }
}

bool Query::SearchL::loadMetadata()
{
  bool is_ok = true;
  if (_verbose)
    std::cout << "Loading metadata of the library" << std::endl;
  std::string load_dirpath(_library_dirpath);

  FILE *metadataf = IO::open_file((load_dirpath + "/metadata").c_str(), is_ok, "rb");
  std::fread(&_k, sizeof(uint16_t), 1, metadataf);
  std::fread(&_h, sizeof(uint16_t), 1, metadataf);
  std::fread(&_b, sizeof(uint16_t), 1, metadataf);
  std::fread(&_capacity_size, sizeof(uint64_t), 1, metadataf);
  std::fread(&_total_batches, sizeof(uint16_t), 1, metadataf);
  std::fread(&_tbatch_size, sizeof(uint32_t), 1, metadataf);
  std::fread(&_num_rows, sizeof(uint64_t), 1, metadataf);
  std::fread(&_num_species, sizeof(uint64_t), 1, metadataf);
  std::fread(&_root_size, sizeof(uint64_t), 1, metadataf);
  _bases_sizes.resize(_num_species);
  std::fread(_bases_sizes.data(), sizeof(std::pair<tT, uint16_t>), _num_species, metadataf);
  _positions.resize(_h);
  std::fread(_positions.data(), sizeof(uint8_t), _positions.size(), metadataf);
  _npositions.resize(_k - _h);
  std::fread(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadataf);

  if (std::ferror(metadataf)) {
    std::puts("I/O error when reading metadata file from the library.\n");
    is_ok = false;
  }
  std::fclose(metadataf);

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded library metadata" << std::endl;
    else
      LOG(ERROR) << "Failed loading library metadata" << std::endl;
  }

  return is_ok;
}

template<typename T>
T Query::getLowestCommonAncestor(T a, T b)
{
  if (!a || !b) // LCA(x,0) = LCA(0,x) = x
    return a ? a : b;
  while ((a != b) && (a != 1) && (b != 1)) {
    if (_tax_depth_vec[a] < _tax_depth_vec[b]) {
      b = _tax_parent_vec[b];
      assert(b != 0);
    } else {
      a = _tax_parent_vec[a];
      assert(a != 0);
    }
  }
  if ((a == 1) || (b == 1))
    return a || b;
  return a;
}

bool Query::SearchL::loadTaxonomy()
{
  bool is_ok = true;
  if (_verbose)
    std::cout << "Loading taxonomy data and records of the library" << std::endl;
  std::string load_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint64_t>> tIDs_taxIDs;

  FILE *taxonomyf = IO::open_file((load_dirpath + "/taxonomy").c_str(), is_ok, "rb");
  std::fread(&_tax_num_input, sizeof(tT), 1, taxonomyf);
  std::fread(&_tax_num_nodes, sizeof(tT), 1, taxonomyf);
  tIDs_taxIDs.resize(_tax_num_nodes);
  _tax_parent_vec.resize(_tax_num_nodes);
  _tax_depth_vec.resize(_tax_num_nodes);
  std::fread(tIDs_taxIDs.data(), sizeof(std::pair<tT, uint64_t>), _tax_num_nodes, taxonomyf);
  std::fread(_tax_parent_vec.data(), sizeof(tT), _tax_num_nodes, taxonomyf);
  std::fread(_tax_depth_vec.data(), sizeof(uint8_t), _tax_num_nodes, taxonomyf);

  for (auto &kv : tIDs_taxIDs)
    _tID_to_taxID[kv.first] = kv.second;

  if (std::ferror(taxonomyf)) {
    std::puts("I/O error when reading taxonomy-record file from the library.\n");
    is_ok = false;
  }
  std::fclose(taxonomyf);

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded taxonomy data and records" << std::endl;
    else
      LOG(ERROR) << "Failed loading taxonomy data and records" << std::endl;
  }

  return is_ok;
}
