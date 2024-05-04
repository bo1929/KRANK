#include "query.h"
#include "common.h"

#define ROOT 1

Query::Query(std::vector<std::string> library_dirpaths,
             const char *output_dirpath,
             const char *query_filepath,
             float tvote_threshold,
             uint8_t max_match_hdist,
             bool save_match_info,
             bool verbose,
             bool log)
  : _library_dirpaths(library_dirpaths)
  , _output_dirpath(output_dirpath)
  , _query_filepath(query_filepath)
  , _tvote_threshold(tvote_threshold)
  , _max_match_hdist(max_match_hdist)
  , _save_match_info(save_match_info)
  , _log(log)
{
  _num_libraries = _library_dirpaths.size();

  for (unsigned int i = 0; i < _num_libraries; ++i) {
    _slib_ptr_v.push_back(std::make_unique<QLibrary>(_library_dirpaths[i].c_str(), _log));
    if (i == 0) {
      _k = _slib_ptr_v[i]->_k;
    } else if (!(_k == _slib_ptr_v[i]->_k)) {
      std::cerr << "All libraries to search in must have the same k-mer length" << std::endl;
      exit(EXIT_FAILURE);
    }
    bool is_compatible = true;
    for (auto &kv : _slib_ptr_v[i]->_tID_to_taxID) {
      is_compatible = is_compatible && ((_parent_inmap.count(kv.second) == 0) ||
                                        (_parent_inmap[kv.second] == _slib_ptr_v[i]->_parent_inmap[kv.second]));
      is_compatible = is_compatible && ((_depth_inmap.count(kv.second) == 0) ||
                                        (_depth_inmap[kv.second] == _slib_ptr_v[i]->_depth_inmap[kv.second]));
      is_compatible = is_compatible && ((_name_inmap.count(kv.second) == 0) ||
                                        (_name_inmap[kv.second] == _slib_ptr_v[i]->_name_inmap[kv.second]));
      is_compatible = is_compatible && ((_rank_inmap.count(kv.second) == 0) ||
                                        (_rank_inmap[kv.second] == _slib_ptr_v[i]->_rank_inmap[kv.second]));
      if (!is_compatible) { // If so, can merge.
        std::cerr << "All libraries must share the same taxonomy" << std::endl;
        exit(EXIT_FAILURE);
      }
      _parent_inmap[kv.second] = _slib_ptr_v[i]->_parent_inmap[kv.second];
      _depth_inmap[kv.second] = _slib_ptr_v[i]->_depth_inmap[kv.second];
      _name_inmap[kv.second] = _slib_ptr_v[i]->_name_inmap[kv.second];
      _rank_inmap[kv.second] = _slib_ptr_v[i]->_rank_inmap[kv.second];
    }
  }

  std::unordered_map<uint32_t, uint32_t> taxIDs_to_libcounts;
  for (unsigned int i = 0; i < _num_libraries; ++i) {
    for (auto &kv : _slib_ptr_v[i]->_tID_to_length) {
      taxIDs_to_libcounts[_slib_ptr_v[i]->_tID_to_taxID[kv.first]]++;
      _taxID_to_length[_slib_ptr_v[i]->_tID_to_taxID[kv.first]] += kv.second;
    }
  }
  for (auto &kv : _taxID_to_length) {
    _taxID_to_length[kv.first] /= taxIDs_to_libcounts[kv.first];
  }

  if (IO::ensureDirectory(_output_dirpath)) {
    std::cout << "Results will be written to " << _output_dirpath << std::endl;
  } else {
    std::cerr << "Directory to output results can not be found " << _output_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }

  if (IO::checkFASTAQ(_query_filepath)) {
    ghc::filesystem::path qpath{_query_filepath};
    qpath.replace_extension();
    std::string fpath = _query_filepath;
    std::string queryID = qpath.generic_string();
    _queryID_to_path[queryID] = fpath;
  } else {
    std::ifstream query_file(_query_filepath);
    if (!query_file.good()) {
      std::cerr << "Error opening " << _query_filepath << std::endl;
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
    query_file.close();
  }
}

void Query::postprocessProfile(std::unordered_map<uint32_t, float> &profile_corrected,
                               std::unordered_map<uint32_t, float> &profile_accumulated)
{
  std::unordered_map<std::string, float> rank_sum, rank_crsum;
  for (auto &kv : profile_accumulated) {
    std::string tmp_rank = _rank_inmap[kv.first];
    tmp_rank = tmp_rank == "no rank" ? tmp_rank + std::to_string(_depth_inmap[kv.first]) : tmp_rank;
    rank_sum[tmp_rank] += kv.second;
    rank_crsum[tmp_rank] += kv.second * _taxID_to_length[kv.first];
  }
  for (auto &kv : profile_accumulated) {
    std::string tmp_rank = _rank_inmap[kv.first];
    tmp_rank = tmp_rank == "no rank" ? tmp_rank + std::to_string(_depth_inmap[kv.first]) : tmp_rank;
    profile_corrected[kv.first] = kv.second * _taxID_to_length[kv.first] / rank_crsum[tmp_rank] * rank_sum[tmp_rank];
  }
}

void Query::profileBatch(std::unordered_map<uint32_t, float> &profile_accumulated, std::vector<tvote_info_t> &total_vinfo_v)
{
  uint32_t curr_taxID;
  std::string curr_rank, curr_name;
  for (auto &vi : total_vinfo_v) {
    curr_taxID = vi.pred_taxID;
    if ((curr_taxID != 0) && (vi.tvote_n > _tvote_threshold)) {
      while (curr_taxID != 0) {
        profile_accumulated[curr_taxID] += 1.0;
        curr_taxID = _parent_inmap[curr_taxID];
      }
    }
  }
}

void Query::classifyBatch(std::vector<tvote_info_t> &total_vinfo_v,
                          vec_str &names_vec,
                          vvec_uint32 &tlca_vec_or,
                          vvec_uint32 &tlca_vec_rc,
                          vvec_uint8 &hdist_vec_or,
                          vvec_uint8 &hdist_vec_rc)
{
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (std::size_t ix = 0; ix < names_vec.size(); ++ix) {
    std::unordered_map<uint32_t, float> tvotes_map_or, tvotes_map_rc;
    uint32_t curr_tID;
    for (std::size_t j = 0; j < tlca_vec_or[ix].size(); ++j) {
      float vote = pow((1 - static_cast<float>(hdist_vec_or[ix][j]) / static_cast<float>(_k)), _k);
      curr_tID = tlca_vec_or[ix][j];
      tvotes_map_or[curr_tID] += vote;
      while (curr_tID != 0) {
        curr_tID = _parent_inmap[curr_tID];
        tvotes_map_or[curr_tID] += vote;
      }
    }
    for (std::size_t j = 0; j < tlca_vec_rc[ix].size(); ++j) {
      float vote = pow((1 - static_cast<float>(hdist_vec_rc[ix][j]) / static_cast<float>(_k)), _k);
      curr_tID = tlca_vec_rc[ix][j];
      tvotes_map_rc[curr_tID] += vote;
      while (curr_tID != 0) {
        curr_tID = _parent_inmap[curr_tID];
        tvotes_map_rc[curr_tID] += vote;
      }
    }
    if ((tlca_vec_or[ix].size() + tlca_vec_rc[ix].size()) > 0) {
      std::unordered_map<uint32_t, float> &tvotes_map =
        (tvotes_map_rc[ROOT] > tvotes_map_or[ROOT]) ? tvotes_map_rc : tvotes_map_or;
      float majority_th = tvotes_map[ROOT] / 2;
      for (auto &kv : tvotes_map) {
        if ((kv.second > majority_th) && (_depth_inmap[kv.first] > _depth_inmap[total_vinfo_v[ix].pred_taxID])) {
          total_vinfo_v[ix].pred_taxID = kv.first;
          total_vinfo_v[ix].tvote_n = kv.second;
          total_vinfo_v[ix].tvote_r = tvotes_map[ROOT];
        }
      }
    }
  }
}

void Query::processBatch(std::vector<sseq_t> &seqBatch,
                         vec_str &names_vec,
                         vvec_uint32 &tlca_vec_or,
                         vvec_uint32 &tlca_vec_rc,
                         vvec_uint8 &hdist_vec_or,
                         vvec_uint8 &hdist_vec_rc)
{
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (std::size_t ix = 0; ix < seqBatch.size(); ++ix) {
    std::vector<uint8_t> hdist_vor;
    std::vector<uint8_t> hdist_vrc;
    std::vector<uint32_t> tlca_vor;
    std::vector<uint32_t> tlca_vrc;
    names_vec[ix] = seqBatch[ix].name;
    std::string kmer_seq;
    unsigned int kix, mi, cc;
    uint32_t rix;
    uint32_t enc32_lr, enc32_bp;
    uint64_t enc64_bp, enc64_lr, cenc64_bp, cenc64_lr;
    uint8_t min_dist, closest_di;
    if (seqBatch[ix].len > _k) {
      for (kix = mi = 0; kix < seqBatch[ix].len; ++kix) {
        cc = seq_nt4_table[static_cast<uint8_t>(seqBatch[ix].nseq[kix])];
        if (cc < 4) { // not an "N" base
          mi++;
          if (mi == _k) { // we find a k-mer
            kmer_seq = seqBatch[ix].nseq.substr(kix - (_k - 1), _k);
            kmerEncodingCompute(kmer_seq.c_str(), enc64_lr, enc64_bp);
          } else if (mi > _k) // updates
          {
            kmer_seq = seqBatch[ix].nseq[kix];
            kmerEncodingUpdate(kmer_seq.c_str(), enc64_lr, enc64_bp);
          }
          cenc64_bp = enc64_bp & _mask_bp;
          cenc64_lr = enc64_lr & _mask_lr;
          for (bool r : {false, true}) {
            if (r) {
              cenc64_bp = revcomp64bp(cenc64_bp, _k);
              cenc64_lr = conv64bp2lr(cenc64_bp, _k);
            }
            bool pm = false;
            for (unsigned int lix = 0; lix < _num_libraries && (mi >= _k); ++lix) {
              min_dist = _k;
              closest_di = 0;
              uint32_t max_rix = std::numeric_limits<uint32_t>::max();
              max_rix = max_rix >> (32 - 2 * _slib_ptr_v[lix]->_h);
              rix = computeValueLSH(cenc64_bp, _slib_ptr_v[lix]->_lsh_vg);
              assert(rix <= max_rix);
              for (uint8_t di = 0; di < _slib_ptr_v[lix]->_ind_arr[rix]; ++di) {
                uint8_t dist;
                if (std::is_same<encT, uint64_t>::value) {
                  dist = computeHammingDistance64(cenc64_lr, _slib_ptr_v[lix]->_enc_arr[rix * _slib_ptr_v[lix]->_b + di]);
                } else if (std::is_same<encT, uint32_t>::value) {
                  drop64Encoding32(_slib_ptr_v[lix]->_npositions, cenc64_bp, cenc64_lr, enc32_bp, enc32_lr);
                  dist = computeHammingDistance32(enc32_lr, _slib_ptr_v[lix]->_enc_arr[rix * _slib_ptr_v[lix]->_b + di]);
                } else {
                  std::puts("Available encoding types are 'uint64_t' and 'uint32_t'\n.");
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
                    tlca_vrc.push_back((
                      _slib_ptr_v[lix]->_tID_to_taxID)[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]]);
                  } else {
                    hdist_vor.push_back(min_dist);
                    tlca_vor.push_back((
                      _slib_ptr_v[lix]->_tID_to_taxID)[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]]);
                  }
                  pm = true;
                } else {
                  if (r) {
                    if (hdist_vrc.back() == min_dist) {
                      tlca_vrc.back() = getLowestCommonAncestor(
                        tlca_vrc.back(),
                        _slib_ptr_v[lix]
                          ->_tID_to_taxID[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]]);
                    } else if (hdist_vrc.back() > min_dist) {
                      hdist_vrc.back() = min_dist;
                      tlca_vrc.back() =
                        (_slib_ptr_v[lix]
                           ->_tID_to_taxID)[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]];
                    }
                  } else {
                    if (hdist_vor.back() == min_dist) {
                      tlca_vor.back() = getLowestCommonAncestor(
                        tlca_vor.back(),
                        _slib_ptr_v[lix]
                          ->_tID_to_taxID[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]]);
                    } else if (hdist_vor.back() > min_dist) {
                      hdist_vor.back() = min_dist;
                      tlca_vor.back() =
                        (_slib_ptr_v[lix]
                           ->_tID_to_taxID)[_slib_ptr_v[lix]->_tlca_arr[rix * _slib_ptr_v[lix]->_b + closest_di]];
                    }
                  }
                }
              }
            }
          }
        } else {
          mi = 0;
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
}

void Query::perform(uint64_t rbatch_size)
{
  std::vector<std::string> keys;
  for (auto &kv : _queryID_to_path) {
    keys.push_back(kv.first);
  }

  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  _mask_bp = u64m >> (32 - _k) * 2;
  _mask_lr = ((u64m >> (64 - _k)) << 32) + ((u64m << 32) >> (64 - _k));

  for (unsigned int i = 0; i < _queryID_to_path.size(); ++i) {
    std::string queryID = keys[i];
    std::fstream ofs_minfo, ofs_clsinfo, ofs_aprofile;
    std::string output_file(_output_dirpath);
    if (_save_match_info) {
      ofs_minfo.open(output_file + "/match_info-" + queryID, std::fstream::out);
    }
    ofs_clsinfo.open(output_file + "/classification_info-" + queryID, std::fstream::out);
    ofs_aprofile.open(output_file + "/abundance_profile-" + queryID, std::fstream::out);
    if (ofs_minfo.fail() || ofs_clsinfo.fail() || ofs_aprofile.fail())
      std::cerr << "Failed to open output file in " << output_file << std::endl;
    kseq_t *reader = IO::getReader(_queryID_to_path[queryID].c_str());
    rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
    std::vector<sseq_t> seqBatch;
    std::unordered_map<uint32_t, float> profile_accumulated;
    IO::readBatch(seqBatch, reader, rbatch_size);
    uint64_t tnum_reads = 0;
    while (!(seqBatch.empty())) {
      std::vector<tvote_info_t> total_vinfo_v(seqBatch.size());
      vec_str names_vec(seqBatch.size());
      vvec_uint8 hdist_vec_or(seqBatch.size());
      vvec_uint8 hdist_vec_rc(seqBatch.size());
      vvec_uint32 tlca_vec_or(seqBatch.size());
      vvec_uint32 tlca_vec_rc(seqBatch.size());
      Query::processBatch(seqBatch, names_vec, tlca_vec_or, tlca_vec_rc, hdist_vec_or, hdist_vec_rc);
      Query::classifyBatch(total_vinfo_v, names_vec, tlca_vec_or, tlca_vec_rc, hdist_vec_or, hdist_vec_rc);
      Query::profileBatch(profile_accumulated, total_vinfo_v);
      if (_save_match_info) {
        for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
          ofs_minfo << names_vec[ix] << std::endl;
          ofs_minfo << "or";
          for (unsigned int mi = 0; mi < tlca_vec_or[ix].size(); ++mi) {
            ofs_minfo << " " << tlca_vec_or[ix][mi] << ":" << std::to_string(hdist_vec_or[ix][mi]);
          }
          ofs_minfo << std::endl;
          ofs_minfo << "rc";
          for (unsigned int mi = 0; mi < tlca_vec_rc[ix].size(); ++mi) {
            ofs_minfo << " " << tlca_vec_rc[ix][mi] << ":" << std::to_string(hdist_vec_rc[ix][mi]);
          }
          ofs_minfo << std::endl;
        }
      }
      float th_ratio;
      std::string curr_rank, curr_name;
      ofs_clsinfo << "SEQ_ID\tRANK\tTAXON_ID\tTAXON_NAME\tPREDICTION_SCORE\tMATCH_SCORE" << std::endl;
      for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
        if ((total_vinfo_v[ix].pred_taxID != 0) && (total_vinfo_v[ix].tvote_n > _tvote_threshold)) {
          curr_rank = _rank_inmap[total_vinfo_v[ix].pred_taxID];
          curr_name = _name_inmap[total_vinfo_v[ix].pred_taxID];
          th_ratio = total_vinfo_v[ix].tvote_n / total_vinfo_v[ix].tvote_r;
          ofs_clsinfo << names_vec[ix];
          ofs_clsinfo << "\t" << curr_rank << "\t" << total_vinfo_v[ix].pred_taxID << "\t" << curr_name;
          ofs_clsinfo << "\t" << std::setprecision(3) << th_ratio << "\t" << total_vinfo_v[ix].tvote_n << std::endl;
        } else {
          ofs_clsinfo << names_vec[ix] << "\tU\tNA\tNA\tNA\t" << total_vinfo_v[ix].tvote_n << std::endl;
        }
      }
      tnum_reads += seqBatch.size();
      IO::readBatch(seqBatch, reader, rbatch_size);
    }
    kseq_destroy(reader);
    gzclose(reader->f->f);

    std::string curr_rank, curr_name;
    ofs_aprofile << "RANK\tTAXON_ID\tTAXON_NAME\tREAD_COUNT\tREAD_ABUNDANCE\tCELL_ABUNDANCE" << std::endl;
    std::vector<std::pair<uint32_t, float>> profile_acc_v(profile_accumulated.begin(), profile_accumulated.end());
    std::sort(profile_acc_v.begin(),
              profile_acc_v.end(),
              [](const std::pair<uint32_t, float> &l, const std::pair<uint32_t, float> &r) { return l.second < r.second; });
    for (auto &kv : profile_accumulated) {
      profile_accumulated[kv.first] = kv.second / tnum_reads;
    }
    std::unordered_map<uint32_t, float> profile_corrected;
    postprocessProfile(profile_corrected, profile_accumulated);
    for (auto &tc : profile_acc_v) {
      ofs_aprofile << _rank_inmap[tc.first] << "\t" << tc.first << "\t" << _name_inmap[tc.first];
      ofs_aprofile << "\t" << std::setprecision(6) << tc.second << "\t" << profile_accumulated[tc.first] << "\t"
                   << profile_corrected[tc.first] << std::endl;
    }

    if (_save_match_info) {
      ofs_minfo.close();
    }
    ofs_clsinfo.close();
    ofs_aprofile.close();
  }
}

Query::QLibrary::QLibrary(const char *library_dirpath, bool log)
  : _library_dirpath(library_dirpath)
  , _log(log)
{
  QLibrary::loadMetadata();
  _lsh_vg = generateMaskLSH(_positions);
  QLibrary::loadTaxonomy();

  if (IO::ensureDirectory(_library_dirpath)) {
    std::cout << "Library will be read from " << _library_dirpath << std::endl;
  } else {
    std::cerr << "Library can not be found at " << _library_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }

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
    /* if (_log) */
    /*   LOG(INFO) << "Allocated memory for the species-count array" << std::endl; */
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
        LOG(INFO) << "Failed to load the library into the memory: " << i << std::endl;
      exit(EXIT_FAILURE);
    } else {
      if (_log)
        LOG(INFO) << "Library has been loaded into the memory: " << i << "/" << _total_batches << std::endl;
    }
  }

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

  if (_log) {
    std::map<uint8_t, uint64_t> hist_map;
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

bool Query::QLibrary::loadMetadata()
{
  bool is_ok = true;
  if (_log)
    LOG(INFO) << "Loading metadata of the library" << std::endl;
  std::string load_dirpath(_library_dirpath);

  std::vector<std::pair<tT, uint64_t>> bases_sizes;
  std::vector<std::pair<tT, uint64_t>> tIDs_sizes;
  std::vector<std::pair<tT, uint32_t>> tIDs_ngenomes;
  std::vector<std::pair<tT, float>> tIDs_lengths;

  FILE *metadata_f = IO::open_file((load_dirpath + "/metadata").c_str(), is_ok, "rb");
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
  tIDs_sizes.resize(_num_nodes);
  tIDs_ngenomes.resize(_num_nodes);
  tIDs_lengths.resize(_num_nodes);
  _positions.resize(_h);
  _npositions.resize(_k - _h);

  std::fread(bases_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_species, metadata_f);
  std::fread(tIDs_sizes.data(), sizeof(std::pair<tT, uint64_t>), _num_nodes, metadata_f);
  std::fread(tIDs_ngenomes.data(), sizeof(std::pair<tT, uint32_t>), _num_nodes, metadata_f);
  std::fread(tIDs_lengths.data(), sizeof(std::pair<tT, float>), _num_nodes, metadata_f);
  std::fread(_positions.data(), sizeof(uint8_t), _positions.size(), metadata_f);
  std::fread(_npositions.data(), sizeof(uint8_t), _npositions.size(), metadata_f);

  if (std::ferror(metadata_f)) {
    std::puts("I/O error when reading metadata file from the library.\n");
    is_ok = false;
  }
  std::fclose(metadata_f);

  /* for (auto kv : bases_sizes) { */
  /*   _basis_to_size[kv.first] = kv.second; */
  /* } */
  /* for (auto kv : tIDs_sizes) { */
  /*   _tID_to_size[kv.first] = kv.second; */
  /* } */
  /* for (auto kv : tIDs_ngenomes) { */
  /*   _tID_to_ngenomes[kv.first] = kv.second; */
  /* } */
  for (auto kv : tIDs_lengths) {
    _tID_to_length[kv.first] = kv.second;
  }

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded library metadata" << std::endl;
    else
      LOG(ERROR) << "Failed loading library metadata" << std::endl;
  }

  return is_ok;
}

bool Query::QLibrary::loadTaxonomy()
{
  bool is_ok = true;
  if (_log)
    LOG(INFO) << "Loading taxonomy data and records of the library" << std::endl;
  std::string load_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint32_t>> tIDs_taxIDs;
  std::vector<std::pair<uint32_t, uint32_t>> taxIDs_parents;
  std::vector<std::pair<uint32_t, uint8_t>> taxIDs_depths;

  FILE *taxonomy_f = IO::open_file((load_dirpath + "/taxonomy").c_str(), is_ok, "rb");
  std::fread(&_tax_num_nodes, sizeof(tT), 1, taxonomy_f);
  std::fread(&_tax_num_input, sizeof(uint32_t), 1, taxonomy_f);

  tIDs_taxIDs.resize(_tax_num_nodes);
  _tax_parent_vec.resize(_tax_num_nodes);
  _tax_depth_vec.resize(_tax_num_nodes);
  taxIDs_parents.resize(_tax_num_nodes);
  taxIDs_depths.resize(_tax_num_nodes);

  std::fread(_tax_parent_vec.data(), sizeof(tT), _tax_num_nodes, taxonomy_f);
  std::fread(_tax_depth_vec.data(), sizeof(uint8_t), _tax_num_nodes, taxonomy_f);
  std::fread(tIDs_taxIDs.data(), sizeof(std::pair<tT, uint32_t>), _tax_num_nodes, taxonomy_f);
  std::fread(taxIDs_parents.data(), sizeof(std::pair<uint32_t, uint32_t>), _tax_num_nodes, taxonomy_f);
  std::fread(taxIDs_depths.data(), sizeof(std::pair<uint32_t, uint8_t>), _tax_num_nodes, taxonomy_f);

  for (auto kv : tIDs_taxIDs)
    _tID_to_taxID[kv.first] = kv.second;
  for (auto kv : taxIDs_parents)
    _parent_inmap[kv.first] = kv.second;
  for (auto kv : taxIDs_depths)
    _depth_inmap[kv.first] = kv.second;

  uint32_t curr_taxID;
  size_t curr_size_str;
  std::string curr_name, curr_rank;
  for (unsigned int i = 0; i < _tax_num_nodes; ++i) {
    std::fread(&curr_taxID, sizeof(uint32_t), 1, taxonomy_f);
    std::fread(&curr_size_str, sizeof(size_t), 1, taxonomy_f);
    curr_rank.resize(curr_size_str);
    std::fread(&curr_rank[0], sizeof(char), curr_size_str, taxonomy_f);
    _rank_inmap[curr_taxID] = curr_rank;
  }
  for (unsigned int i = 0; i < _tax_num_nodes; ++i) {
    std::fread(&curr_taxID, sizeof(uint32_t), 1, taxonomy_f);
    std::fread(&curr_size_str, sizeof(size_t), 1, taxonomy_f);
    curr_name.resize(curr_size_str);
    std::fread(&curr_name[0], sizeof(char), curr_size_str, taxonomy_f);
    _name_inmap[curr_taxID] = curr_name;
  }

  if (std::ferror(taxonomy_f)) {
    std::puts("I/O error when reading taxonomy-record file from the library.\n");
    is_ok = false;
  }
  std::fclose(taxonomy_f);

  if (_log) {
    if (is_ok)
      LOG(NOTICE) << "Successfully loaded taxonomy data and records" << std::endl;
    else
      LOG(ERROR) << "Failed loading taxonomy data and records" << std::endl;
  }

  return is_ok;
}
