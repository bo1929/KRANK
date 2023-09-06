#include "query.h"

Query::Query(const char* library_dirpath,
             const char* output_dirpath,
             const char* query_filepath,
             uint8_t max_match_hdist,
             bool save_match_info)
  : _library_dirpath(library_dirpath)
  , _output_dirpath(output_dirpath)
  , _query_filepath(query_filepath)
  , _max_match_hdist(max_match_hdist)
  , _save_match_info(save_match_info)
{
  Query::loadMetadata();
  _lsh_vg = generateMaskLSH(_positions);
  Query::loadTaxonomy();

  if (IO::ensureDirectory(_library_dirpath)) {
    std::cout << "Library will be read from " << _library_dirpath << std::endl;
  } else {
    std::cerr << "Library can not be found at " << _library_dirpath << std::endl;
    exit(EXIT_FAILURE);
  }
  if (IO::ensureDirectory(_output_dirpath)) {
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
    std::string queryID, path;
    if (!(std::getline(iss, queryID, '\t') && std::getline(iss, path, '\t'))) {
      std::cerr << "Failed to read file for query ID to query path map." << std::endl;
      exit(EXIT_FAILURE);
    }
    _queryID_to_path[queryID] = path;
  }

  try {
    _enc_arr = new encT[_num_rows * _b];
    std::fill(_enc_arr, _enc_arr + _num_rows * _b, 0);
    _tlca_arr = new tT[_num_rows * _b];
    std::fill(_tlca_arr, _tlca_arr + _num_rows * _b, 0);
    /* scount_arr = new scT[_num_rows * _b]; */
    /* std::fill(scount_arr, scount_arr + _num_rows * _b, 0); */
    _ind_arr = new uint8_t[_num_rows];
    std::fill(_ind_arr, _ind_arr + _num_rows, 0);
  } catch (std::bad_alloc& ba) {
    std::cerr << "Failed to allocate memory for the library." << ba.what() << std::endl;
  }
  std::cout << "Memory allocation for the library is completed successfully." << std::endl;

  std::string save_dirpath(_library_dirpath);

#pragma omp parallel for num_threads(num_threads),                                                 \
  schedule(static) shared(_enc_arr, _tlca_arr, _ind_arr)
  for (unsigned int i = 1; i <= _total_batches; ++i) {
    bool is_ok = true;
    FILE* encf =
      IO::open_file((save_dirpath + "/enc_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(encf)) {
      std::puts("I/O error when reading encoding array from the library.\n");
      is_ok = false;
    } else
      std::fread(_enc_arr + (i - 1) * _tbatch_size * _b, sizeof(encT), _tbatch_size * _b, encf);
    std::fclose(encf);
    FILE* tlcaf =
      IO::open_file((save_dirpath + "/tlca_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(tlcaf)) {
      std::puts("I/O error when reading taxon-LCA array from the library.\n");
      is_ok = false;
    } else
      std::fread(_tlca_arr + (i - 1) * _tbatch_size * _b, sizeof(tT), _tbatch_size * _b, tlcaf);
    std::fclose(tlcaf);
    FILE* indf =
      IO::open_file((save_dirpath + "/ind_arr-" + std::to_string(i)).c_str(), is_ok, "r");
    if (std::ferror(indf)) {
      std::puts("I/O error when reading indicator array from the library.\n");
      is_ok = false;
    } else
      std::fread(_ind_arr + (i - 1) * _tbatch_size, sizeof(uint8_t), _tbatch_size, indf);
    std::fclose(indf);
    if (!is_ok)
      exit(EXIT_FAILURE);
  }
  std::cout << "Library is loaded. Ready for performing queries." << std::endl;
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

void
Query::run(uint64_t rbatch_size)
{
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * _h);
  std::vector<std::string> keys;
  for (auto& kv : _queryID_to_path) {
    keys.push_back(kv.first);
  }
  for (unsigned int i = 0; i < _queryID_to_path.size(); ++i) {
    std::string queryID = keys[i];
    std::fstream ofs_match_info;
    if (_save_match_info) {
      std::string output_file(_output_dirpath);
      output_file = output_file + "/match_info-" + queryID;
      ofs_match_info.open(output_file, std::fstream::out);
      if (!ofs_match_info.is_open())
        std::cout << "Failed to open " << output_file << std::endl;
    }
    kseq_t* reader = IO::getReader(_queryID_to_path[queryID].c_str());
    rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
    std::vector<sseq_t> seqBatch = IO::readBatch(reader, rbatch_size);
    uint64_t tnum_reads = 0;
    while (!(seqBatch.empty())) {
      std::vector<std::string> names_vec(seqBatch.size());
      std::vector<std::vector<uint8_t>> hdist_vec_or(seqBatch.size());
      std::vector<std::vector<uint8_t>> hdist_vec_rc(seqBatch.size());
      std::vector<std::vector<tT>> tlca_vec_or(seqBatch.size());
      std::vector<std::vector<tT>> tlca_vec_rc(seqBatch.size());
#pragma omp parallel for num_threads(num_threads), schedule(static)
      for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
        names_vec[ix] = seqBatch[ix].name;
        std::string or_read(seqBatch[ix].nseq.begin(), seqBatch[ix].nseq.end());
        std::string rc_read(seqBatch[ix].nseq.rbegin(), seqBatch[ix].nseq.rend());
        for (bool r : { false, true }) {
          std::istringstream iss_read;
          if (r)
            iss_read.str(rc_read);
          else
            iss_read.str(or_read);
          std::string curr_cont_read;
          while (getline(iss_read, curr_cont_read, 'N')) {
            std::string kmer_str;
            encT enc_lr;
            encT enc_bp;
            uint32_t rix;
            uint8_t min_dist;
            uint8_t closest_di;
            if (curr_cont_read.length() >= _k) {
              for (unsigned int ki = 0; ki < curr_cont_read.length(); ki++) {
                min_dist = _k;
                closest_di = 0;
                if (ki == 0) {
                  kmer_str = curr_cont_read.substr(ki, _k);
                  if (r)
                    kmerEncodingComputeC(kmer_str.c_str(), enc_lr, enc_bp);
                  else
                    kmerEncodingCompute(kmer_str.c_str(), enc_lr, enc_bp);
                  ki = _k - 1;
                } else {
                  kmer_str = curr_cont_read.substr(ki, 1);
                  if (r)
                    kmerEncodingUpdateC(kmer_str.c_str(), enc_lr, enc_bp);
                  else
                    kmerEncodingUpdate(kmer_str.c_str(), enc_lr, enc_bp);
                }
                rix = computeValueLSH(enc_bp, _lsh_vg);
                assert(rix <= max_rix);
                for (uint8_t di = 0; di < _ind_arr[rix]; ++di) {
                  uint8_t dist;
                  if (std::is_same<encT, uint64_t>::value) {
                    dist = computeHammingDistance64(enc_lr, _enc_arr[rix * _b + di]);
                  } else if (std::is_same<encT, uint32_t>::value) {
                    dist = computeHammingDistance32(enc_lr, _enc_arr[rix * _b + di]);
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
                  if (r) {
                    hdist_vec_rc[ix].push_back(min_dist);
                    tlca_vec_rc[ix].push_back(_tlca_arr[rix * _b + closest_di]);
                  } else {
                    hdist_vec_or[ix].push_back(min_dist);
                    tlca_vec_or[ix].push_back(_tlca_arr[rix * _b + closest_di]);
                  }
                }
              }
            }
          }
        }
      }
      if (_save_match_info) {
        for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
          ofs_match_info << names_vec[ix] << std::endl;
          ofs_match_info << "or";
          for (unsigned int mi = 0; mi < tlca_vec_or[ix].size(); ++mi) {
            ofs_match_info << " " << _tID_to_taxID[tlca_vec_or[ix][mi]] << ":"
                           << std::to_string(hdist_vec_or[ix][mi]);
          }
          ofs_match_info << std::endl;
          ofs_match_info << "rc";
          for (unsigned int mi = 0; mi < tlca_vec_rc[ix].size(); ++mi) {
            ofs_match_info << " " << _tID_to_taxID[tlca_vec_rc[ix][mi]] << ":"
                           << std::to_string(hdist_vec_rc[ix][i]);
          }
          ofs_match_info << std::endl;
        }
      }
      tnum_reads += seqBatch.size();
      if (seqBatch.size() == rbatch_size)
        seqBatch = IO::readBatch(reader, rbatch_size);
      else
        seqBatch.clear();
    }
    if (_save_match_info)
      ofs_match_info.close();
  }
}

bool
Query::loadMetadata()
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);

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

  return is_ok;
}

bool
Query::loadTaxonomy()
{
  bool is_ok = true;
  std::string save_dirpath(_library_dirpath);
  std::vector<std::pair<tT, uint64_t>> tIDs_taxIDs;

  FILE* taxonomyf = IO::open_file((save_dirpath + "/taxonomy").c_str(), is_ok, "rb");
  std::fread(&_tax_num_input, sizeof(tT), 1, taxonomyf);
  std::fread(&_tax_num_nodes, sizeof(tT), 1, taxonomyf);
  tIDs_taxIDs.resize(_tax_num_nodes);
  _tax_parent_vec.resize(_tax_num_nodes);
  _tax_depth_vec.resize(_tax_num_nodes);
  std::fread(tIDs_taxIDs.data(), sizeof(std::pair<tT, uint64_t>), _tax_num_nodes, taxonomyf);
  std::fread(_tax_parent_vec.data(), sizeof(tT), _tax_num_nodes, taxonomyf);
  std::fread(_tax_depth_vec.data(), sizeof(uint8_t), _tax_num_nodes, taxonomyf);

  for (auto& kv : tIDs_taxIDs)
    _tID_to_taxID[kv.first] = kv.second;

  if (std::ferror(taxonomyf)) {
    std::puts("I/O error when reading taxonomy-record file from the library.\n");
    is_ok = false;
  }
  std::fclose(taxonomyf);

  return is_ok;
}