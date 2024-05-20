#include "io.h"

#define BUFF_SIZE 1024 * 1024

struct vhash
{
  inline std::size_t operator()(const std::pair<int, int> &v) const { return v.second; }
};

unsigned int gp_hash(const std::string &str)
{
  unsigned int b = 378551;
  unsigned int a = 63689;
  unsigned int hash = 0;

  for (std::size_t i = 0; i < str.length(); i++) {
    hash = hash * a + str[i];
    a = a * b;
  }
  return (hash & 0x7FFFFFFF);
}

template<typename T>
inline void sortColumns(vvec<T> &table)
{
  uint32_t num_rows = table.size();
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!table[rix].empty())
      std::sort(table[rix].begin(), table[rix].end());
  }
}

template<typename encT>
bool inputHandler<encT>::saveInput(const char *dirpath, tT trID_key, uint16_t total_batches, uint32_t tbatch_size)
{
  bool is_ok = true;
  auto vec_begin = lsh_enc_vec.begin();
  auto vec_end = lsh_enc_vec.end();
  for (int i = 1; i <= total_batches; ++i) {
    std::string batch_dirpath = dirpath;
    batch_dirpath += +"/batch" + std::to_string(i);
    std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
    FILE *vec_f = IO::open_file(disk_path.c_str(), is_ok, "wb");
    auto vec_p =
      std::upper_bound(vec_begin, vec_end, i * tbatch_size - 1, [](uint32_t value, const std::pair<uint32_t, encT> &p) {
        return p.first > value;
      });
    if (vec_begin != lsh_enc_vec.end()) {
      uint64_t num_elements = std::distance(vec_begin, vec_p);
      std::fwrite(&(*vec_begin), sizeof(std::pair<uint32_t, encT>), num_elements, vec_f);
      if (std::ferror(vec_f)) {
        std::puts("I/O error when writing LSH-value and encoding pairs in the below function.");
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cerr << "Error: " << strerror(errno);
        std::fclose(vec_f);
        is_ok = false;
        break;
      }
    }
    std::fclose(vec_f);
    vec_begin = vec_p;
  }

  std::string rcounts_dirpath = dirpath;
  rcounts_dirpath += "/rcounts";
  std::string rcounts_fpath = (rcounts_dirpath + "/" + std::to_string(trID_key));
  if (!ghc::filesystem::exists(rcounts_fpath)) {
    std::vector<std::pair<encT, uint64_t>> rcounts_vec(rcounts.begin(), rcounts.end());
    FILE *rcounts_f = IO::open_file(rcounts_fpath.c_str(), is_ok, "wb");
    std::fwrite(rcounts_vec.data(), sizeof(std::pair<encT, uint64_t>), rcounts_vec.size(), rcounts_f);
    if (std::ferror(rcounts_f)) {
      std::puts("I/O error when writing the genome counts for the shared k-mers.\n");
      is_ok = false;
    }
    std::fclose(rcounts_f);
  }

  return is_ok;
}

template<typename encT>
bool inputHandler<encT>::checkInput(const char *dirpath, tT trID_key, uint16_t total_batches)
{
  bool is_ready = true;
  for (int i = 1; i <= total_batches; ++i) {
    std::string batch_dirpath = dirpath;
    batch_dirpath += +"/batch" + std::to_string(i);
    std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
    if (!exists_test(disk_path.c_str())) {
      is_ready = false;
      break;
    }
  }
  return is_ready;
}

template<typename encT>
bool inputHandler<encT>::loadInput(const char *dirpath, tT trID_key, uint16_t total_batches)
{
  bool is_ok = true;
  lsh_enc_vec.clear();
  for (int i = 1; i <= total_batches; ++i) {
    std::string batch_dirpath = dirpath;
    batch_dirpath += +"/batch" + std::to_string(i);
    std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
    std::ifstream vec_ifs = IO::open_ifstream(disk_path.c_str(), is_ok);
    char buf[BUFF_SIZE];
    vec_ifs.rdbuf()->pubsetbuf(buf, BUFF_SIZE);
    size_t sr = ghc::filesystem::file_size(disk_path) / sizeof(std::pair<uint32_t, encT>);
    if (vec_ifs.good() && sr > 0) {
      size_t lsize = lsh_enc_vec.size();
      lsh_enc_vec.resize(lsize + sr);
      vec_ifs.read(reinterpret_cast<char *>(&(*(lsh_enc_vec.begin() + lsize))), sr * sizeof(std::pair<uint32_t, encT>));
    }
    if ((sr > 0) && (vec_ifs.fail() || (vec_ifs.peek() != EOF))) {
      std::puts("I/O error when reading LSH-value and encoding pairs in the below function.");
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cerr << "Error: " << strerror(errno);
      break;
    } else if (vec_ifs.eof()) {
      is_ok = true;
    }
    vec_ifs.close();
  }
  return is_ok;
}

template<typename encT>
void inputHandler<encT>::resetInput()
{
  l_rix = 0;
  curr_vix = 0;
}

template<typename encT>
void inputHandler<encT>::clearInput()
{
  lsh_enc_vec.clear();
  rcounts.clear();
  l_rix = 0;
  curr_vix = 0;
}

template<typename encT>
void inputStream<encT>::loadBatch(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec, unsigned int curr_batch)
{
  bool is_ok = true;
  std::string batch_dirpath = dirpath;
  batch_dirpath += +"/batch" + std::to_string(curr_batch);
  std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
  std::ifstream batch_ifs = IO::open_ifstream(disk_path.c_str(), is_ok);
  if (!is_ok)
    exit(EXIT_FAILURE);
  char buf[BUFF_SIZE];
  batch_ifs.rdbuf()->pubsetbuf(buf, BUFF_SIZE);
  if (batch_ifs.good()) {
    lsh_enc_vec.resize(ghc::filesystem::file_size(disk_path) / sizeof(std::pair<uint32_t, encT>));
    batch_ifs.read(reinterpret_cast<char *>(lsh_enc_vec.data()), lsh_enc_vec.size() * sizeof(std::pair<uint32_t, encT>));
  }
  if (batch_ifs.fail() || (batch_ifs.peek() != EOF)) {
    std::puts("I/O eror when reading LSH-value and encoding pairs.\n");
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "Error: " << strerror(errno);
  }
  batch_ifs.close();
}

template<typename encT>
void inputStream<encT>::loadCounts(std::unordered_map<encT, uint32_t> &rcounts)
{
  bool is_ok = true;
  std::string rcounts_path = dirpath;
  rcounts_path += "/rcounts/" + std::to_string(trID_key);
  std::vector<std::pair<encT, uint32_t>> rcounts_vec;
  uint64_t num_kmers = ghc::filesystem::file_size(rcounts_path) / sizeof(std::pair<encT, uint32_t>);
  if (num_kmers > 0) {
    FILE *rcounts_f = IO::open_file((rcounts_path).c_str(), is_ok, "rb");
    rcounts_vec.resize(num_kmers);
    std::fread(rcounts_vec.data(), sizeof(std::pair<encT, uint32_t>), num_kmers, rcounts_f);
    if (std::ferror(rcounts_f)) {
      std::puts("I/O error when reading the genome counts for the shared k-mers.\n");
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cerr << "Error: " << strerror(errno);
      is_ok = false;
    }
    std::fclose(rcounts_f);
  }
  std::copy(rcounts_vec.begin(), rcounts_vec.end(), std::inserter(rcounts, rcounts.begin()));
}

template<typename encT>
uint64_t inputStream<encT>::retrieveBatch(vvec<encT> &td, uint32_t tbatch_size, unsigned int curr_batch, bool shared_table)
{
  uint32_t bix = (curr_batch - 1) * tbatch_size;
  uint32_t eix = curr_batch * tbatch_size;
  uint32_t toff_rix = (curr_batch - 1) * tbatch_size;
  assert(td.size() >= (eix - bix));
  bool is_ok = true;
  std::string batch_dirpath = dirpath;
  batch_dirpath += +"/batch" + std::to_string(curr_batch);
  std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
  uint64_t num_kmers = ghc::filesystem::file_size(disk_path) / sizeof(std::pair<uint32_t, encT>);
  std::ifstream batch_ifs = IO::open_ifstream(disk_path.c_str(), is_ok);
  if (!is_ok)
    exit(EXIT_FAILURE);
  char buf[BUFF_SIZE];
  batch_ifs.rdbuf()->pubsetbuf(buf, BUFF_SIZE);
  uint64_t num_retrieved = 0;
  while (!batch_ifs.eof() && batch_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    batch_ifs.read((char *)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if ((lsh_enc.first >= bix) && (lsh_enc.first < eix) && !batch_ifs.eof()) {
      if (shared_table)
#pragma omp critical
        td[lsh_enc.first - toff_rix].push_back(lsh_enc.second);
      else
        td[lsh_enc.first - toff_rix].push_back(lsh_enc.second);
      num_retrieved++;
    }
    if (lsh_enc.first >= eix)
      break;
  }
  if (batch_ifs.fail() && !batch_ifs.eof()) {
    std::puts("I/O eror when reading LSH-value and encoding pairs.\n");
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "Error: " << strerror(errno);
  }
  assert(num_retrieved == num_kmers);
  batch_ifs.close();
  return num_retrieved;
}

template<typename encT>
void inputStream<encT>::removeBatch(unsigned int curr_batch)
{
  std::string batch_dirpath = dirpath;
  batch_dirpath += +"/batch" + std::to_string(curr_batch);
  std::string disk_path = batch_dirpath + "/lsh_enc_vec-" + std::to_string(trID_key);
  ghc::filesystem::remove(disk_path);
}

template<typename encT>
void inputHandler<encT>::extractMers(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec_f,
                                     std::string filepath,
                                     uint64_t rbatch_size)
{
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * h);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  uint64_t mask_bp = u64m >> (32 - k) * 2;
  uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  auto url_regexp = std::regex(
    R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
  unsigned int i, l, c;
  bool is_url = std::regex_match(filepath, url_regexp);
  if (is_url) {
#pragma omp critical(pd)
    std::cout << "Downloading: " << filepath << std::endl;
    filepath = IO::downloadURL(filepath);
  } else {
#pragma omp critical(pf)
    std::cout << "Reading: " << filepath << std::endl;
  }
  kseq_t *reader = IO::getReader(filepath.c_str());
  std::vector<sseq_t> seqBatch;
  IO::readBatch(seqBatch, reader, rbatch_size);
  while (!(seqBatch.empty())) {
    for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
      if (seqBatch[ix].len >= k) {
        uint64_t enc64_bp, enc64_lr, cenc64_bp, cenc64_lr;
        uint32_t cenc32_lr, cenc32_bp;
        uint32_t rix;
        uint8_t ldiff = w - k + 1;
        std::string kmer_seq;
        uint8_t kix = 0;
        size_t wix = lsh_enc_vec_f.size();
        std::vector<std::pair<uint32_t, encT>> lsh_enc_win(ldiff);
        lsh_enc_vec_f.resize(wix + seqBatch[ix].len - k + 1);
        for (i = l = 0; i < seqBatch[ix].len; ++i) {
          c = seq_nt4_table[static_cast<uint8_t>(seqBatch[ix].nseq[i])];
          if (c < 4) { // not an "N" base
            l++;
            if (l == k) { // we find a k-mer
              kmer_seq = seqBatch[ix].nseq.substr(i - (k - 1), k);
              kmerEncodingCompute(kmer_seq.c_str(), enc64_lr, enc64_bp);
            } else if (l > k) // updates
            {
              kmer_seq = seqBatch[ix].nseq[i];
              kmerEncodingUpdate(kmer_seq.c_str(), enc64_lr, enc64_bp);
            }
            if (l >= k) {
              cenc64_bp = enc64_bp & mask_bp;
              cenc64_lr = enc64_lr & mask_lr;
#ifdef CANONICAL
              if (cenc64_bp < revcomp64bp(cenc64_bp, k)) {
                cenc64_bp = revcomp64bp(cenc64_bp, k);
                cenc64_lr = conv64bp2lr(cenc64_bp, k);
              }
#endif
              rix = computeValueLSH(cenc64_bp, *(ptr_lsh_vg));
              assert(rix <= max_rix);
              if (std::is_same<encT, uint64_t>::value) {
                if (ldiff > 1) {
                  lsh_enc_win[kix % ldiff] = std::make_pair(rix, cenc64_lr);
                  kix++;
                } else {
                  lsh_enc_vec_f[wix] = std::make_pair(rix, cenc64_lr);
                  wix++;
                }
              } else if (std::is_same<encT, uint32_t>::value) {
                drop64Encoding32(*ptr_npositions, cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
                if (ldiff > 1) {
                  lsh_enc_win[kix % ldiff] = std::make_pair(rix, cenc32_lr);
                  kix++;
                } else {
                  lsh_enc_vec_f[wix] = std::make_pair(rix, cenc32_lr);
                  wix++;
                }
              } else {
                std::puts("Available encoding types are 'uint64_t' and 'uint32_t'.\n");
                exit(EXIT_FAILURE);
              }
            }
            if ((l >= w || ((i == seqBatch[ix].len - 1) && l >= k)) && ldiff > 1) {
              lsh_enc_vec_f[wix] = *std::min_element(
                lsh_enc_win.begin(), lsh_enc_win.end(), [](std::pair<uint32_t, encT> lhs, std::pair<uint32_t, encT> rhs) {
                  return murmur64(lhs.second) < murmur64(rhs.second);
                  /* return lhs.second < rhs.second; */
                });
              wix++;
            }
          } else
            l = 0;
        }
        lsh_enc_vec_f.resize(wix);
      }
    }
    IO::readBatch(seqBatch, reader, rbatch_size);
  }
  kseq_destroy(reader);
  gzclose(reader->f->f);
  lsh_enc_vec_f.shrink_to_fit();
  std::sort(
    lsh_enc_vec_f.begin(), lsh_enc_vec_f.end(), [](const std::pair<uint32_t, encT> &l, const std::pair<uint32_t, encT> &r) {
      return (l.first == r.first) ? l.second < r.second : l.first < r.first;
    });
  auto rone_it = std::unique(lsh_enc_vec_f.begin(), lsh_enc_vec_f.end());
  lsh_enc_vec_f.erase(rone_it, lsh_enc_vec_f.end());
  if (is_url)
    std::remove(filepath.c_str());
}

template<typename encT>
float inputHandler<encT>::extractInput(uint64_t rbatch_size)
{
  size_t last_ix = lsh_enc_vec.size();
  float total_genome_len = 0.0;
  unsigned int inc = std::min(static_cast<unsigned int>(filepath_v.size()), static_cast<unsigned int>(GENOME_BATCH_SIZE));
  for (unsigned int fix = 0; fix < filepath_v.size(); fix += inc) {
    std::vector<std::vector<std::pair<uint32_t, encT>>> lsh_enc_vec_fv(inc);
#pragma omp parallel
    {
#pragma omp single
      {
        for (unsigned int bix = 0; bix < inc; ++bix) {
#pragma omp task untied
          {
            extractMers(lsh_enc_vec_fv[bix], filepath_v[fix + bix], rbatch_size);
          }
        }
#pragma omp taskwait
        last_ix = lsh_enc_vec_fv[0].size();
        for (unsigned int bix = 0; bix < inc; ++bix) {
          lsh_enc_vec.insert(lsh_enc_vec.end(), lsh_enc_vec_fv[bix].begin(), lsh_enc_vec_fv[bix].end());
          total_genome_len += static_cast<float>(lsh_enc_vec_fv[bix].size());
        }
      }
    }
    std::sort(lsh_enc_vec.begin() + last_ix,
              lsh_enc_vec.end(),
              [](const std::pair<uint32_t, encT> &l, const std::pair<uint32_t, encT> &r) {
                return (l.first == r.first) ? l.second < r.second : l.first < r.first;
              });
    std::inplace_merge(lsh_enc_vec.begin(), lsh_enc_vec.begin() + last_ix, lsh_enc_vec.end());
    auto nuniq_it = std::unique(lsh_enc_vec.begin(), lsh_enc_vec.end());
    for (auto it = nuniq_it; it != lsh_enc_vec.end(); ++it) {
      rcounts[it->second]++;
    }
    lsh_enc_vec.erase(nuniq_it, lsh_enc_vec.end());
    last_ix = lsh_enc_vec.size();
  }
  tnum_kmers = lsh_enc_vec.size();
  return total_genome_len;
}

template<typename encT>
void inputHandler<encT>::readMers(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec_f,
                                  std::string filepath,
                                  uint64_t rbatch_size)
{
  uint64_t tnum_kmers_sum = 0;
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * h);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  uint64_t mask_bp = u64m >> (32 - k) * 2;
  uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
  kseq_t *reader = IO::getReader(filepath.c_str());
  std::vector<sseq_t> seqBatch;
  IO::readBatch(seqBatch, reader, rbatch_size);
  while (!(seqBatch.empty())) {
    uint8_t ldiff = w - k + 1;
    uint8_t kix = 0;
    std::vector<std::pair<uint32_t, encT>> lsh_enc_win(ldiff);
    uint32_t wix = lsh_enc_vec_f.size();
    lsh_enc_vec_f.resize(wix + seqBatch.size());
    for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
      uint64_t enc64_bp;
      uint64_t enc64_lr;
      uint64_t cenc64_bp;
      uint64_t cenc64_lr;
      uint32_t cenc32_lr;
      uint32_t cenc32_bp;
      uint32_t rix;
      if (seqBatch[ix].len != w) {
        std::puts("An input reference k-mer length conflicts with given k& w parameters.\n");
        exit(EXIT_FAILURE);
      }
      for (unsigned int kix = 0; kix < ldiff; ++kix) {
        if (kix == 0)
          kmerEncodingCompute(seqBatch[ix].nseq.substr(kix, k).c_str(), enc64_lr, enc64_bp);
        else
          kmerEncodingUpdate(seqBatch[ix].nseq.substr(kix + k - 1, 1).c_str(), enc64_lr, enc64_bp);
        cenc64_bp = enc64_bp & mask_bp;
        cenc64_lr = enc64_lr & mask_lr;
#ifdef CANONICAL
        if (cenc64_bp < revcomp64bp(cenc64_bp, k)) {
          cenc64_bp = revcomp64bp(cenc64_bp, k);
          cenc64_lr = conv64bp2lr(cenc64_bp, k);
        }
#endif
        rix = computeValueLSH(cenc64_bp, *(ptr_lsh_vg));
        assert(rix <= max_rix);
        if (std::is_same<encT, uint64_t>::value) {
          if (ldiff > 1)
            lsh_enc_win[kix] = std::make_pair(rix, cenc64_lr);
          else
            lsh_enc_vec_f[wix + ix] = std::make_pair(rix, cenc64_lr);
        } else if (std::is_same<encT, uint32_t>::value) {
          drop64Encoding32(*ptr_npositions, cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
          if (ldiff > 1)
            lsh_enc_win[kix] = std::make_pair(rix, cenc32_lr);
          else
            lsh_enc_vec_f[wix + ix] = std::make_pair(rix, cenc32_lr);
        } else {
          std::puts("Available encoding types are 'uint64_t' and 'uint32_t'.\n");
          exit(EXIT_FAILURE);
        }
      }
      if (ldiff > 1)
        lsh_enc_vec_f[wix + ix] = *std::min_element(
          lsh_enc_win.begin(), lsh_enc_win.end(), [](std::pair<uint32_t, encT> lhs, std::pair<uint32_t, encT> rhs) {
            return murmur64(lhs.second) < murmur64(rhs.second);
            /* return lhs.second < rhs.second; */
          });
    }
    tnum_kmers_sum += seqBatch.size();
    IO::readBatch(seqBatch, reader, rbatch_size);
  }
  kseq_destroy(reader);
  gzclose(reader->f->f);
  lsh_enc_vec_f.shrink_to_fit();
  std::sort(
    lsh_enc_vec_f.begin(), lsh_enc_vec_f.end(), [](const std::pair<uint32_t, encT> &l, const std::pair<uint32_t, encT> &r) {
      return (l.first == r.first) ? l.second < r.second : l.first < r.first;
    });
  auto rone_it = std::unique(lsh_enc_vec_f.begin(), lsh_enc_vec_f.end());
  lsh_enc_vec_f.erase(rone_it, lsh_enc_vec_f.end());
}

template<typename encT>
float inputHandler<encT>::readInput(uint64_t rbatch_size)
{
  size_t last_ix = lsh_enc_vec.size();
  float total_genome_len = 0.0;
  unsigned int inc = std::min(static_cast<unsigned int>(filepath_v.size()), static_cast<unsigned int>(GENOME_BATCH_SIZE));
  for (unsigned int fix = 0; fix < filepath_v.size(); fix += inc) {
    std::vector<std::vector<std::pair<uint32_t, encT>>> lsh_enc_vec_fv(inc);
#pragma omp parallel
    {
#pragma omp single
      {
        for (unsigned int bix = 0; bix < inc; ++bix) {
#pragma omp task untied
          {
            readMers(lsh_enc_vec_fv[bix], filepath_v[fix + bix], rbatch_size);
          }
        }
#pragma omp taskwait
        last_ix = lsh_enc_vec_fv[0].size();
        for (unsigned int bix = 0; bix < inc; ++bix) {
          lsh_enc_vec.insert(lsh_enc_vec.end(), lsh_enc_vec_fv[bix].begin(), lsh_enc_vec_fv[bix].end());
          total_genome_len += static_cast<float>(lsh_enc_vec_fv[bix].size());
        }
      }
    }
    std::sort(lsh_enc_vec.begin() + last_ix,
              lsh_enc_vec.end(),
              [](const std::pair<uint32_t, encT> &l, const std::pair<uint32_t, encT> &r) {
                return (l.first == r.first) ? l.second < r.second : l.first < r.first;
              });
    std::inplace_merge(lsh_enc_vec.begin(), lsh_enc_vec.begin() + last_ix, lsh_enc_vec.end());
    auto nuniq_it = std::unique(lsh_enc_vec.begin(), lsh_enc_vec.end());
    for (auto it = nuniq_it; it != lsh_enc_vec.end(); ++it) {
      rcounts[it->second]++;
    }
    lsh_enc_vec.erase(nuniq_it, lsh_enc_vec.end());
    last_ix = lsh_enc_vec.size();
  }
  tnum_kmers = lsh_enc_vec.size();
  return total_genome_len;
}

template<typename encT>
std::map<uint8_t, uint64_t> inputHandler<encT>::histRowSizes()
{
  std::map<uint32_t, uint8_t> row_sizes;
  for (auto kv : lsh_enc_vec) {
    row_sizes[kv.first]++;
  }
  std::map<uint8_t, uint64_t> hist_map;
  for (auto kv : row_sizes) {
    hist_map[kv.second]++;
  }
  return hist_map;
}

bool IO::ensureDirectory(const char *dirpath)
{
  bool is_ok;
  struct stat info;
  if (stat(dirpath, &info) != 0) {
    std::cerr << "Cannot access to given directory path " << dirpath << std::endl;
    is_ok = false;
  } else if (info.st_mode & S_IFDIR) {
    is_ok = true;
  } else {
    std::cerr << "Given path is not a directory " << dirpath << std::endl;
    is_ok = false;
  }
  return is_ok;
}

size_t IO::writeData(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t written = fwrite(ptr, size, nmemb, stream);
  return written;
}

std::string IO::downloadURL(std::string url)
{
  char tmp_input_path[FILENAME_MAX] = "/tmp/seq";
  const char *sx = std::to_string(gp_hash(url)).c_str();
  strcat(tmp_input_path, sx);
  strcat(tmp_input_path, ".XXXXXX");
  int tmp_fd = mkstemp(tmp_input_path);
  CURL *curl;
  FILE *fp;
  CURLcode resb;
  curl = curl_easy_init();
  if (curl) {
    fp = fopen(tmp_input_path, "wb");
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, IO::writeData);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    resb = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    fclose(fp);
  }
  return tmp_input_path;
}

kseq_t *IO::getReader(const char *fpath)
{
  gzFile fp;
  fp = gzopen(fpath, "rb");
  if (fp == nullptr) {
    std::cerr << "Failed to open file at " << fpath << std::endl;
    exit(1);
  }

  kseq_t *kseq;
  kseq = kseq_init(fp);

  return kseq;
}

uint64_t IO::adjustBatchSize(uint64_t batch_size, uint8_t num_threads) { return (batch_size / num_threads) * num_threads; }

bool IO::checkFASTAQ(const char *filepath)
{
  gzFile fp;
  fp = gzopen(filepath, "rb");
  bool is_fastaq = false;
  if (fp == nullptr) {
    std::cerr << "Failed to open file at " << filepath << std::endl;
    exit(1);
  }
  kseq_t *reader;
  reader = kseq_init(fp);
  if (kseq_read(reader) >= 0)
    is_fastaq = true;
  else
    is_fastaq = false;
  kseq_destroy(reader);
  gzclose(reader->f->f);
  return is_fastaq;
}

void IO::readBatch(std::vector<sseq_t> &seqRead, kseq_t *kseq, uint64_t batch_size)
{
  int l;
  unsigned int i = 0;
  seqRead.clear();
  seqRead.resize(batch_size);
  while ((i < batch_size) && (l = kseq_read(kseq)) >= 0) {
    sseq_t tmp_seq;
    tmp_seq.nseq.assign(kseq->seq.s);
    tmp_seq.name.assign(kseq->name.s);
    tmp_seq.len = (kseq->seq.l);
    seqRead[i] = tmp_seq;
    i++;
  }
  seqRead.resize(i);
}

FILE *IO::open_file(const char *filepath, bool &is_ok, const char *mode)
{
  FILE *f;
  f = std::fopen(filepath, mode);
  if (!f) {
    std::cerr << "File opening failed! " << filepath << std::endl;
    is_ok = false;
  }
  return f;
}

std::ifstream IO::open_ifstream(const char *filepath, bool is_ok)
{
  std::ifstream ifs;
  ifs.open(filepath, std::ios::binary | std::ios::in);
  if (!ifs.is_open()) {
    std::fprintf(stderr, "File opening failed! %s", filepath);
    is_ok = false;
  } else if (!ifs.good()) {
    std::fprintf(stderr, "Input stream is open but not good! %s", filepath);
    is_ok = false;
  }
  return ifs;
}

template float inputHandler<uint32_t>::readInput(uint64_t rbatch_size);

template float inputHandler<uint64_t>::readInput(uint64_t rbatch_size);

template float inputHandler<uint32_t>::extractInput(uint64_t rbatch_size);

template float inputHandler<uint64_t>::extractInput(uint64_t rbatch_size);

template void inputHandler<uint32_t>::extractMers(std::vector<std::pair<uint32_t, uint32_t>> &lsh_enc_vec_f,
                                                  std::string filepath,
                                                  uint64_t rbatch_size);

template void inputHandler<uint64_t>::extractMers(std::vector<std::pair<uint32_t, uint64_t>> &lsh_enc_vec_f,
                                                  std::string filepath,
                                                  uint64_t rbatch_size);

template void inputHandler<uint32_t>::readMers(std::vector<std::pair<uint32_t, uint32_t>> &lsh_enc_vec_f,
                                               std::string filepath,
                                               uint64_t rbatch_size);

template void inputHandler<uint64_t>::readMers(std::vector<std::pair<uint32_t, uint64_t>> &lsh_enc_vec_f,
                                               std::string filepath,
                                               uint64_t rbatch_size);

template bool
inputHandler<uint64_t>::saveInput(const char *dirpath, tT trID_key, uint16_t total_batches, uint32_t tbatch_size);

template bool
inputHandler<uint32_t>::saveInput(const char *dirpath, tT trID_key, uint16_t total_batches, uint32_t tbatch_size);

template bool inputHandler<uint64_t>::loadInput(const char *dirpath, tT trID_key, uint16_t total_batches);

template bool inputHandler<uint32_t>::loadInput(const char *dirpath, tT trID_key, uint16_t total_batches);

template bool inputHandler<uint64_t>::checkInput(const char *dirpath, tT trID_key, uint16_t total_batches);

template bool inputHandler<uint32_t>::checkInput(const char *dirpath, tT trID_key, uint16_t total_batches);

template void inputHandler<uint32_t>::clearInput();

template void inputHandler<uint64_t>::clearInput();

template void inputHandler<uint32_t>::resetInput();

template void inputHandler<uint64_t>::resetInput();

template void
inputStream<uint32_t>::loadBatch(std::vector<std::pair<uint32_t, uint32_t>> &lsh_enc_vec, unsigned int curr_batch);

template void
inputStream<uint64_t>::loadBatch(std::vector<std::pair<uint32_t, uint64_t>> &lsh_enc_vec, unsigned int curr_batch);

template void inputStream<uint64_t>::loadCounts(std::unordered_map<uint64_t, uint32_t> &rcounts);

template void inputStream<uint32_t>::loadCounts(std::unordered_map<uint32_t, uint32_t> &rcounts);

template uint64_t
inputStream<uint32_t>::retrieveBatch(vvec<uint32_t> &td, uint32_t tbatch_size, unsigned int curr_batch, bool shared_table);

template uint64_t
inputStream<uint64_t>::retrieveBatch(vvec<uint64_t> &td, uint32_t tbatch_size, unsigned int curr_batch, bool shared_table);

template void inputStream<uint32_t>::removeBatch(unsigned int curr_batch);

template void inputStream<uint64_t>::removeBatch(unsigned int curr_batch);

template std::map<uint8_t, uint64_t> inputHandler<uint64_t>::histRowSizes();

template std::map<uint8_t, uint64_t> inputHandler<uint32_t>::histRowSizes();
