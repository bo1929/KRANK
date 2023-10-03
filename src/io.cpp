#include "io.h"

/* #define CANONICAL */

template<typename T>
inline void
sortColumns(vvec<T>& table)
{
  uint32_t num_rows = table.size();
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!table[rix].empty())
      std::sort(table[rix].begin(), table[rix].end());
  }
}

template<typename encT>
bool
StreamIM<encT>::save(const char* filepath)
{
  bool is_ok = true;
  if (StreamIM<encT>::lsh_enc_vec.empty()) {
    std::puts("The LSH-value and encoding pair vector is empty, nothing to save!\n");
    is_ok = false;
    return is_ok;
  }
  FILE* vec_f = IO::open_file(filepath, is_ok, "wb");
  std::fwrite(StreamIM<encT>::lsh_enc_vec.data(),
              sizeof(std::pair<uint32_t, encT>),
              StreamIM<encT>::lsh_enc_vec.size(),
              vec_f);
  if (std::ferror(vec_f)) {
    std::puts("I/O error when writing LSH-value and encoding pairs.\n");
    is_ok = false;
  }
  std::fclose(vec_f);
  return is_ok;
}

template<typename encT>
bool
StreamIM<encT>::load(const char* filepath)
{
  bool is_ok = false;
  StreamIM<encT>::lsh_enc_vec.clear();
  std::ifstream vec_ifs = IO::open_ifstream(filepath, is_ok);
  while (!vec_ifs.eof() && vec_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if (!vec_ifs.eof())
      StreamIM<encT>::lsh_enc_vec.push_back(lsh_enc);
  }
  if (vec_ifs.fail() && !vec_ifs.eof()) {
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (vec_ifs.eof()) {
    is_ok = true;
  }
  vec_ifs.close();
  return is_ok;
}

template<typename encT>
void
StreamIM<encT>::resetStream()
{
  l_rix = 0;
  curr_vix = 0;
}

template<typename encT>
void
StreamIM<encT>::clearStream()
{
  StreamIM<encT>::lsh_enc_vec.clear();
  l_rix = 0;
  curr_vix = 0;
}

template<typename encT>
uint64_t
StreamIM<encT>::getBatch(vvec<encT>& batch_table, uint32_t tbatch_size)
{
  assert(batch_table.size() >= tbatch_size);
  uint64_t num_kmers = 0;
  while (StreamIM<encT>::curr_vix < StreamIM<encT>::lsh_enc_vec.size()) {
    std::pair<uint32_t, encT> lsh_enc = StreamIM<encT>::lsh_enc_vec[StreamIM<encT>::curr_vix];
    if (lsh_enc.first < (StreamIM<encT>::l_rix + tbatch_size)) {
      batch_table[lsh_enc.first - StreamIM<encT>::l_rix].push_back(lsh_enc.second);
      num_kmers++;
    } else {
      break;
    }
    StreamIM<encT>::curr_vix++;
  }
  StreamIM<encT>::l_rix += tbatch_size;
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
void
StreamOD<encT>::openStream()
{
  is_open = true;
  vec_ifs = IO::open_ifstream(StreamOD::filepath, is_open);
  if (is_open)
    curr_pos = vec_ifs.tellg();
}

template<typename encT>
void
StreamOD<encT>::closeStream()
{
  is_open = false;
  vec_ifs.close();
}

template<typename encT>
uint64_t
StreamOD<encT>::getBatch(vvec<encT>& batch_table, uint32_t tbatch_size, bool contd)
{
  if (!contd)
    assert(batch_table.size() >= tbatch_size);
  uint64_t num_kmers = 0;
  uint64_t row_ix = 0;
  uint32_t curr_rix = StreamOD<encT>::curr_rix;
  uint32_t f_rix = StreamOD<encT>::f_rix;
  const size_t bufsize = 1024 * 1024;
  char buf[bufsize];
  vec_ifs.rdbuf()->pubsetbuf(buf, bufsize);
  while ((row_ix < tbatch_size) && (!StreamOD::vec_ifs.eof() && StreamOD::vec_ifs.good())) {
    std::pair<uint32_t, encT> lsh_enc;
    StreamOD::vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    curr_rix = lsh_enc.first;
    row_ix = curr_rix - f_rix;
    if (!StreamOD::vec_ifs.eof() && (row_ix < tbatch_size)) {
      if (!contd) {
        assert(row_ix < batch_table.size());
        batch_table[row_ix].push_back(lsh_enc.second);
      }
      StreamOD<encT>::curr_rix = curr_rix + 1;
      StreamOD<encT>::curr_pos = StreamOD<encT>::vec_ifs.tellg();
      num_kmers++;
    }
  }
  if (StreamOD<encT>::vec_ifs.fail() && !StreamOD<encT>::vec_ifs.eof()) {
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (StreamOD<encT>::vec_ifs.eof()) {
    StreamOD<encT>::vec_ifs.close();
  } else {
    StreamOD<encT>::vec_ifs.seekg(curr_pos);
    StreamOD<encT>::f_rix = f_rix + tbatch_size;
  }
  if (!contd)
    sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
uint64_t
StreamIM<encT>::extractInput(uint64_t rbatch_size)
{
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * h);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  uint64_t mask_bp = u64m >> (32 - k) * 2;
  uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  unsigned int i, l, c;
  for (std::string& filepath : filepath_v) {
    kseq_t* reader = IO::getReader(filepath.c_str());
    std::vector<sseq_t> seqBatch = IO::readBatch(reader, rbatch_size);
    while (!(seqBatch.empty())) {
      for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
        uint64_t enc64_bp;
        uint64_t enc64_lr;
        uint64_t cenc64_bp;
        uint64_t cenc64_lr;
        uint32_t cenc32_lr;
        uint32_t cenc32_bp;
        uint32_t rix;
        uint8_t ldiff = w - k + 1;
        std::string kmer_seq;
        uint8_t kix = 0;
        uint32_t wix = lsh_enc_vec.size();
        std::vector<std::pair<uint32_t, encT>> lsh_enc_win(ldiff);
        lsh_enc_vec.resize(wix + seqBatch[ix].len - k + 1);
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
                  lsh_enc_vec[wix] = std::make_pair(rix, cenc64_lr);
                  wix++;
                }
              } else if (std::is_same<encT, uint32_t>::value) {
                drop64Encoding32(*ptr_npositions, cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
                if (ldiff > 1) {
                  lsh_enc_win[kix % ldiff] = std::make_pair(rix, cenc32_lr);
                  kix++;
                } else {
                  lsh_enc_vec[wix] = std::make_pair(rix, cenc32_lr);
                  wix++;
                }
              } else {
                std::puts("Available encoding types are 'uint64_t' and 'uint32_t'.\n");
                exit(EXIT_FAILURE);
              }
            }
            if (l >= w && ldiff > 1) {
              lsh_enc_vec[wix] =
                *std::min_element(lsh_enc_win.begin(),
                                  lsh_enc_win.end(),
                                  [](std::pair<uint32_t, encT> lhs, std::pair<uint32_t, encT> rhs) {
                                    return murmur64(lhs.second) < murmur64(rhs.second);
                                    /* return lhs.second < rhs.second; */
                                  });
              wix++;
            }
          } else
            l = 0;
        }
        lsh_enc_vec.resize(wix);
      }
      if (seqBatch.size() == rbatch_size)
        seqBatch = IO::readBatch(reader, rbatch_size);
      else
        seqBatch.clear();
    }
  }
  lsh_enc_vec.shrink_to_fit();
  std::sort(lsh_enc_vec.begin(),
            lsh_enc_vec.end(),
            [](const std::pair<uint32_t, encT>& l, const std::pair<uint32_t, encT>& r) {
              if (l.first == r.first)
                return l.second < r.second;
              else
                return l.first < r.first;
            });
  lsh_enc_vec.erase(std::unique(lsh_enc_vec.begin(), lsh_enc_vec.end()), lsh_enc_vec.end());
  tnum_kmers = lsh_enc_vec.size();
  return tnum_kmers;
}

template<typename encT>
uint64_t
StreamIM<encT>::processInput(uint64_t rbatch_size)
{
  uint64_t tnum_kmers_sum = 0;
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * h);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  uint64_t mask_bp = u64m >> (32 - k) * 2;
  uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
  for (std::string& filepath : filepath_v) {
    kseq_t* reader = IO::getReader(filepath.c_str());
    std::vector<sseq_t> seqBatch = IO::readBatch(reader, rbatch_size);
    while (!(seqBatch.empty())) {
      uint8_t ldiff = w - k + 1;
      uint8_t kix = 0;
      std::vector<std::pair<uint32_t, encT>> lsh_enc_win(ldiff);
      uint32_t wix = lsh_enc_vec.size();
      lsh_enc_vec.resize(wix + seqBatch.size());
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
            kmerEncodingUpdate(
              seqBatch[ix].nseq.substr(kix + k - 1, 1).c_str(), enc64_lr, enc64_bp);
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
              lsh_enc_vec[wix + ix] = std::make_pair(rix, cenc64_lr);
          } else if (std::is_same<encT, uint32_t>::value) {
            drop64Encoding32(*ptr_npositions, cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
            if (ldiff > 1)
              lsh_enc_win[kix] = std::make_pair(rix, cenc32_lr);
            else
              lsh_enc_vec[wix + ix] = std::make_pair(rix, cenc32_lr);
          } else {
            std::puts("Available encoding types are 'uint64_t' and 'uint32_t'.\n");
            exit(EXIT_FAILURE);
          }
        }
        if (ldiff > 1)
          lsh_enc_vec[wix + ix] =
            *std::min_element(lsh_enc_win.begin(),
                              lsh_enc_win.end(),
                              [](std::pair<uint32_t, encT> lhs, std::pair<uint32_t, encT> rhs) {
                                return murmur64(lhs.second) < murmur64(rhs.second);
                                /* return lhs.second < rhs.second; */
                              });
      }
      tnum_kmers_sum += seqBatch.size();
      if (seqBatch.size() == rbatch_size)
        seqBatch = IO::readBatch(reader, rbatch_size);
      else
        seqBatch.clear();
    }
  }
  lsh_enc_vec.shrink_to_fit();
  std::sort(lsh_enc_vec.begin(),
            lsh_enc_vec.end(),
            [](const std::pair<uint32_t, uint64_t>& l, const std::pair<uint32_t, uint64_t>& r) {
              if (l.first == r.first)
                return l.second < r.second;
              else
                return l.first < r.first;
            });
  lsh_enc_vec.erase(std::unique(lsh_enc_vec.begin(), lsh_enc_vec.end()), lsh_enc_vec.end());
  tnum_kmers = lsh_enc_vec.size();
  if (tnum_kmers_sum != tnum_kmers)
    std::puts("Duplicate k-mers exist in the given input k-mer set.");
  return tnum_kmers;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
StreamIM<encT>::histRowSizes()
{
  std::unordered_map<uint32_t, uint8_t> row_sizes;
  for (auto kv : StreamIM<encT>::lsh_enc_vec) {
    row_sizes[kv.first]++;
  }
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (auto kv : row_sizes) {
    hist_map[kv.second]++;
  }
  return hist_map;
}

bool
IO::ensureDirectory(const char* dirpath)
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

kseq_t*
IO::getReader(const char* fpath)
{
  gzFile fp;
  fp = gzopen(fpath, "r");
  if (fp == nullptr) {
    std::cerr << "Failed to open file at " << fpath << std::endl;
    exit(1);
  }

  kseq_t* kseq;
  kseq = kseq_init(fp);

  return kseq;
}

uint64_t
IO::adjustBatchSize(uint64_t batch_size, uint8_t num_threads)
{
  return (batch_size / num_threads) * num_threads;
}

std::vector<sseq_t>
IO::readBatch(kseq_t* kseq, uint64_t batch_size)
{
  int l;
  unsigned int i = 0;
  std::vector<sseq_t> seqRead;
  seqRead.reserve(batch_size);
  while ((i < batch_size) && (l = kseq_read(kseq)) >= 0) {
    sseq_t tmp_seq;
    tmp_seq.nseq.assign(kseq->seq.s);
    tmp_seq.name.assign(kseq->name.s);
    tmp_seq.len = (kseq->seq.l);
    seqRead.push_back(tmp_seq);
    i++;
  }
  if ((i < batch_size) && (!seqRead.empty())) {
    kseq_destroy(kseq);
    gzclose(kseq->f->f);
  }
  return seqRead;
}

FILE*
IO::open_file(const char* filepath, bool& is_ok, const char* mode)
{
  FILE* f;
  f = std::fopen(filepath, mode);
  if (!f) {
    std::cerr << "File opening failed!" << filepath << std::endl;
    is_ok = false;
  }
  return f;
}

std::ifstream
IO::open_ifstream(const char* filepath, bool is_ok)
{
  std::ifstream ifs;
  ifs.open(filepath, std::ios::binary | std::ios::in);
  if (!ifs.is_open()) {
    std::fprintf(stderr, "File opening failed! %s", filepath);
    is_ok = false;
  }
  return ifs;
}

void
IO::read_encinput(std::string disk_path, std::vector<std::pair<uint32_t, encT>>& lsh_enc_vec)
{
  bool is_ok = false;
  std::ifstream vec_ifs = open_ifstream(disk_path.c_str(), is_ok);
  while (!vec_ifs.eof() && vec_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if (!vec_ifs.eof())
      lsh_enc_vec.push_back(lsh_enc);
  }
  if (vec_ifs.fail() && !vec_ifs.eof()) {
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (vec_ifs.eof()) {
    is_ok = true;
  }
  vec_ifs.close();
}

template uint64_t
StreamIM<uint32_t>::processInput(uint64_t rbatch_size);

template uint64_t
StreamIM<uint64_t>::processInput(uint64_t rbatch_size);

template uint64_t
StreamIM<uint32_t>::extractInput(uint64_t rbatch_size);

template uint64_t
StreamIM<uint64_t>::extractInput(uint64_t rbatch_size);

template bool
StreamIM<uint64_t>::save(const char* filepath);

template bool
StreamIM<uint32_t>::save(const char* filepath);

template bool
StreamIM<uint64_t>::load(const char* filepath);

template bool
StreamIM<uint32_t>::load(const char* filepath);

template void
StreamIM<uint32_t>::clearStream();

template void
StreamIM<uint64_t>::clearStream();

template void
StreamIM<uint32_t>::resetStream();

template void
StreamIM<uint64_t>::resetStream();

template uint64_t
StreamIM<uint64_t>::getBatch(vvec<uint64_t>& batch_table, uint32_t tbatch_size);

template uint64_t
StreamIM<uint32_t>::getBatch(vvec<uint32_t>& batch_table, uint32_t tbatch_size);

template void
StreamOD<uint64_t>::closeStream();

template void
StreamOD<uint32_t>::closeStream();

template void
StreamOD<uint64_t>::openStream();

template void
StreamOD<uint32_t>::openStream();

template uint64_t
StreamOD<uint64_t>::getBatch(vvec<uint64_t>& batch_table, uint32_t tbatch_size, bool cont);

template uint64_t
StreamOD<uint32_t>::getBatch(vvec<uint32_t>& batch_table, uint32_t tbatch_size, bool cont);

template std::unordered_map<uint8_t, uint64_t>
StreamIM<uint64_t>::histRowSizes();

template std::unordered_map<uint8_t, uint64_t>
StreamIM<uint32_t>::histRowSizes();
