#include "kestane.h"
#include <cstdint>
#include <iostream>

int
main(int argc, char** argv)
{
  srand(0);
  /* srand(time(NULL)); */
  unsigned int num_threads = 1;
  bool classified_out = false;
  bool unclassified_out = false;

  int cf_tmp;

  while (1) {
    static struct option long_options[] = {
      { "num-threads", 0, 0, 't' },
      { "unclassified-out", 0, 0, UNCLASSIFIED_OUT_OPT },
      { "classified-out", 0, 0, CLASSIFIED_OUT_OPT },
      { 0, 0, 0, 0 },
    };

    int option_ID = 0;
    cf_tmp = getopt_long(argc, argv, "t:", long_options, &option_ID);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    else {
      switch (cf_tmp) {
        case 't':
          num_threads = atoi(optarg); // Default is 1.
          break;
        case CLASSIFIED_OUT_OPT:
          classified_out = optarg;
        case UNCLASSIFIED_OUT_OPT:
          unclassified_out = optarg;
          break;
        case ':':
          printf("Missing option for '-%s'.\n", argv[optind - 2]);
          if (long_options[option_ID].has_arg == 1) {
            return 1;
          }
          break;
        case '?':
          if (optopt == 'i')
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
          else if (optopt == 'q')
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
          else if (isprint(optopt))
            fprintf(stderr, "Unknown option '-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
          return 1;
        default:
          abort();
      }
    }
  }

  uint8_t k = 32;
  uint8_t h = 15;
  uint8_t b = 3;
  unsigned int read_batch_size = (unsigned int)DEFAULT_BATCH_SIZE;
  num_threads = 128;
  omp_set_dynamic(0);               // Explicitly disable dynamic teams
  omp_set_num_threads(num_threads); // Use num_threads threads

  std::cout << "-----------------" << std::endl;
  std::cout << "Length of k-mers : " << (int)k << std::endl;
  std::cout << "Number of positions for LSH : " << (int)h << std::endl;
  std::cout << "Number of threads : " << num_threads << std::endl;
  std::cout << "Read batch size : " << read_batch_size << std::endl;
  std::cout << "-----------------" << std::endl;

  maskLSH lsh_vg = generateMaskLSH(k, h);

  char* fpathA1 = (char*)"./test/counts-A1.fa";
  char* fpathB1 = (char*)"./test/counts-B1.fa";
  char* fpathA2 = (char*)"./test/counts-A2.fa";
  char* fpathB2 = (char*)"./test/counts-B2.fa";

  tableC<uint64_t> tcA1(k, h, &lsh_vg);
  tableC<uint64_t> tcB1(k, h, &lsh_vg);
  tableC<uint64_t> tcA2(k, h, &lsh_vg);
  tableC<uint64_t> tcB2(k, h, &lsh_vg);

  uint64_t total_num_kmersA1 = tcA1.fillVec(fpathA1, read_batch_size);
  uint64_t total_num_kmersB1 = tcB1.fillVec(fpathB1, read_batch_size);
  uint64_t total_num_kmersA2 = tcA2.fillVec(fpathA2, read_batch_size);
  uint64_t total_num_kmersB2 = tcB2.fillVec(fpathB2, read_batch_size);

  std::unordered_map<uint8_t, uint64_t> hist_mapA1 = tcA1.histNumCols();
  std::unordered_map<uint8_t, uint64_t> hist_mapB1 = tcB1.histNumCols();
  std::unordered_map<uint8_t, uint64_t> hist_mapA2 = tcA2.histNumCols();
  std::unordered_map<uint8_t, uint64_t> hist_mapB2 = tcB2.histNumCols();

  tcA1.saveVec("./test/lsh_enc-A1");
  tcB1.saveVec("./test/lsh_enc-B1");
  tcA2.saveVec("./test/lsh_enc-A2");
  tcB2.saveVec("./test/lsh_enc-B2");

  uint32_t num_batch_rows = pow(2, 20);
  tableD<uint64_t> tdA(k, h, num_batch_rows, &lsh_vg);
  tableD<uint64_t> tdB(k, h, num_batch_rows, &lsh_vg);

  tableS<uint64_t> tsA1("./test/lsh_enc-A1");
  tableS<uint64_t> tsB1("./test/lsh_enc-B1");
  tableS<uint64_t> tsA2("./test/lsh_enc-A2");
  tableS<uint64_t> tsB2("./test/lsh_enc-B2");

  uint64_t read_num_kmersA1 = tsA1.getBatch(tdA.enc_vvec, num_batch_rows);
  uint64_t read_num_kmersB1 = tsB1.getBatch(tdB.enc_vvec, num_batch_rows);

  std::cout << "-----------------D----------------" << std::endl;
  tdA.updateSize();
  tdB.updateSize();
  std::cout << "Num-kmers in B after reading 1 :" << tdA.num_kmers << std::endl;
  std::cout << "Num-kmers in A after reading 1 :" << tdB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  uint64_t read_num_kmersA2 = tsA2.getBatch(tdA.enc_vvec, num_batch_rows);
  uint64_t read_num_kmersB2 = tsB2.getBatch(tdB.enc_vvec, num_batch_rows);

  std::cout << "-----------------" << std::endl;
  tdA.updateSize();
  tdB.updateSize();
  std::cout << "Num-kmers in B after reading 2 :" << tdA.num_kmers << std::endl;
  std::cout << "Num-kmers in A after reading2 :" << tdB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "Are rows sorted in A : " << tdA.areColumnsSorted() << std::endl;
  std::cout << "Are rows sorted in B : " << tdB.areColumnsSorted() << std::endl;

  tdA.sortColumns();
  tdB.sortColumns();
  std::cout << "After sorting, are rows sorted in A : " << tdA.areColumnsSorted() << std::endl;
  std::cout << "After sorting, are rows sorted in B : " << tdB.areColumnsSorted() << std::endl;
  std::cout << "-----------------" << std::endl;

  tdA.makeUnique();
  tdB.makeUnique();

  std::cout << "-----------------" << std::endl;
  tdA.updateSize();
  tdB.updateSize();
  std::cout << "Num-kmers in B after makeUnique :" << tdA.num_kmers << std::endl;
  std::cout << "Num-kmers in A after makeUnique :" << tdB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  std::cout << "-----------------" << std::endl;
  std::unordered_map<uint8_t, uint64_t> histA = tdA.histNumCols();
  for (auto kv : histA) {
    if (kv.second > 0)
      std::cout << "D-A) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }
  std::unordered_map<uint8_t, uint64_t> histB = tdB.histNumCols();
  for (auto kv : histB) {
    if (kv.second > 0)
      std::cout << "D-B) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }
  std::cout << "-----------------" << std::endl;

  std::chrono::steady_clock::time_point begin_pt = std::chrono::steady_clock::now();
  /* tdA.trimColumns(b); */
  /* tdB.trimColumns(b); */
  tdA.pruneColumns(b);
  tdB.pruneColumns(b);
  std::chrono::steady_clock::time_point end_pt = std::chrono::steady_clock::now();
  std::cout << "-----------------" << std::endl;
  std::cout << "Total trimming/pruning time " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_pt - begin_pt).count() << "[ms]" << std::endl;
  std::cout << "After trimming/prunning columns;" << std::endl;
  std::unordered_map<uint8_t, uint64_t> histfA = tdA.histNumCols();
  for (auto kv : histfA) {
    if (kv.second > 0)
      std::cout << "D-A) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }
  std::unordered_map<uint8_t, uint64_t> histfB = tdB.histNumCols();
  for (auto kv : histfB) {
    if (kv.second > 0)
      std::cout << "D-B) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }

  std::cout << "After pruning, are rows sorted in A : " << tdA.areColumnsSorted() << std::endl;
  std::cout << "After pruning, are rows sorted in B : " << tdB.areColumnsSorted() << std::endl;

  tdA.updateSize();
  tdB.updateSize();
  std::cout << "Num-kmers in B after pruning :" << tdA.num_kmers << std::endl;
  std::cout << "Num-kmers in A after pruning :" << tdB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  tableF<uint64_t> tfA(k, h, b, num_batch_rows, &lsh_vg);
  tdA.transformTableF(tfA);
  std::cout << "tableD-A, size is " << tdA.table_size << std::endl;
  std::cout << "tableF-A, size is " << tfA.table_size << std::endl;

  tableF<uint64_t> tfB(k, h, b, num_batch_rows, &lsh_vg);
  tdB.transformTableF(tfB);
  std::cout << "tableD-B, size is " << tdB.table_size << std::endl;
  std::cout << "tableF-B, size is " << tfB.table_size << std::endl;

  std::cout << "-----------------" << std::endl;
  std::unordered_map<uint8_t, uint64_t> histtA = tfA.histNumCols();
  for (auto kv : histtA) {
    if (kv.second > 0)
      std::cout << "F-A) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }
  std::unordered_map<uint8_t, uint64_t> histtB = tfB.histNumCols();
  for (auto kv : histtB) {
    if (kv.second > 0)
      std::cout << "F-B) Number of columns : count = " << (int)kv.first << " :" << kv.second << std ::endl;
  }

  std::cout << "---------D--------" << std::endl;
  std::cout << "Are rows sorted in A : " << tdA.areColumnsSorted() << std::endl;
  std::cout << "Are rows sorted in B : " << tdB.areColumnsSorted() << std::endl;

  std::cout << "-----------------F----------------" << std::endl;
  tfA.updateSize();
  tfB.updateSize();
  std::cout << "Num-kmers in B :" << tfA.num_kmers << std::endl;
  std::cout << "Num-kmers in A :" << tfB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  std::cout << "-----------------" << std::endl;
  std::cout << "Are rows sorted in A : " << tfA.areColumnsSorted() << std::endl;
  std::cout << "Are rows sorted in B : " << tfB.areColumnsSorted() << std::endl;

  tfA.sortColumns();
  tfB.sortColumns();
  std::cout << "After sorting, are rows sorted in A : " << tfA.areColumnsSorted() << std::endl;
  std::cout << "After sorting, are rows sorted in B : " << tfB.areColumnsSorted() << std::endl;
  std::cout << "-----------------" << std::endl;

  tfA.makeUnique();
  tfB.makeUnique();

  std::cout << "-----------------" << std::endl;
  tfA.updateSize();
  tfB.updateSize();
  std::cout << "Num-kmers in B after makeUnique : " << tfA.num_kmers << std::endl;
  std::cout << "Num-kmers in A after makeUnique : " << tfB.num_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  std::cout << "---------D--------" << std::endl;
  tdA.mergeRows(tdB);
  /* tdA.unionRows(tdB); */
  tdA.updateSize();
  std::cout << "Num-kmers in tdA after mergeRows : " << tdA.num_kmers << std::endl;
  tdA.makeUnique();
  tdA.updateSize();
  std::cout << "Num-kmers in tdA after makeUnique : " << tdA.num_kmers << std::endl;
  tdA.pruneColumns(b);
  tdA.updateSize();
  std::cout << "Num-kmers in tdA after pruneColumns : " << tdA.num_kmers << std::endl;

  std::cout << "---------F--------" << std::endl;
  tfA.mergeRows(tfB);
  /* tfA.unionRows(tfB); */
  tfA.updateSize();
  std::cout << "Num-kmers in tfA after mergeRows : " << tfA.num_kmers << std::endl;
  tfA.makeUnique();
  tfA.updateSize();
  std::cout << "Num-kmers in tfA after makeUnique : " << tfA.num_kmers << std::endl;

  return 0;
}
