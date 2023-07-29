#include "kestane.h"

int
main(int argc, char** argv)
{
  srand(0);
  /* srand(time(NULL)); */
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
  uint8_t h = 8;
  uint8_t b = 15;
  unsigned int read_batch_size = (unsigned int)DEFAULT_BATCH_SIZE;
  num_threads = 16;
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(1); // Use num_threads threads if explicitly stated

  std::cout << "-----------------" << std::endl;
  std::cout << "Length of k-mers : " << (int)k << std::endl;
  std::cout << "Number of positions for LSH : " << (int)h << std::endl;
  std::cout << "Number of threads : " << num_threads << std::endl;
  std::cout << "Read batch size : " << read_batch_size << std::endl;
  std::cout << "-----------------" << std::endl;

  maskLSH lsh_vg = generateMaskLSH(k, h);

  std::string str_pathA = "./test/test-A_and_B/counts-A";
  std::string str_pathB = "./test/test-A_and_B/counts-B";
  unsigned int num_reprs = 20;

  uint32_t num_batch_rows = pow(2, 16);
  HTd<uint64_t> tdA(k, h, num_batch_rows, &lsh_vg);
  HTd<uint64_t> tdB(k, h, num_batch_rows, &lsh_vg);

  uint64_t read_num_kmersA = 0;
  uint64_t read_num_kmersB = 0;

  for (unsigned int i = 1; i <= num_reprs; ++i) {
    /* std::cout << str_pathA + std::to_string(i) + ".fa" */
    /*           << " & " << str_pathB + std::to_string(i) + ".fa" << std::endl; */
    StreamIM<uint64_t> tcAi(k, h, &lsh_vg);
    StreamIM<uint64_t> tcBi(k, h, &lsh_vg);

    uint64_t total_num_kmersAi = tcAi.readBatch((str_pathA + std::to_string(i) + ".fa").c_str(), read_batch_size);
    uint64_t total_num_kmersBi = tcBi.readBatch((str_pathB + std::to_string(i) + ".fa").c_str(), read_batch_size);

    HTd<uint64_t> tdAi(k, h, num_batch_rows, &lsh_vg);
    HTd<uint64_t> tdBi(k, h, num_batch_rows, &lsh_vg);

    read_num_kmersA += tcAi.getBatch(tdAi.enc_vvec, num_batch_rows);
    read_num_kmersB += tcBi.getBatch(tdBi.enc_vvec, num_batch_rows);
    tdAi.updateSize();
    tdBi.updateSize();
    tdAi.initCounts();
    tdBi.initCounts();

    tcAi.save(("./test/test-A_and_B/lsh_enc-A" + std::to_string(i)).c_str());
    tcBi.save(("./test/test-A_and_B/lsh_enc-B" + std::to_string(i)).c_str());

    tdA.mergeRows(tdAi, true);
    tdB.mergeRows(tdBi, true);
    tdA.sortColumns();
    tdB.sortColumns();
    /* std::cout << "After merge (A): " << tdA.num_kmers << "/" << read_num_kmersA << std::endl; */
    /* std::cout << "After merge (B): " << tdB.num_kmers << "/" << read_num_kmersB << std::endl; */
    tdA.makeUnique(true);
    tdB.makeUnique(true);
    /* std::cout << "After make unique (A): " << tdA.num_kmers << "/" << read_num_kmersA << std::endl; */
    /* std::cout << "After make unique (B): " << tdB.num_kmers << "/" << read_num_kmersB << std::endl; */
    /* std::cout << std::endl; */
  }

  for (unsigned int i = 1; i <= num_reprs; ++i) {
    /* std::cout << str_pathA + std::to_string(i) + ".fa" */
    /*           << " & " << str_pathB + std::to_string(i) + ".fa" << std::endl; */
    StreamOD<uint64_t> tsAi(("./test/test-A_and_B/lsh_enc-A" + std::to_string(i)).c_str());
    StreamOD<uint64_t> tsBi(("./test/test-A_and_B/lsh_enc-B" + std::to_string(i)).c_str());

    HTd<uint64_t> tdAi(k, h, num_batch_rows, &lsh_vg);
    HTd<uint64_t> tdBi(k, h, num_batch_rows, &lsh_vg);
    read_num_kmersA += tsAi.getBatch(tdAi.enc_vvec, num_batch_rows);
    read_num_kmersB += tsBi.getBatch(tdBi.enc_vvec, num_batch_rows);
    tdAi.updateSize();
    tdBi.updateSize();
    tdAi.initCounts();
    tdBi.initCounts();

    tdA.unionRows(tdAi, true);
    tdB.unionRows(tdBi, true);
    /* std::cout << "After union (A): " << tdA.num_kmers << "/" << read_num_kmersA << std::endl; */
    /* std::cout << "After union (B): " << tdB.num_kmers << "/" << read_num_kmersB << std::endl; */
    tdA.makeUnique(true);
    tdB.makeUnique(true);
    /* std::cout << "After make unique (A): " << tdA.num_kmers << "/" << read_num_kmersA << std::endl; */
    /* std::cout << "After make unique (B): " << tdB.num_kmers << "/" << read_num_kmersB << std::endl; */
    /* std::cout << std::endl; */
  }

  tdA.updateSize();
  tdB.updateSize();

  /* std::cout << "Final: " << tdA.num_kmers << ", # of species: " << tdA.num_species << std::endl; */
  /* std::cout << "Final: " << tdB.num_kmers << ", # of species: " << tdB.num_species << std::endl; */

  HTs<uint64_t> tfA(k, h, b, num_batch_rows, &lsh_vg);
  tdA.transformHTs(tfA);

  HTs<uint64_t> tfB(k, h, b, num_batch_rows, &lsh_vg);
  tdB.transformHTs(tfB);

  tfA.mergeRows(tfB);
  tfA.mergeRows(tfB);
  tfA.makeUnique(true);

  // std::cout << "Af+Bf after make unique: (F) " << tfA.num_kmers << std::endl;

  // for (int i = 0; i < tfA.num_rows; ++i) {
  //   for (int j = 0; j < tfA.ind_arr[i]; ++j) {
  //     std::cout << tfA.scounts_arr[i * tfA.b + j] << std::endl;
  //   }
  // }

  tdA.mergeRows(tdB);
  tdA.mergeRows(tdB);
  tdA.makeUnique(true);

  std::cout << "Ad+Bd after make unique: (D) " << tdA.num_kmers << std::endl;

  for (int i = 0; i < tdA.num_rows; ++i) {
    for (auto val : tdA.scounts_vvec[i]) {
      std::cout << val << std::endl;
    }
  }

  return 0;
}
