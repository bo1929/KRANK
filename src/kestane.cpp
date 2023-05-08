#include "kestane.h"

int
main(int argc, char** argv)
{
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

  /* char* fpath = (char*)"/expanse/lustre/projects/uot138/asapci/CONSULT-II/data_-_RECOMB-CG/kmers-k35C_minimized32C.fa.gz"; */
  char* fpath = (char*)"./test/counts.fa";
  uint8_t k = 32;
  uint8_t h = 12;
  unsigned int batch_size = (unsigned int)DEFAULT_BATCH_SIZE;
  num_threads = 128;
  double th_cdist = 2.1;

  std::cout << "-----------------" << std::endl;
  std::cout << "Length of k-mers : " << (int)k << std::endl;
  std::cout << "Number of positions for LSH : " << (int)h << std::endl;
  std::cout << "Number of threads : " << (int)num_threads << std::endl;
  std::cout << "Batch size : " << (int)batch_size << std::endl;
  std::cout << "Distance threshold for hierarchical clustering: " << th_cdist << std::endl;
  std::cout << "-----------------" << std::endl;

  std::vector<std::vector<uint32_t>> index_table((pow(2, 2 * h))); // (pow(2, h));
  uint32_t num_kmers = fillIndexTable(fpath, index_table, k, h, batch_size, num_threads);

  std::vector<uint64_t> vecHist = computeHistNumCols(index_table);
  std::cout << "-----------------" << std::endl;
  std::cout << "Number of k-mers : Count" << std::endl;
  for (auto i = 0; i < vecHist.size(); ++i) {
    if (vecHist[i] != 0) {
      std::cout << i << " : " << vecHist[i] << std::endl;
    }
  }
  std::cout << "-----------------" << std::endl;

  uint64_t* enc_arr;
  retrieveEncodings(fpath, enc_arr, num_kmers, batch_size, num_threads);

  uint64_t num_total_kmers = 0;
  uint64_t num_total_clusters = 0;

  std::cout << "-----------------" << std::endl;
  std::cout << "# of clusters / # of k-mers" << std::endl;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (int ix_row = 0; ix_row < index_table.size(); ++ix_row) {
    std::vector<uint32_t> row_index_table = index_table[ix_row];
    int num_kmers = row_index_table.size();
    if (num_kmers > 1) {
      uint64_t* row_enc_arr = new uint64_t[num_kmers];
      double* dist_matrix = new double[(num_kmers * (num_kmers - 1)) / 2];

      int i, j;
      int h_ij = 0;
      for (i = 0; i < num_kmers; i++) {
        row_enc_arr[i] = enc_arr[row_index_table[i]];
      }
      for (i = 0; i < num_kmers; i++) {
        for (j = i + 1; j < num_kmers; j++) {
          // Compute the Hamming distance between encoding row_index_table[i] and row_index_table[j].
          dist_matrix[h_ij] = (double)computeHammingDistance(row_enc_arr[i], row_enc_arr[j]);
          h_ij++;
        }
      }

      int* merge = new int[2 * (num_kmers - 1)];
      double* height = new double[num_kmers - 1];
      int* labels = new int[num_kmers];

      hclust_fast(num_kmers, dist_matrix, HCLUST_METHOD_COMPLETE, merge, height);
      // stop clustering at step with custer distance >= cdist
      cutree_cdist(num_kmers, merge, height, th_cdist, labels);

      int num_clusters = 0;
      for (int i = 0; i < num_kmers; ++i) {
        if (num_clusters < (labels[i] + 1)) {
          num_clusters = (labels[i] + 1);
        }
      }
#pragma omp critical
      {
        std::cout << num_clusters << "/" << num_kmers << std::endl;
        num_total_kmers += num_kmers;
        num_total_clusters += num_clusters;
      }

      delete[] dist_matrix;
      delete[] merge;
      delete[] height;
      delete[] labels;
    } else {
#pragma omp critical
      {
        num_total_kmers += num_kmers;
        num_total_clusters += num_kmers;
      }
    }
  }
  delete[] enc_arr;
  std::cout << "Final summation ratio (# of clusters / # of k-mers) : " << num_total_clusters << "/" << num_total_kmers << std::endl;
  std::cout << "-----------------" << std::endl;

  return 0;
}
