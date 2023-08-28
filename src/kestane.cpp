#include "kestane.h"
#include "table.h"

void
testLCAtaxID(TaxonomyRecord<uint16_t> record, uint64_t t1, uint64_t t2)
{
  uint16_t t = record.getLowestCommonAncestor(record.changeIDt(t1), record.changeIDt(t2));
  std::cout << "LCA(" << t1 << "," << t2 << ")";
  std::cout << " = " << record.changeIDtax(t) << std::endl;
}

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
  uint8_t h = 12;
  uint8_t b = 7;
  uint32_t tbatch_size = pow(2, 2 * h) / 128;
  uint64_t capacity = pow(2, 2 * h) * b / 2;
  num_threads = 32;
  maskLSH lsh_vg = generateMaskLSH(k, h);
  Library lib("./test/library/", "./test/Nodes10000.dmp", "./test/test_inputpath.tsv", k, h, b, capacity, tbatch_size, true);
  /* HTd<encT> td_root(1, k, h, tbatch_size, &lsh_vg, large_scount); */
  /* lib.getBatchHTd(td_root); */
  HTs<encT> ts_root(1, k, h, b, tbatch_size, &lsh_vg, random_kmer);
  lib.getBatchHTs(ts_root, 0, 3);

  return 0;
}
