#include "krank.h"

int main(int argc, char **argv)
{
  PRINT_VERSION
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace);

  CLI::App app{"Memory-bound & accurate taxonomic classification and profiling."};
  app.set_help_flag("--help");
  bool log = false;
  app.add_flag("--log,!--no-log", log, "Extensive logging, might be too much and helpful for troubleshooting.");
  bool verbose = true;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report.");
  app.require_subcommand();
  int seed = 0;
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness.");
  app.callback([&]() {
    if (app.count("--seed"))
      gen.seed(seed);
  });

  CLI::App *sub_build =
    app.add_subcommand("build", "Builds a reference library with given k-mers sets or reference genomes.");
  std::string library_dir;
  sub_build->add_option("-l,--library-dir", library_dir, "Path to the directory containing the library.")->required();
  std::string tax_dir;
  sub_build
    ->add_option(
      "-t,--taxonomy-dir", tax_dir, "Path to the directory containing the taxonomy files (nodes.dmp and names.dmp).")
    ->required()
    ->check(CLI::ExistingDirectory);
  std::string input_file;
  sub_build->add_option("-i,--input-file", input_file, "Path to the file containing paths and taxon IDs of references.")
    ->required()
    ->check(CLI::ExistingFile);
  bool from_library = false;
  sub_build->add_flag("--from-library,!--from-scratch",
                      from_library,
                      "Are k-mers already encoded and stored in the library? "
                      "Default: --from-scratch, and it reads k-mer sets or sequences from given input paths. "
                      "If --from-library is given, KRANK will try to read k-mers from an already-initialized library.");
  bool input_kmers = false;
  sub_build->add_flag(
    "--input-kmers,!--input-sequences",
    input_kmers,
    "Are given input files k-mers sets (extracted with some external tool) or sequences (genomes, contigs etc.)? "
    "If sequences, k-mers sets will be extracted internally. "
    "Ignored if --from-library given. "
    "Default: --input-sequences.");
  uint8_t k = 29;
  uint8_t w = k + 3;
  sub_build->add_option("-k,--kmer-length", k, "Length of k-mers. Default: 29.");
  sub_build->add_option("-w,--window-length", w, "Length of minimizer window. Default: k+3.");
  uint64_t max_memory = 12000;
  sub_build->add_option(
    "-m,--max-memory",
    max_memory,
    "The maximum memory limit for the final library (in megabyte - MB). "
    "This will override -h (--num-positions), -b (--num-columns) and -k (--kmer-length). "
    "KRANK will pick a good configuration of k, h and b values for you. "
    "Note that this is only the final library size. "
    "Memory requirements may exceed this limit but you can always decrease needed memory per batch by increasing -s.");
  uint8_t h = 13;
  sub_build->add_option("-h,--num-positions", h, "Number of positions for the LSH. Default: 13.");
  uint8_t b = 16;
  sub_build->add_option("-b,--num-columns", b, "Number of columns of the table. Default: 16.");
  uint8_t batch_bsize = 6;
  sub_build->add_option("-s,--batch-size",
                        batch_bsize,
                        "Number of bits to divide the table into batches. "
                        "Default: 6, i.e., 64 batches.");
  uint16_t target_batch = 0;
  bool only_init = false;
  sub_build->add_option(
    "--target-batch",
    target_batch,
    "The specific library batch to be built. "
    "If 0, all batches will be processed one by one. "
    "If not given, the library will only be initialized after reading the input data and encoding k-mers.");
  std::map<std::string, RankingMethod> map_ranking{{"random", random_kmer}, {"representative", representative_kmer}};
  RankingMethod ranking_method = representative_kmer;
  sub_build
    ->add_option("--kmer-ranking",
                 ranking_method,
                 "Which strategy will be used for k-mer ranking? (0: random_kmer, 1: representative_kmer) "
                 "Default: representative_kmer, selected based on coverage heuristic.")
    ->transform(CLI::CheckedTransformer(map_ranking, CLI::ignore_case));
  std::map<std::string, LabelsLCA> map_labels{{"hard", hard_lca}, {"soft", soft_lca}};
  LabelsLCA labels_lca = soft_lca;
  sub_build
    ->add_option("--lca",
                 labels_lca,
                 "This option determines LCA computation method for k-mer labels? (0: hard_lca, 1: soft_lca) "
                 "Default: soft_lca, computed using CONSULT-II's heuristic.")
    ->transform(CLI::CheckedTransformer(map_labels, CLI::ignore_case));
  bool adaptive_size = false;
  sub_build->add_flag(
    "--adaptive-size,!--free-size", adaptive_size, "Use size constraint heuristic while gradually building the library.");
  sub_build->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP-based parallelism.");
  bool fast_mode = true;
  sub_build->add_flag(
    "--fast-mode,!--selection-mode",
    fast_mode,
    "The optional mode is --selection-mode which traverses the taxonomy and selects k-mers accordingly. "
    "When --selection-mode is not given, tree traversal will be skipped, and the final library will be built at the root. "
    "With --kmer-ranking random_kmer, this is equivalent to CONSULT-II. "
    "If --selection-mode is not given, --adaptive-size will be ignored and have no effect. "
    "Note  --fast-mode is significantly faster.");
  bool update_annotations = false;
  sub_build->add_flag(
    "--update-annotations,!--build-tables",
    update_annotations,
    "When --update-annotations option is given, KRANK tries to update soft LCAs of k-mers by going over reference genomes. "
    "If the intermediate files are deleted, attempting to update annotations will result in an error. "
    /* "This will be done without rebuilding the tables, hence it will be very fast. " */
    /* "This might be particularly useful when parameters for soft LCA are changed. " */
    /* "Without a target batch given (using --target-batch), either options would be ignored. " */
    /* "KRANK would just initialize the library. " */
    "Default --build-tables selects k-mers, builds tables, and also computes soft LCAs.");
  bool remove_intermediate = true;
  sub_build->add_flag("--remove-intermediate,!--keep-intermediate",
                      remove_intermediate,
                      "When --keep-intermediate is given, KRANK will not delete batch data stored on the disk. "
                      "You may need to manually remove many small files (potentially tens of thousands). "
                      "This may be desired if one wants to use --update-annotations later.");
  sub_build->callback([&]() {
    if (sub_build->count("--fast-mode"))
      ranking_method = map_ranking["random_kmer"];
    if (!sub_build->count("--target-batch"))
      only_init = true;
    if ((sub_build->count("-m") + sub_build->count("--max-memory"))) {
      float coef_m = 17.3168587;
      float log_mm = log2(static_cast<float>(max_memory)) + coef_m;
      h = (log_mm - 3) / 2;
      b = pow(2, log_mm - 2 * h);
      k = h + 16;
      std::cout << "Configuration for the given maximum memory: k,h,b=" << +k << "," << +h << "," << +b << std::endl;
    }
    if (!(sub_build->count("-w") + sub_build->count("--window-length")))
      w = k + 3;
  });

  CLI::App *sub_query = app.add_subcommand("query", "Query given sequences with respect to given reference libraries.");
  std::vector<std::string> library_dir_v;
  sub_query
    ->add_option("-l,--library-dir",
                 library_dir_v,
                 "Path(s) to the directory containing the library. "
                 "Note that multiple libraries could be given to this option.")
    ->required();
  std::string output_dir = "./";
  sub_query
    ->add_option(
      "-o,--output-dir", output_dir, "Path to the directory to output result files. Default: the current working directory.")
    ->check(CLI::ExistingDirectory);
  std::string query_file;
  sub_query
    ->add_option(
      "-q,--query-file", query_file, "Path to the tab-separated file containing paths and IDs of query FASTA/FASTQ files.")
    ->required()
    ->check(CLI::ExistingFile);
  std::vector<float> tvote_threshold_v;
  sub_query->add_option(
    "--total-vote-threshold,--tvote-threshold",
    tvote_threshold_v,
    "The minimum total vote threshold(s) to classify, the order should match the order of the given libraries. Default: 0.05.");
  sub_query->callback([&]() {
    if (!(sub_query->count("--tvote-threshold") + sub_query->count("--total-vote-threshold"))) {
      tvote_threshold_v.resize(library_dir_v.size());
      for (unsigned int i = 0; i < tvote_threshold_v.size(); ++i) {
        tvote_threshold_v[i] = 0.05;
      }
    }
  });
  uint8_t max_match_hdist = 5;
  sub_query->add_option("--max-match-distance,--max-match-hdist",
                        max_match_hdist,
                        "The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.");
  bool save_match_info = false;
  sub_query->add_flag("--save-match-info,!--no-match-info",
                      save_match_info,
                      "Save matching information to --output-dir for each query, this flag is not given by default. "
                      "There is no practical need to give this flag. "
                      "This is for debugging purposes and maybe for alternative down-stream analyses.");
  sub_query->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP-based parallelism.");

  CLI11_PARSE(app, argc, argv);

  if (w < k)
    std::cerr << "The minimizer window size (-w) cannot be smaller than k-mer length (-k)." << std::endl;
  if (b < 2)
    std::cerr << "The number of columns of the hash table (-b) cannot besmaller than 2." << std::endl;
  if (h < 2)
    std::cerr << "The number of LSH positions (-h) cannot be smaller than 2." << std::endl;
  if (h > 16)
    std::cerr << "The number of LSH positions (-h) cannot be greater than or equal to 16." << std::endl;
  if (k > 32)
    std::cerr << "The maximum allowed k-mer length (-k) is 32." << std::endl;
  if (h >= k)
    std::cerr << "The number of LSH positions (-h) cannot be greater than or equal to k-mer length (-k)." << std::endl;
  if ((batch_bsize < 1) || (batch_bsize > (2 * h - 1)))
    std::cerr << "The number of batching-bits (-s) must be smaller than 2h and greater than 0." << std::endl;
  if (target_batch > pow(2, batch_bsize))
    std::cerr << "The given target batch index cannot be greater than the total number of batches (2^s)." << std::endl;
  if (max_match_hdist > k)
    std::cerr << "Maximum Hamming distance for a match cannot be greater than k-mer length." << std::endl;
  if ((w < k) || (b < 2) || (h < 2) || (h > 16) || (k > 32) || (h >= k) || (batch_bsize > (2 * h - 1)) ||
      (target_batch > pow(2, batch_bsize)) || (max_match_hdist > k)) {
    exit(EXIT_FAILURE);
  }
  if (tvote_threshold_v.size() != library_dir_v.size()) {
    std::cerr << "The number of total threshold values must be equal to the number of KRANK libraries." << std::endl;
    exit(EXIT_FAILURE);
  }
#ifndef SHORT_TABLE
  if ((k - h) > 16) {
    std::cerr << "In order to use compact k-mer encodings, k and h must satisfy (k-h) <= 16." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  if (adaptive_size && fast_mode)
    std::puts(
      "Since flag --fast-mode has been given, --adaptive-size will be ignored, fast mode cannot enforce a size constraint.");
  if (sub_build->count("--kmer-ranking") && fast_mode)
    std::puts(
      "Since flag --fast-mode has been given, --kmer-ranking will be ignored, k-mer selection will be inherently random.");
  if (only_init && (sub_build->count("--update-annotations") || sub_build->count("--build-tables")))
    std::puts("Since no target batch is given (using --target-batch), --build-tables/--update-annotations will be ignored.");
  if (batch_bsize > 10)
    std::puts("Having too many batches may cause issues since a separate directory will be created for each batch.");

  if (sub_build->parsed()) {
    Library l(library_dir.c_str(),
              tax_dir.c_str(),
              input_file.c_str(),
              k,
              w,
              h,
              b,
              ranking_method,
              labels_lca,
              adaptive_size,
              pow(2, 2 * h) * b,           // capacity
              pow(2, 2 * h - batch_bsize), // tbatch_size
              from_library,
              input_kmers,
              target_batch,
              only_init,
              update_annotations,
              fast_mode,
              remove_intermediate,
              verbose || log,
              log);
  }
  if (sub_query->parsed()) {
    Query q(library_dir_v,
            output_dir.c_str(),
            query_file.c_str(),
            tvote_threshold_v,
            max_match_hdist,
            save_match_info,
            verbose || log,
            log);
    q.perform(DEFAULT_BATCH_SIZE);
  }

  return 0;
}
