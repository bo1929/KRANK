#include "krank.h"

int main(int argc, char **argv)
{
  PRINT_VERSION
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace);

  CLI::App app{"Memory-bound and accurate taxonomic classification and profiling."};
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
  std::string taxonomy_dir;
  sub_build->add_option("-t,--taxonomy-dir", taxonomy_dir, "Path to the directory containing the taxonomy files.")
    ->required()
    ->check(CLI::ExistingDirectory);
  std::string input_file;
  sub_build
    ->add_option("-i,--input-file", input_file, "Path to the file containing paths and taxon IDs of reference k-mer sets.")
    ->required()
    ->check(CLI::ExistingFile);
  bool from_library = false;
  sub_build->add_flag("--from-library,!--from-scratch",
                      from_library,
                      "Are k-mers already encoded and stored in the library? "
                      "Default: --from-scratch, and it reads k-mer sets or sequences from given input paths.");
  bool input_kmers = false;
  sub_build->add_flag("--input-kmers,!--input-sequences",
                      input_kmers,
                      "Are given input files k-mers sets or sequences? "
                      "If sequences, k-mers sets will be extracted internally. "
                      "Ignored if --from-library given. "
                      "Default: --input-sequences.");
  uint8_t k = 29;
  uint8_t w = k + 3;
  sub_build->add_option("-k,--kmer-length", k, "Length of k-mers. Default: 29.");
  sub_build->add_option("-w,--window-length", w, "Length of minimizer window. Default: k+3.");
  uint8_t h = 13;
  sub_build->add_option("-h,--num-positions", h, "Number of positions for the LSH. Default: 13.");
  uint8_t b = 16;
  sub_build->add_option("-b,--num-columns", b, "Number of columns of the table. Default: 16.");
  uint8_t batch_size = 7;
  sub_build->add_option("-s,--batch-size",
                        batch_size,
                        "Number of bits to divide the table into batches. "
                        "Default: 7, i.e., 128 batches.");
  uint16_t target_batch = 0;
  bool only_init = false;
  sub_build->add_option(
    "--target-batch",
    target_batch,
    "The specific library batch to be built. "
    "If 0, all batches will be processed one by one. "
    "If not given, the library will only be initialized after reading the input data and encoding k-mers.");
  std::map<std::string, RankingMethod> map_ranking{{"random_kmer", random_kmer},
                                                   {"representative_kmer", representative_kmer}};
  RankingMethod ranking_method = representative_kmer;
  sub_build
    ->add_option(
      "--kmer-ranking", ranking_method, "Which strategy will be used for k-mer ranking? (random_kmer, representative_kmer)")
    ->transform(CLI::CheckedTransformer(map_ranking, CLI::ignore_case));
  bool adaptive_size = false;
  sub_build->add_flag(
    "--adaptive-size,!--free-size", adaptive_size, "Use size constraint heuristic while gradually building the library.");
  sub_build->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP-based parallelism.");
  bool fast_mode = false;
  sub_build->add_flag(
    "--fast-mode,!--selection-mode",
    fast_mode,
    "The default mode is --selection-mode which traverses the taxonomy and selects k-mers accordingly. "
    "When --fast-mode is given, tree traversal will be skipped, and the final library will be built at the root. "
    "With --kmer-ranking random_kmer, this is equivalent to CONSULT-II. "
    "If --fast-mode is given, --adaptive-size will be ignored and have no effect. "
    "Note  --fast-mode is significantly faster.");
  bool update_annotations = false;
  sub_build->add_flag(
    "--update-annotations,!--build-tables",
    update_annotations,
    "When --update-annotations option is given, KRANK tries to update soft LCAs of k-mers by going over reference genomes. "
    "This will be done without rebuilding the tables, hence it would be quite fast. "
    "This might be particularly useful when parameters for soft LCA are changed. "
    "Without a target batch given (using --target-batch), both options would be ignored. "
    "Then, KRANK would only initialize the library. "
    "Default --build-tables selects k-mers, builds tables, and also computes soft LCAs.");
  sub_build->callback([&]() {
    if (sub_build->count("--fast-mode"))
      ranking_method = map_ranking["random_kmer"];
    if (!sub_build->count("--target-batch"))
      only_init = true;
    if (!(sub_build->count("-w") + sub_build->count("--window-length")))
      w = k + 3;
  });

  CLI::App *sub_query = app.add_subcommand("query", "Performs querying with respect to a given reference library.");
  std::vector<std::string> library_dir_v;
  sub_query->add_option("-l,--library-dir", library_dir_v, "Path(s) to the directory containing the library.")->required();
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
  float tvote_threshold = 0.03;
  sub_query->add_option("--total-vote-threshold,--tvote-threshold",
                        tvote_threshold,
                        "The minimum total vote to classify, can be considered as a confidence threshold. Default: 0.03.");
  uint8_t max_match_hdist = 5;
  sub_query->add_option("--max-match-distance,--max-match-hdist",
                        max_match_hdist,
                        "The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.");
  bool save_match_info = false;
  sub_query->add_flag("--save-match-info,!--no-match-info",
                      save_match_info,
                      "Save matching information to --output-dir for each query, this flag is not given by default.");
  sub_query->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP-based parallelism.");

  CLI11_PARSE(app, argc, argv);

  if (w < k)
    std::cerr << "The minimizer window size w can not be smaller than k-mer length." << std::endl;
  if (b < 2)
    std::cerr << "The number of columns of the hash table, b, can not be smaller than 2." << std::endl;
  if (h < 2)
    std::cerr << "The number of positions for LSH, h, can not be smaller than 2." << std::endl;
  if (h > 16)
    std::cerr << "The number of positions for LSH, h, can not be greater than or equal to 16." << std::endl;
  if (k > 32)
    std::cerr << "The maximum allowed k-mer length is 32." << std::endl;
  if (h >= k)
    std::cerr << "The number of positions for LSH, h, can not be greater than or equal to k-mer length." << std::endl;
  if ((batch_size < 1) || (batch_size > (2 * h - 1)))
    std::cerr << "The number of bits to determine the number of batches must be smaller than 2h and greater than 0."
              << std::endl;
  if (target_batch > pow(2, batch_size))
    std::cerr << "The given target batch index (starts from 1) is greater than the number of batches." << std::endl;
  if (max_match_hdist > k)
    std::cerr << "Maximum Hamming distance for a match can not be greater than k-mer length." << std::endl;
  if ((w < k) || (b < 2) || (h < 2) || (h > 16) || (k > 32) || (h >= k) || (batch_size > (2 * h - 1)) ||
      (target_batch > pow(2, batch_size)) || (max_match_hdist > k)) {
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
      "Since flag --fast-mode has been given, --adaptive-size will be ignored, fast mode can not enforce a size constraint.");
  if (sub_build->count("--kmer-ranking") && fast_mode)
    std::puts(
      "Since flag --fast-mode has been given, --kmer-ranking will be ignored, k-mer selection will be inherently random.");
  if (only_init && (sub_build->count("--update-annotations") || sub_build->count("--build-tables")))
    std::puts("Since no target batch is given (using --target-batch), --build-tables/--update-annotations will be ignored.");

  if (sub_build->parsed()) {
    if (only_init)
      std::cout << "Initializing the library..." << std::endl;
    else
      std::cout << "Building the library..." << std::endl;
    Library l(library_dir.c_str(),
              taxonomy_dir.c_str(),
              input_file.c_str(),
              k,
              w,
              h,
              b,
              ranking_method,
              adaptive_size,
              pow(2, 2 * h) * b,          // capacity
              pow(2, 2 * h - batch_size), // tbatch_size
              from_library,
              input_kmers,
              target_batch,
              only_init,
              update_annotations,
              fast_mode,
              verbose,
              log);
  }

  if (sub_query->parsed()) {
    std::cout << "Querying the given sequences..." << std::endl;
    Query q(
      library_dir_v, output_dir.c_str(), query_file.c_str(), tvote_threshold, max_match_hdist, save_match_info, verbose, log);
    q.perform(DEFAULT_BATCH_SIZE);
    std::cout << "Results for the input queries have been saved" << std::endl;
  }

  return 0;
}
