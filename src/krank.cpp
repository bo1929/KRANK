#include "krank.h"

int main(int argc, char **argv)
{
  PRINT_VERSION
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace);

  CLI::App app{"KRANK: a memory-bound taxonomic identification tool."};
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

  CLI::App *sub_build = app.add_subcommand("build", "Builds a referenece library with given k-mers sets.");
  std::string library_dir;
  sub_build->add_option("-l,--library-dir", library_dir, "Path to the directory containing the library.")
    ->required()
    ->check(CLI::ExistingDirectory);
  std::string taxonomy_dmp;
  sub_build->add_option("-t,--taxonomy-dmp", taxonomy_dmp, "Path to the file containing the taxonomy dmp.")
    ->required()
    ->check(CLI::ExistingFile);
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
  uint8_t k = 28;
  uint8_t w = k + 3;
  sub_build->add_option("-k,--kmer-length", k, "Length of k-mers. Default: 28.");
  sub_build->add_option("-w,--window-length", w, "Length of minimizer window. Default: k+3.");
  sub_build->callback([&]() {
    if (!(sub_build->count("-w") + sub_build->count("--window-length")))
      w = k + 3;
  });
  uint8_t h = 12;
  sub_build->add_option("-h,--num-positions", h, "Number of positions for the LSH. Default: 12.");
  uint8_t b = 16;
  sub_build->add_option("-b,--num-columns", b, "Number of columns of the table. Default: 16.");
  uint8_t batch_size = 2;
  sub_build->add_option("-s,--batch-size",
                        batch_size,
                        "Number of bits to divide the table into batches. "
                        "Default: 2, i.e., 4 batches.");
  uint8_t target_batch = 0;
  bool only_init = false;
  sub_build->add_option("--target-batch",
                        target_batch,
                        "The specific library batch to be built. "
                        "If 0, all batches will be processed one by one. "
                        "If not given, library will only be initialized after reading the input data and encoding k-mers.");
  sub_build->callback([&]() {
    if (!sub_build->count("--target-batch"))
      only_init = true;
  });
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
  sub_build->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP based parallelism.");

  CLI::App *sub_query = app.add_subcommand("query", "Performs query with respect to a given referenece library.");
  std::vector<std::string> library_dir_v;
  sub_query->add_option("-l,--library-dir", library_dir_v, "Path(s) to the directory containing the library.")->required();
  std::string output_dir = "./";
  sub_query
    ->add_option(
      "-o,--output-dir", output_dir, "Path to the directory to output result files. Default: the current working directory.")
    ->required()
    ->check(CLI::ExistingDirectory);
  std::string query_file;
  sub_query
    ->add_option(
      "-q,--query-file", query_file, "Path to the tab-seperated file containing paths and IDs of query FASTA/FASTQ files.")
    ->required()
    ->check(CLI::ExistingFile);
  uint8_t max_match_hdist = 5;
  sub_query->add_option("--max-match-distance,--max-match-hdist",
                        max_match_hdist,
                        "The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.");
  bool save_match_info = false;
  sub_query->add_flag("--save-match-info,!--no-match-info",
                      save_match_info,
                      "Save macthing information to --output-dir for each query, this flag is given by default.");
  sub_query->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP based parallelism.");

  CLI11_PARSE(app, argc, argv);

  if (w < k)
    std::cerr << "The minimizer window size w can not be smaller than k-mer length." << std::endl;
  if (b < 2)
    std::cerr << "The number of columns of the hash table, b, can not be smaller than 2." << std::endl;
  if (h < 2)
    std::cerr << "The number of positions for LSH, h, can not be smaller than 2." << std::endl;
  if (h >= k)
    std::cerr << "The number of positions for LSH, h, can not be greater than or equal to k-mer length." << std::endl;
  if (batch_size > (2 * h - 1))
    std::cerr << "The number of bits to determine the number of batches must be smaller than 2h." << std::endl;
  if (target_batch > pow(2, batch_size))
    std::cerr << "The given target batch index (starts from 1) is greater than the number of batches." << std::endl;
  if (max_match_hdist > k)
    std::cerr << "Maximum Hamming distance for a match can not be greater than k-mer length." << std::endl;
  if ((w < k) || (b < 2) || (h < 2) || (h >= k) || (batch_size > (2 * h - 1)) || (target_batch > pow(2, batch_size)) ||
      (max_match_hdist > k)) {
    exit(EXIT_FAILURE);
  }
#ifndef SHORT_TABLE
  if ((k - h) > 16) {
    std::cerr << "In order to use compact k-mer encodings, k and h must satisfy (k-h) <= 16." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  if (sub_build->parsed()) {
    std::cout << "Building the library..." << std::endl;
    Library l(library_dir.c_str(),
              taxonomy_dmp.c_str(),
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
              verbose,
              log);
    std::cout << "Library has been built and saved" << std::endl;
  }

  if (sub_query->parsed()) {
    std::cout << "Querying the given sequences..." << std::endl;
    Query q(library_dir_v, output_dir.c_str(), query_file.c_str(), max_match_hdist, save_match_info, verbose, log);
    q.run(DEFAULT_BATCH_SIZE);
    std::cout << "Results for the input queries have been saved" << std::endl;
  }

  return 0;
}
