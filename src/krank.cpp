#include "krank.h"

int main(int argc, char **argv)
{
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace);

  CLI::App app{"KRANK: a memory-bound taxonomic identification tool."};
  app.set_help_flag("--help", "Note to self: WABM.");
  bool log = false;
  app.add_flag("--log,!--no-log", log, "Increased verbosity and extensive logging.");
  app.require_subcommand();

  CLI::App *sub_build = app.add_subcommand("build", "Builds a referenece library with given k-mers sets.");
  std::string library_dir;
  sub_build->add_option("-l,--library-dir", library_dir, "Path to the directory containing the library.")->required();
  std::string taxonomy_dmp;
  sub_build->add_option("-t,--taxonomy-dmp", taxonomy_dmp, "Path to the file containing the taxonomy dmp.")->required();
  std::string input_file;
  sub_build
    ->add_option("-i,--input-file", input_file, "Path to the file containing paths and taxon IDs of reference k-mer sets.")
    ->required();
  bool on_disk = true;
  sub_build->add_flag("--on-disk,!--in-memory",
                      on_disk,
                      "Are streams of k-mer sets going to be on the disk? "
                      "Default: --on-disk, note that --in-memory is memory hungry.");
  bool from_library = false;
  sub_build->add_flag("--from_library,!--from-scratch",
                      from_library,
                      "Are k-mers already encoded and stored in the library? "
                      "Default: --from-scratch, and it reads k-mers sets or sequences from given input paths.");
  bool from_kmers = false;
  sub_build->add_flag("--from-kmers,!--from-sequences",
                      from_kmers,
                      "Are given input files k-mers sets or sequences? "
                      "If sequences, k-mers sets will be extracted internally. "
                      "Default: --from-sequences.");
  uint8_t k = 32;
  sub_build->add_option("-k,--kmer-length", k, "Length of k-mers. Default: 32.");
  uint8_t w = k;
  sub_build->add_option("-w,--window-length", w, "Length of minimizer window. Default: k.");
  uint8_t h = 13;
  sub_build->add_option("-h,--num-positions", h, "Number of positions for the LSH. Default: 13.");
  uint8_t b = 8;
  sub_build->add_option("-b,--num-columns", b, "Number of columns of the table. Default: 8.");
  uint8_t batch_size = 3;
  sub_build->add_option("-s,--batch-size",
                        batch_size,
                        "Number of bits to divide the table into batches. "
                        "Default: 3, i.e., 8 batches.");
  uint8_t target_batch = 0;
  sub_build->add_option("--target-batch",
                        target_batch,
                        "The specific library batch to be built. "
                        "Default: 0, i.e., all batches will be processed one by one.");
  std::map<std::string, RankingMethod> map_ranking{{"random_kmer", random_kmer},
                                                   {"representative_kmer", representative_kmer}};
  RankingMethod ranking_method = {representative_kmer};
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
    ->add_option("-o,--output-dir",
                 output_dir,
                 "Path to the directory to output result files. "
                 "Default: the current working directory.")
    ->required();
  std::string query_file;
  sub_query
    ->add_option(
      "-q,--query-file", query_file, "Path to the tab-seperated file containing paths and IDs of query FASTA/FASTQ files.")
    ->required();
  uint8_t max_match_hdist = 5;
  sub_query->add_option("--max-match-distance,--max-match-hdist",
                        max_match_hdist,
                        "The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.");
  bool save_match_info = true;
  sub_query->add_flag("--save-match-info,!--no-match-info",
                      save_match_info,
                      "Save macthing information to --output-dir for each query, this flag is given by default.");
  sub_query->add_option("--num-threads", num_threads, "Number of threads to use for OpenMP based parallelism.");

  CLI11_PARSE(app, argc, argv);

  // Parameter check and display a report.
  if (!(sub_build->count("-w") + sub_build->count("--window-length"))) {
    w = k;
  }
  assert(on_disk || !from_library);

  if (sub_build->parsed()) {
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
              on_disk,
              from_kmers,
              target_batch,
              log);
  }

  if (sub_query->parsed()) {
    Query q(library_dir_v, output_dir.c_str(), query_file.c_str(), max_match_hdist, save_match_info, log);
  }

  return 0;
}
