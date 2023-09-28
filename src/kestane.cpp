#include "kestane.h"

int
main(int argc, char** argv)
{
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace);

  CLI::App app{ "kestane: a memory-efficient taxonomic identification tool" };
  app.set_help_flag("--help", "Note to self: WABM.");
  bool log = false;
  app.add_flag("--log,!--no-log", log, "Increased verbosity and extensive logging.");
  app.require_subcommand();

  CLI::App* sub_build =
    app.add_subcommand("build", "Builds a referenece library with given k-mers sets.");
  std::string library_dir;
  sub_build
    ->add_option("-l,--library-dir", library_dir, "Path to the directory containing the library.")
    ->required();
  std::string taxonomy_dmp;
  sub_build
    ->add_option("-t,--taxonomy-dmp", taxonomy_dmp, "Path to the file containing the taxonomy dmp.")
    ->required();
  std::string input_file;
  sub_build
    ->add_option("-i,--input-file",
                 input_file,
                 "Path to the file containing paths and taxon IDs of "
                 "reference k-mer sets.")
    ->required();
  uint8_t k = 32;
  sub_build->add_option("-k,--kmer-length", k, "Length of k-mers. Default is 32.");
  uint8_t h = 11;
  sub_build->add_option("-h,--num-positions", h, "Number of positions for the LSH. Default is 14.");
  uint8_t b = 8;
  sub_build->add_option("-b,--num-columns", b, "Number of columns of the table. Default is 8.");
  uint8_t switch_depth = 4;
  sub_build->add_option("--switch-depth",
                        switch_depth,
                        "The depth in the taxonomic tree at which the type of the table changes "
                        "from static to dynamic. Default is 4.");
  uint8_t batch_size = 3;
  sub_build->add_option("-s,--batch-size",
                        batch_size,
                        "Number of bits to divide the table into batches. "
                        "Default is 3, i.e., 8 batches.");
  uint8_t specify_batch = 0;
  sub_build->add_option("--specify-batch",
                        specify_batch,
                        "The specific library batch to be built. Default is 0, "
                        "i.e., all batches will be processed.");
  bool in_library = false;
  sub_build->add_flag("--in-library,!--raw-input",
                      in_library,
                      "Are k-mers already processed and in the library? "
                      "Default is --raw-input, "
                      "reads FASTA k-mer sets.");
  bool on_disk = true;
  sub_build->add_flag("--on-disk,!--in-memory",
                      on_disk,
                      "Are stream of k-mers will be on the disk? If so, k-mers will be read in "
                      "batches from the disk. Default is true, --in-memory is memory hungry.");
  std::map<std::string, RankingMethod> map_ranking{ { "random_kmer", random_kmer },
                                                    { "species_count", species_count },
                                                    { "information_score", information_score },
                                                    { "taxa_count", taxa_count } };
  std::pair<RankingMethod, RankingMethod> ranking_methods = { random_kmer, random_kmer };
  sub_build
    ->add_option("--kmer-ranking",
                 ranking_methods,
                 "Which strategy will be used for k-mer ranking? Requires a pair of "
                 "(random_kmer, species_count, information_score, taxa_count) "
                 "respectively before and after --switch-depth.")
    ->transform(CLI::CheckedTransformer(map_ranking, CLI::ignore_case));
  sub_build->add_option(
    "--num-threads", num_threads, "Number of threads to use for OpenMP based parallelism.");
  bool run_build = true;
  sub_build->add_flag("--run,!--only-init",
                      run_build,
                      "Run library building in addition to initialization, given by default.");

  CLI::App* sub_query =
    app.add_subcommand("query", "Performs query with respect to a given referenece library.");
  std::vector<std::string> library_dir_v;
  sub_query
    ->add_option(
      "-l,--library-dir", library_dir_v, "Path(s) to the directory containing the library.")
    ->required();
  std::string output_dir = "./";
  sub_query
    ->add_option("-o,--output_dir",
                 output_dir,
                 "Path to the directory to output result files. Default is "
                 "the current working directory.")
    ->required();
  std::string query_file;
  sub_query
    ->add_option("-q,--query-file",
                 query_file,
                 "Path to the tab-seperated file containing paths and IDs of "
                 "query FASTA/FASTQ files.")
    ->required();
  uint8_t max_match_hdist = 5;
  sub_query->add_option("--max-match-distance,--max-match-hdist",
                        max_match_hdist,
                        "The maximum Hamming distance for a k-mer to be "
                        "considered as a match. Default is 5.");
  bool save_match_info = true;
  sub_query->add_flag("--save-match-info,!--no-match-info",
                      save_match_info,
                      "Save macthing information to --output-dir for each "
                      "query, this flag is given by default.");
  sub_query->add_option(
    "--num-threads", num_threads, "Number of threads to use for OpenMP based parallelism.");
  bool run_query = true;
  sub_query->add_flag("--run,!--only-init",
                      run_query,
                      "Run queries in addition to instance initialization, given by default.");

  CLI11_PARSE(app, argc, argv);

  if (sub_build->parsed()) {
    Library lib(library_dir.c_str(),
                taxonomy_dmp.c_str(),
                input_file.c_str(),
                k,
                h,
                b,
                ranking_methods.first,
                ranking_methods.second,
                pow(2, 2 * h) * b,          // capacity
                pow(2, 2 * h - batch_size), // tbatch_size
                in_library,
                on_disk,
                specify_batch,
                log);
    if (run_build)
      lib.run(switch_depth);
  }

  if (sub_query->parsed()) {
    Query query(
      library_dir_v, output_dir.c_str(), query_file.c_str(), max_match_hdist, save_match_info, log);
    if (run_query)
      query.run();
  }

  return 0;
}
