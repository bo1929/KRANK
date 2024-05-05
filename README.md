# KRANK

## Quickstart
KRANK consists of two main subprograms: `build` and `query`.
The subprogram `krank build [OPTIONS]` reads input reference sequences and builds a KRANK library.
Once you have a library (or multiple libraries),  `krank query [OPTIONS]` performs queries against it (or against multiple libraries at the same time).
If you would rather use available libraries, you can download relevant ones and skip building your custom library.
Running each subprogram with `--help` will display available subprograms together with brief descriptions regarding the program and its arguments ([see](#Usage)).

## Installation
KRANK currently does not have pre-built binaries. You need to compile with make. Requirements are as follows:
- libcurl library (this is only needed for downloading reference sequences using URLs, otherwise not needed),
- OpenMP for parallelism,
- C++11 support,
- zlib library.

That is all you need. Simply compile by running `make`.
```bash
git clone https://github.com/bo1929/KRANK.git
make -C KRANK
KRANK/krank build --help
KRANK/krank query --help
```

### Querying sequences against a KRANK library
KRANK can query against a single or multiple libraries.
Each KRANK library is essentially a directory containing a searchable data structure in binary format.
See [a list of available libraries](#available-libraries) to download one or multiple based on your needs and quickly get started with your queries.
In [the next section](#building-a-krank-library), we discuss the library building procedure for custom libraries, this might be useful if you have some specific new genomes that you want to be considered during the search.
Below, a basic command to query against one library (stored in the directory `$LIBRARY_DIRECTORY`) is given.
```bash
krank query \
  -l $LIBRARY_DIRECTORY -o $OUTPUT_DIRECTORY -q $QUERY_PATHS \
  --num-threads $NUM_THREADS
```
You can substitute `$NUM_THREADS` simply with the number of cores available.
The value of `$OUTPUT_DIRECTORY` for option  `-o` (or `--output-dir`)  should be the path to the directory where you wish the results to be saved.
The option `-q` (or `--query-file`) can be one of the following: i) a path to the FASTA/FASTQ file containing the query sequences/reads, ii) a path to a tab-separated file containing a two-column mapping from a query ID to the corresponding query filepath to the FASTA/FASTQ file.
```
query_1	./query1.fna
query_2	./query2.fq
```

As mentioned earlier, you can use multiple libraries to expand your reference set by specifying corresponding directory paths after `-l` (or `--library-dir`, e.g., `-l $LIBRARY_DIRECTORY1 $LIBRARY_DIRECTORY2`).

#### Outputs: sequence classifications and relative abundances
Querying against a file generates two output files by default: a taxonomic label assignment for each sequence (i.e., short/long read or contig) and a summary of the taxonomic composition of the entire query file, i.e., relative abundances of taxa at different ranks.

The naming convention for output reports consists of a prefix + query ID -- or prefix + filename without the extension if a single query file is given: `classification_info-*` and `abundance_profile-*`. The first couple of lines of the classification report should look like the below example.
```
SEQ_ID  RANK  TAXON_ID  TAXON_NAME  PREDICTION_SCORE  MATCH_SCORE
AIJN01000024.1-771547 species 201 Campylobacter lari  0.744 0.101
AIJN01000024.1-369645 species 201 Campylobacter lari  0.978 11.4
```
`SEQ_ID` simply corresponds to the sequence identifier in FASTA/FASTQ.
If `PREDICTION_SCORE` is not close to 1, it means that KRANK also found some *k*-mers related to other taxa.
If it is close to 1, the prediction is not ambiguous, i.e., there is no other potential assignment.
However, this does not necessarily imply that the assignment is *confident*.
`MATCH_SCORE` quantifies the confidence, and the highest possible value is the number of *k*-mers in the sequence (e.g., 122 for short reads in the default parameter configuration).
KRANK internally filters assignments with `MATCH_SCORE` smaller than 0.03.
You can disable this by setting `--tvote-threshold 0`.
All values greater than 1 could be considered sufficiently high.
Below 1 is still worth considering especially for novel/distant queries that are of interest.

For relative abundances, an example for the output report would look as below.
```
RANK  TAXON_ID  TAXON_NAME  READ_COUNT  READ_ABUNDANCE  CELL_ABUNDANCE
class 117743  Flavobacteriia  3 0.012 0.0256103
class 1236  Gammaproteobacteria 6 0.024 0.0495194
class 183924  Thermoprotei  2 0.008 0.00916362
class 186801  Clostridia  4 0.016 0.0277098
```
`READ_COUNT` columns are simply the counts of query sequences/reads that are assigned to a particular taxon.
`READ_ABUNDANCE` is the `READ_COUNT` value divided by the total number of sequences given in the query/sample.
`CELL_ABUNDANCE` is the corrected abundance value that incorporates genome sizes on top of `READ_ABUNDANCE`, and this is the value that you should probably consider.
Note that the neither `READ_ABUNDANCE` nor `CELL_ABUNDANCE` columns do not sum up to 1 for a fixed taxonomic rank (e.g., phylum).
This is because some sequences can be left unassigned to any taxon.
You can further normalize these values for each taxonomic level.
If you would like to benchmark KRANK against some other tools, you can use [OPAL](https://github.com/CAMI-challenge/OPAL), it has an option (`-n`) which automatically does this normalization and enables comparison between normalized and un-normalized profiles.

#### Available libraries
More and more diverse libraries will be added soon.
- [WoL-v1 dataset - archaeal and bacterial genomes - lightweight (6.25Gb) and w/ ranking v0.3.2](https://ter-trees.ucsd.edu/data/krank/library-WoLv1-k29w32h13b16s7-reps_adpt.tar.gz)

#### Recommendations for choosing the right set of libraries
Soon.

### Building a KRANK library
To build a reference library, KRANK requires two input files: a taxonomy and a mapping from taxon IDs to filepaths (or URLs for FTP) of input reference sequences.
The path to the directory in which the krank library will be created (or updated) should be specified with `--library-dir` or `-l`. Note that you have to create the empty directory first.

#### Input reference sequences
The input reference sequences could be any type: assembled sequences (genomes, contigs, scaffolds, etc.) or sets of *k*-mers[^1] (but not a mixture).
[^1]: Although for most users it is neither practical nor useful, it is possible to use an external tool to extract *k*-mer sets, such as Jellyfish, and give *k*-mer sets directly as the input data.
These reference sequences might be in both FASTA and FASTQ formats.
Both `gzip` compressed or raw files are allowed.
Whatever input type you prefer, you need to provide filepaths or URLs (or a mixture of them).
All you have to provide is a mapping between taxon IDs (i.e., species ID) and paths/URLs, which will be a tab-separated file.
Some lines from such a file would look similar to this.
```
562	/path/to/file/escherichia_coli-0001.fq
562	/path/to/file/escherichia_coli-0002.fq.gz
54736 /path/to/file/salmonell_abongori-0001.fa
2287  https://ftp.ncbi.nlm.nih.gov/path/to/genome.fna.gz
```
Note that, you can specify more than one file per taxon ID, but not vice versa.
Full paths are always preferred., filenames are not relevant and can be anything.
This file must be given with the option `--input-file /path/to/file/mapping` argument[^2].
[^2]: If the input references are sets of *k*-mers, then use the `--input-kmers` flag (the default is the complement option `--input-sequences`).
Use `-k` or `--kmer-length` to specify the *k*-mer lengths, and `-w` or `--window-length` for minimizer window length.
See [this section](#parameters-and-library-size) for a discussion of these parameters.
Simply set -k and -w to the same value if you do not want to use minimizers[^3].
[^3]: If `--input-kmers` is given, the length of each line must be equal to -w.
Otherwise, it is undefined behavior.

#### Taxonomy
The taxonomy consists of two files, namely `nodes.dmp` and `names.dmp`, and the path to the directory containing both must be given to the option `-t` (or `--taxonomy-dir`).
The format of these files is given [here](https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/NCBI/sourcerepresentation.html#file3).
The latest version of the NCBI taxonomy can be found [here](https://ftp.ncbi.nih.gov/pub/taxonomy/).
The keys (the first column) in the `--input-file` must appear in the taxonomy, i.e., the first field in the `nodes.dmp` and `names.dmp`.

#### Example commands and usage
Building a library is a relatively expensive but one-time operation.
It consists of two steps: library initialization and batch building.
The subprogram `krank build` can either initialize a library or construct all/some batches of the library.
To initialize a library from reference sequences using default parameters (6.5Gb lightweight mode), run the below command.
```bash
 krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --batch-size 7 --from-scratch --num-threads $NUM_THREADS
```
You can benefit from parallel processing to a great extent by setting `--num-threads` to the number of available cores you have.
Note that the overhead is very little, the speed-up will increase with the number of cores you use.
The option `--from-scratch` specifies that this command is intended to initialize a non-existing library.
Hence, KRANK will not be looking for an already initialized library at `$LIBRARY_DIRECTORY`.
The option `--batch-size` sets the number of batches that the table will be split in log scale.
For `--batch-size 7`, it is $2^7=128$.
This parameter is a bit nuanced and the optimal value will vary.
Luckily, it does not affect the classification performance, but using a value that is too low may result in an explosion in memory usage.
A value between 4 and 9 should work just fine.
If you have the computational resources available to build multiple batches in parallel (e.g., on different cluster nodes or on the same large-memory machine), higher values (e.g., 7 to 9) are recommended.
The next paragraph discusses how library batches are built in parallel or one by one.

After initializing the library, you will need to process each batch independently.
This can be achieved by running separate `krank build` commands by setting the `--target-batch`.
For instance, to construct the library for the first batch out of 128, run the below command.
```bash
krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch 1 \
  --num-threads $NUM_THREADS
```
Notice that, when an initialized library is targeted, the `--from-library` flag is given.
Otherwise (initialization with `--from-scratch`, default), input reference *k*-mers would be processed again, and hash keys would be overwritten, possibly computed with a new random seed.
Hence, the `--from-library` option must be given to build an already-initialized library.

Alternatively, you can just build all batches with a single command by setting the `--target-batch` to 0.
```bash
krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch 0 \
  --num-threads $NUM_THREADS
```
However, the running time would be 128 times more compared to the previous command used for a single batch.

Suppose the memory on a machine (a single cluster node or a personal computer) allows to run 4 batches in parallel.
Then, we could process batches in parallel by using `xargs` as below.
```bash
printf '%d\n' {1..128} | xargs -I{} -P4  krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch {} \
  --num-threads $NUM_THREADS
```

The most ideal scenario would be having access to a cluster, and submitting independent jobs for each batch, letting the scheduler do the job.
But first, you need to initialize a library, without giving the `--target-batch` option.

Note that, KRANK first reads input and stores *k*-mer encodings and hash keys in the library as separate files and loads from those files in batches during the table building.
Thus, when the final library is built, you can delete them (i.e., files with the prefix `lsh_enc_vec-*`) by running `find $LIBRARY_DIRECTORY -type f -name "*lsh_enc_vec*" -delete`.

#### Parameters and library size
Two parameters define the shape of the final hash table: `-h` and `-b` (`--num-positions` and `--num-columns`, respectively).
KRANK uses $h$ (out of $k$) of positions (i.e., bases) of each *k*-mer to compute the hash value.
The recommended value for $h$ is equal to $k-16$.
This value determines the number of rows of the table ($2^{2h}$).
Since $b$ is the number of columns, the total number of *k*-mers stored in a table is given by $2^{2h}b$.
So increasing $h$ and $b$ will directly translate into more memory use.
If you are not sure about these values, just use `-b 16` and set `-h` to $k-16$.
The resulting reference library will use $2^{2h}16(4+2)+2^{2h}$ bytes.
For $h=14$, $k=30$, and $b=16$, KRANK will construct a 26Gb library with high sensitivity.
Alternatively, a lighter-weight 6.5Gb version could be built by setting $h=13$ and $k=29$.

The other two flags are related to *k*-mer selection and gradual *k*-mer filtering.
Setting `--kmer-ranking 0` (default is 1) changes *k*-mer selection to random; and KRANK does not use its heuristic to filter *k*-mers.
The other flag is `--free-size`, which disables our sampling strategy (option `--adaptive-size` is on by default).
If you want to use multiple libraries to query against, using `--free-size` for one of them might improve the performance.

To emulate CONSULT-II and to entirely skip the hierarchical *k*-mer selection step (which is only suggested if you have a large number of reference genomes), use `--fast-mode`.
Otherwise, to gradually filter some *k*-mers to fit your input reference to the library size specified, use `--selection-mode`.

## Usage
Running `krank`, `krank build`, and `krank query` with `--help` will give you the list of options, their description, and default values.

### `krank`
```
KRANK version: v0.3.2
Memory-bound & accurate taxonomic classification and profiling.
Usage: ./krank [OPTIONS] SUBCOMMAND

Options:
  --help
  --log,--no-log{false}       Extensive logging, might be too much and helpful for troubleshooting.
  --verbose,--no-verbose{false}
                              Increased verbosity and progress report.
  --seed INT                  Random seed for the LSH and other parts that require randomness.

Subcommands:
  build                       Builds a reference library with given k-mers sets or reference genomes.
  query                       Query given sequences with respect to given reference libraries.
```
The option `--log` is likely to be more than needed and is probably not useful unless you have an error.
It would be helpful to include logs if you need to create a new issue.
The default verbosity should be sufficient (`--verbose`).
Note that these options (`--help`, `--verbose`, `--log` and `--seed`) cannot be given after the subcommand (e.g., `krank query --seed 0 [OPTIONS]` will throw an error).

### `krank build`
```
KRANK version: v0.3.2
Builds a reference library with given k-mers sets or reference genomes.
Usage: ./krank build [OPTIONS]

Options:
  --help
  -l,--library-dir TEXT REQUIRED
                              Path to the directory containing the library.
  -t,--taxonomy-dir TEXT:DIR REQUIRED
                              Path to the directory containing the taxonomy files (nodes.dmp and names.dmp).
  -i,--input-file TEXT:FILE REQUIRED
                              Path to the file containing paths and taxon IDs of references.
  --from-library,--from-scratch{false}
                              Are k-mers already encoded and stored in the library?
							  Default: --from-scratch, and it reads k-mer sets or sequences from given input paths.
							  If --from-library is given, KRANK will try to read k-mers from an already-initialized library.
  --input-kmers,--input-sequences{false}
                              Are given input files k-mers sets (extracted with some external tool) or sequences (genomes, contigs etc.)?
							  If sequences, k-mers sets will be extracted internally.
							  Ignored if --from-library given.
							  Default: --input-sequences.
  -k,--kmer-length UINT       Length of k-mers. Default: 29.
  -w,--window-length UINT     Length of minimizer window. Default: k+3.
  -h,--num-positions UINT     Number of positions for the LSH. Default: 13.
  -b,--num-columns UINT       Number of columns of the table. Default: 16.
  -s,--batch-size UINT        Number of bits to divide the table into batches. Default: 6, i.e., 64 batches.
  --target-batch UINT         The specific library batch to be built.
							  If 0, all batches will be processed one by one.
							  If not given, the library will only be initialized after reading the input data and encoding k-mers.
  --kmer-ranking ENUM:value in {random_kmer->0,representative_kmer->1} OR {0,1}
                              Which strategy will be used for k-mer ranking? (0: random_kmer, 1: representative_kmer)
  --adaptive-size,--free-size{false}
                              Use size constraint heuristic while gradually building the library.
  --num-threads UINT          Number of threads to use for OpenMP-based parallelism.
  --fast-mode,--selection-mode{false}
                              The default mode is --selection-mode which traverses the taxonomy and selects k-mers accordingly.
							  When --fast-mode is given, tree traversal will be skipped, and the final library will be built at the root.
							  With --kmer-ranking random_kmer, this is equivalent to CONSULT-II.
							  If --fast-mode is given, --adaptive-size will be ignored and have no effect.
							  Note  --fast-mode is significantly faster.
  --update-annotations,--build-tables{false}
                              When --update-annotations option is given, KRANK tries to update soft LCAs of k-mers by going over reference genomes.
							  Default --build-tables selects k-mers, builds tables, and also computes soft LCAs.
```

### `krank query`
```
KRANK version: v0.3.2
Query given sequences with respect to given reference libraries.
Usage: ./krank query [OPTIONS]

Options:
  --help
  -l,--library-dir TEXT ... REQUIRED
                              Path(s) to the directory containing the library. Note that multiple libraries could be given to this option.
  -o,--output-dir TEXT:DIR    Path to the directory to output result files. Default: the current working directory.
  -q,--query-file TEXT:FILE REQUIRED
                              Path to the tab-separated file containing paths and IDs of query FASTA/FASTQ files.
  --total-vote-threshold,--tvote-threshold FLOAT
                              The minimum total vote to classify, can be considered as a confidence threshold. Default: 0.03.
  --max-match-distance,--max-match-hdist UINT
                              The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.
  --save-match-info,--no-match-info{false}
                              Save matching information to --output-dir for each query, this flag is not given by default. There is no practical need to give this flag. This is for debugging purposes and maybe for alternative down-stream analyses.
  --num-threads UINT          Number of threads to use for OpenMP-based parallelism.
```
