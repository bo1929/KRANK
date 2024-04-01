# KRANK
:bangbang: This tutorial and examples are outdated, and will be updated soon. :bangbang:

## Quickstart
KRANK consists of two main subprograms: `build` and `query`.
The subprogram `krank build [OPTIONS]` reads input reference sequences and builds a KRANK library.
Once you have a library (or multiple libraries),  `krank query [OPTIONS]` performs queries against it (or against multiple libraries at the same time).
If you would rather use available libraries, you can download relevant ones and skip building your custom library.
Running each subprogram with `--help` will display available subprograms together with brief descriptions regarding the program and its arguments ([see](#Usage)).

### Querying sequences against a KRANK library
KRANK can query against a single or multiple libraries.
Each KRANK library is essentially a directory containing a searchable data structure in binary format.
See [Available libraries](#available-libraries) for ready-to-go libraries, download one or multiple based on your needs.
In [the next section](#building-a-krank-library), we discuss the library building procedure for custom libraries, this might be useful if you have some specific new genomes that you want to be considered during the search.
Below, a command to query against one library (stored at directory `$LIBRARY_DIRECTORY`) is given.
```bash
krank query \
  -l $LIBRARY_DIRECTORY -o $OUTPUT_DIRECTORY -q $QUERY_PATHS \
  --num-threads $NUM_THREADS
```
You can substitute `$NUM_THREADS` simply with the number of cores available.
The value of `$OUTPUT_DIRECTORY` for option  `-o` (or `--output-dir`)  should be the path to the directory which you wish the results to be saved.
The option `-q` (or `--query-file`) can be one of the following: i) a path to the FASTA/FASTQ file containing the query sequences/reads, ii) a path to a tab-separated file containing two column mapping from a query ID to the corresponding query filepath to the FASTA/FASTQ file.
```
query_1	./query1.fna
query_2	./query2.fq
```

As mentioned earlier, you can use multiple libraries to expand your reference set by specifying corresponding directory paths after `-l` (or `--library-dir`, e.g., `-l $LIBRARY_DIRECTORY1 $LIBRARY_DIRECTORY2`).

#### Outputs: sequence classifications and relative abundances
Querying against a file generate two output files by default: a taxonomic classification for each sequence (short/long read or contig) and a summary of the taxonomic composition of the entire query file, i.e., relative abundances of taxa at different levels.

Classification file is named as prefix + query ID (or prefix + filename without the extension if a single fine is given) : `classification_info-*`. The first couple lines should look like the following example.
```
SEQ_ID  RANK  TAXON_ID  TAXON_NAME  PREDICTION_SCORE  MATCH_SCORE
AIJN01000024.1-771547 species 201 Campylobacter lari  0.744 0.101
AIJN01000024.1-369645 species 201 Campylobacter lari  0.978 11.4
```
`SEQ_ID` is simply the identifier of the sequence in FASTA/FASTQ.
If `PREDICTION_SCORE` is not close to 1, it means that KRANK also found some *k*-mers related to other taxa.
If it is close to 1, the prediction is not ambiguous, i.e., there is no other potential assignment.
However, this does not necessarily imply that the assignment is *confident*.
`MATCH_SCORE` quantifies the confidence, and the highest possible value is the number of *k*-mers in the sequence (e.g., 122 for short reads in the default parameter configuration).
KRANK internally filters assignments with `MATCH_SCORE` smaller than 0.03.
You can disable this by setting `--tvote-threshold 0`.
Values greater than 1 are sufficiently high.
Below 1 is still worth considering especially novel/distant queries are of interest.

For relative abundances abundance, output files are prefixed with `abundance_profile-*`.
An example output would look as below.
```
RANK  TAXON_ID  TAXON_NAME  READ_COUNT  READ_ABUNDANCE  CELL_ABUNDANCE
class 117743  Flavobacteriia  3 0.012 0.0256103
class 1236  Gammaproteobacteria 6 0.024 0.0495194
class 183924  Thermoprotei  2 0.008 0.00916362
class 186801  Clostridia  4 0.016 0.0277098
```
`READ_COUNT` columns is simply the count of query sequences that are assigned to a particular taxon.
`READ_ABUNDANCE` is the `READ_COUNT` value divided by the total number of reads given in the query/sample.
`CELL_ABUNDANCE` is the corrected abundance value which incorporates genomes sizes on top of `READ_ABUNDANCE`, and this is the value which you should probably consider.
Note that the neither `READ_ABUNDANCE` nor `CELL_ABUNDANCE` columns does not sum up to 1 for a fixed taxonomic rank (e.g., phylum).
This is because it is possible for some sequences to be left unassigned to any taxon.
You can further normalize these values for each taxonomic level.
If you would like to benchmark KRANK against some other tools, you can use [OPAL](https://github.com/CAMI-challenge/OPAL), it has an option (`-n`) which automatically does this normalization and enables comparison between normalized and un-normalized profiles.

#### Available libraries
Soon.

### Building a KRANK library
To build a reference library, KRANK requires two input files: a taxonomy and a mapping from taxon IDs to filepaths (or URLs for FTP) of input reference sequences.
The path to the directory in which the krank library will be created (or updated) should be specified with `--library-dir` or `-l`.

#### Input reference sequences
The input reference sequences could be any type: assembled sequences (genomes, contigs, scaffolds, etc.) or sets of *k*-mers (but not a mixture).
Although for most users it is neither practical nor useful, it is possible to use an external tool to extract *k*-mer sets, such as Jellyfish, and give *k*-mer sets directly as the input data.
These reference sequences might be in both FASTA and FASTQ format.
Both `gzip` compressed or raw files are allowed.
Whatever input type you prefer, you need to provide filepaths or URLs (or a mixture of them).
All you have to provide is a mapping between taxon IDs (i.e., species ID) and paths/URLs, which will be a tab-separated file.
Lines from such a file should look similar to this.
```
562	/path/to/file/escherichia_coli-0001.fq
562	/path/to/file/escherichia_coli-0002.fq.gz
54736 /path/to/file/salmonell_abongori-0001.fa
2287  https://ftp.ncbi.nlm.nih.gov/path/to/genome.fna.gz
```
Note that, you can specify more than one file per taxon ID, but not vice versa.
Full paths are always preferred., filenames are not relevant and can be anything.
This file must be given with the option `--input-file /path/to/file/mapping` argument.
If the input references are sets of *k*-mers, then use the `--input-kmers` flag (the default is the complement option `--input-sequences`).
Use `-k` or `--kmer-length` to specify the *k*-mer lengths, and `-w` or `--window-length` for minimizer window length.
See [this section](#parameters-and-library-size) for a discussion of these parameters.
Simply set $k$ and $w$ to the same value if you do not want to use minimizers.
If `--input-kmers` is given, the length of each line must be equal to $w$.
Otherwise, it is undefined behavior.

#### Taxonomy
The taxonomy consists of two files, namely `nodes.dmp` and `names.dmp`, and the path to the directory containing both must be given to the option `-t` (or `--taxonomy-dir`).
The format of these file is given [here](https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/NCBI/sourcerepresentation.html#file3).
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
Hence, KRANK will not be looking for a already initialized library at `$LIBRARY_DIRECTORY`.
The option `--batch-size` sets the number of batches that the table will be split in log scale.
For `--batch-size`, it is $2^7=128$.
It is a bit nuanced and the optimal value will vary.
Luckily, it does not effect the classification performance, but using a value that is too low may result in explosion in memory usage.
A value between 4 and 9 should work fine.
If you have the computational resources available to build multiple batches in parallel (e.g., on different cluster nodes or on the same large-memory machine), higher values (e.g., 7 to 9) are recommended.
Next paragraph discusses how library batches are built in parallel or one-by-one.

After initializing the library, you will need to process each batch independently.
This can be achieved by running separate `krank build` commands by setting the `--target-batch`.
For instance to construct the library for the first batch out of 128, run the below command.
```bash
krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch 1 \
  --num-threads $NUM_THREADS
```
Notice that, when an initialized library is targeted, the `--from-library` flag is given.
Otherwise (initialization with `--from-scratch`), k-mers would be processed again, and LSH keys would be overwritten, possibly computed with a new random seed.
Hence, `--from-library` option must be given to build the initialized library.

Alternatively, you can just build all batches with a single command one by setting the `--target-batch` to 0.
```bash
krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch 1 \
  --num-threads $NUM_THREADS
```
But the running time would be 128 times more compared to the previous command used for a single batch.

Suppose the memory on a machine (a single cluster node or a personal computer) allows to run 4 batch in parallel.
Then, we could do this by using `xargs` as below.
```bash
printf '%d\n' {1..128} | xargs -I{} -P4  krank build \
   -l $LIBRARY_DIRECTORY -t $TAXONOMY_DIRECTORY -i $MAPPING_FILE \
   --from-library --batch-size 7 --target-batch {} \
  --num-threads $NUM_THREADS
```

The most desired scenario would be having access to a cluster, and submitting independent jobs for each batch, letting the scheduler do the job.
First, you need to initialize a library, without giving the `--target-batch` option.

Note that, KRANK first reads input and stores *k*-mer encodings and LSH keys in the library as separate files and reads from there in batches during the table building.
Hence, when the final library is built, you can delete them (i.e., files with the prefix `lsh_enc_vec-*`) by running `find $LIBRARY_DIRECTORY -type f -name "*lsh_enc_vec*" -delete`.

#### Parameters and library size
Two parameters define the shape of the final hash table: `-h` and `-b` (`--num-positions` and `--num-columns`, respectively).
KRANK uses $h$ many of positions (i.e., bases) of each *k*-mer to compute the hash value.
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

## Usage
Running `krank`, `krank build`, and `krank query` with `--help` will give you the list of options, their description, and default values.

### `krank`
```
Memory-bound and accurate taxonomic classification and profiling.
Usage: ./krank [OPTIONS] SUBCOMMAND

Options:
  --help
  --log,--no-log{false}       Extensive logging, might be too much and helpful for troubleshooting.
  --verbose,--no-verbose{false}
                              Increased verbosity and progress report.
  --seed INT                  Random seed for the LSH and other parts that require randomness.

Subcommands:
  build                       Builds a referenece library with given k-mers sets or reference genomes.
  query                       Performs query with respect to a given referenece library.

```
The option `--log` may output considerably many log messages and is probably not useful unless you have an error.
It would be helpful to include logs while creating a new issue.
The default verbosity should be sufficient.

### `krank build`
```
Usage: ./krank build [OPTIONS]

Options:
  --help
  -l,--library-dir TEXT:DIR REQUIRED
                              Path to the directory containing the library.
  -t,--taxonomy-dir TEXT:DIR REQUIRED
                              Path to the directory containing the taxonomy files.
  -i,--input-file TEXT:FILE REQUIRED
                              Path to the file containing paths and taxon IDs of reference k-mer sets.
  --from-library,--from-scratch{false}
                              Are k-mers already encoded and stored in the library? Default: --from-scratch, and it reads k-mer sets or sequences from given input paths.
  --input-kmers,--input-sequences{false}
                              Are given input files k-mers sets or sequences? If sequences, k-mers sets will be extracted internally. Ignored if --from-library given. Default: --input-sequences.
  -k,--kmer-length UINT       Length of k-mers. Default: 28.
  -w,--window-length UINT     Length of minimizer window. Default: k+3.
  -h,--num-positions UINT     Number of positions for the LSH. Default: 12.
  -b,--num-columns UINT       Number of columns of the table. Default: 16.
  -s,--batch-size UINT        Number of bits to divide the table into batches. Default: 2, i.e., 4 batches.
  --target-batch UINT         The specific library batch to be built. If 0, all batches will be processed one by one. If not given, library will only be initialized after reading the input data and encoding k-mers.
  --kmer-ranking ENUM:value in {random_kmer->0,representative_kmer->1} OR {0,1}
                              Which strategy will be used for k-mer ranking? (random_kmer, representative_kmer)
  --adaptive-size,--free-size{false}
                              Use size constraint heuristic while gradually building the library.
  --num-threads UINT          Number of threads to use for OpenMP based parallelism.
  --fast-mode,--selection-mode{false}
                              The default mode is --selection-mode which traverses the taxonomy, and selects k-mers accordingly.
                              When --fast-mode is given, tree traversal will be skipped, and the final library will be built at the root.
                              With --kmer-ranking random_kmer, this is equivalent to CONSULT-II.
                              If --fast-mode is given, --adaptive-size will be ignored and has no effect.
                              Note  --fast-mode is significantly faster.
  --update-annotations,--build-tables{false}
                              When --update-annotations is given, KRANK tries to update soft LCAs of k-mers by going over reference genomes.
                              This will be done without rebuilding the tables, hence it would be quite fast.
                              This might be particularly useful when parameters for soft LCA is changed.
                              Without a target batch given (using --target-batch), both options would be ignored.
                              Then, KRANK would only initialize the library.
                              Default --build-tables selects k-mers, builds tables, and also computes soft LCAs.
```
To emulate CONSULT-II and to skip the hierarchical *k*-mer selection step (which is only suggested if you have a large number of reference genomes), use `--fast-mode`.
Otherwise, to gradually filter some *k*-mers to fit your input to the library size specified, use `--selection-mode`.

### `krank query`
```
Performs query with respect to a given referenece library.
Usage: ./krank query [OPTIONS]

Options:
  --help
  -l,--library-dir TEXT ... REQUIRED
                              Path(s) to the directory containing the library.
  -o,--output-dir TEXT:DIR REQUIRED
                              Path to the directory to output result files. Default: the current working directory.
  -q,--query-file TEXT:FILE REQUIRED
                              Path to the tab-seperated file containing paths and IDs of query FASTA/FASTQ files.
  --total-vote-threshold,--tvote-threshold FLOAT
                              The minimum total vote to classify, can be considered as a confidence threshold. Default: 0.03.
  --max-match-distance,--max-match-hdist UINT
                              The maximum Hamming distance for a k-mer to be considered as a match. Default: 5.
  --save-match-info,--no-match-info{false}
                              Save macthing information to --output-dir for each query, this flag is not given by default.
  --num-threads UINT          Number of threads to use for OpenMP based parallelism.
```
