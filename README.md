# KRANK
:bangbang: This tutorial and examples are outdated, and will be updated soon. :bangbang:

## Quickstart
KRANK consists of two main subprograms: `build` and `query`.
The subprogram `krank build [OPTIONS]` reads input reference sequences and builds a KRANK library.
Once you have a library (or multiple libraries),  `krank query [OPTIONS]` performs queries against it (or against multiple libraries at the same time).
If you would rather use available libraries, you can download relevant ones and skip building your custom library.  
Running each subprogram with `--help` will display available subprograms together with brief descriptions regarding the program and its arguments ([see](#Usage)).

### Querying sequences against a KRANK library

### Building a KRANK library
KRANK requires two input files: a taxonomy and a mapping from taxon IDs to filepaths (or FTP URLs) of input sequences.
The path to the directory in which the krank library will be created, (or updated) must be stated with `--library-dir` or `-l`.

#### Input sequences
The input reference sequences could be any type: assembled sequences (genomes, contigs, scaffolds, etc.) or sets of *k*-mers (but not a mixture of all these them).
Although for most users it is not practical or useful, it is possible to use an external tool to extract *k*-mer sets, such as Jellyfish, and give *k*-mer sets directly as the input data.
These references might be in both FASTA and FASTQ format. Both `gzip` compressed or raw files are allowed.
Whatever input type you prefer, you need to provide filepaths or FTP URLs (or a mixture of them).
All you need is a mapping between taxon IDs (i.e., species ID) and paths/URLs, which will be a tab-separated file.
Each line should look similar to this.
```
562	/path/to/file/escherichia_coli-0001.fq
562	/path/to/file/escherichia_coli-0002.fq.gz
54736 /path/to/file/salmonell_abongori-0001.fa
2287  https://ftp.ncbi.nlm.nih.gov/path/to/genome.fna.gz
```
Note that, you can specify more than one file per taxon ID, but not vice versa.
Full paths are always preferred., filenames are not relevant and can be anything.
This file must be given with the option `--input-file /path/to/file/mapping` argument.
If the input references are sets of *k*-mers, then use the `--input-kmers` flag (the default is the complementary option `--input-sequences`).
Use `-k` or `--kmer-length` to specify the *k*-mer lengths, and `-w` or `--window-length` for minimizer window length.
The default value for $w$ is $k+3$, and $k=28$ by default.
You might want to increase $k$ to 31, $k=28$ would result in lighter-weight libraries with slightly worse recall.
Simply set $k$ and $w$ to the same value if you do not want to use minimizers.
If `--input-kmers` is given, the length of each line must be equal to $w$.
Otherwise, it is undefined behavior.

#### Taxonomy
The taxonomy consists of two files, namely `nodes.dmp` and `names.dmp`, and the path to the directory containing both must be given to the option `-t` (or `--taxonomy-dir`).
The format of these file is given [here](https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/NCBI/sourcerepresentation.html#file3).
The latest version of the NCBI taxonomy can be found [here](https://ftp.ncbi.nih.gov/pub/taxonomy/).
The keys (the first column) in the `--input-file` must appear in the taxonomy, i.e., the first field in the `nodes.dmp` and `names.dmp`.

#### Example commands and usage
KRANK is a computationally intensive algorithm, but it is possible to build a library relatively fast if you have enough computational resources.
You can benefit from parallel processing a lot by setting `--num-threads` to the number of available cores you have.
Note that the overhead is very little, the speed-up will grow with the number of cores you utilize.

In order to keep memory usage feasible, there are two ways: use the `--on-disk` flag (default) and increase the number of batches using `--batch-size`.
Unless you have a very small reference dataset (e.g., a few hundred genomes), do not use `--in-memory` flag, and keep using the default (`--on-disk).
With `--on-disk`, KRANK first reads input and stores *k*-mer encodings and LSH keys in the library as separate files and reads from there in batches during the table building.
When the final library is built, you can delete them (files with the prefix `lsh_enc_vec-*`).

The default value for `--batch-size` is $3$, meaning that $3$ bits will be used for batching, and hence, there will be $2^3=8$ batches.
Increasing `--batch-size` will split the table rows into more batches, and memory usage will decrease exponentially.
If you have a cluster with many nodes, you can process each of these batches in parallel on different nodes, using multiple cores too (by submitting each batch's job to a different node).
Note that you will need to process each batch separately, running `krank build` as independent jobs for each batch.
Alternatively, you can just build each batch one by one.
You can configure this behavior using `--target-batch` and `--from-library` options as shown below.

First, you need to initialize a library, without giving the `--target-batch` option.
```bash
./krank build \
  --library-dir /path/to/directory --taxonomy-dmp /path/to/nodes.dmp --input-files /path/to/mapping.tsv \
  -k 28 -w 31 -h 12 -b 16 --batch-size 3 \
  --kmer-ranking 1 --adaptive_size \
   --on-disk --input-sequences \
  --num-threads 32
```
This will initialize the library, and save the metadata & taxonomy file together with encoding and LSH keys (since `--on-disk` is given).

Then, you can process 8 batches one by one using 32 cores as follows:
```bash
./krank build \
  --library-dir /path/to/directory --taxonomy-dmp /path/to/nodes.dmp --input-files /path/to/mapping.tsv \
  -k 28 -w 31 -h 12 -b 16 --batch-size 3 \
  --kmer-ranking 1 --adaptive_size \
  --from-library --on-disk --input-sequences --target-batch 0 \
  --num-threads 32
```
Note that, using 0 as the target process all batches one by one.
Alternatively, you can process a specific batch, and even in parallel as mentioned.
```bash
printf '%d\n' {1..8} | xargs -I{} -P8 ./krank build \
    --library-dir /path/to/directory --taxonomy-dmp /path/to/nodes.dmp --input-files /path/to/mapping.tsv \
    -k 28 -w 31 -h 12 -b 16 --batch-size 3 \
    --kmer-ranking 1 --adaptive_size \
    --from-library --on-disk --input-sequences --target-batch {} \
    --num-threads 1
```
Notice that, when an initialized library is targeted, the `--from-library` flag is given.
Otherwise (default: `--from-scratch`), k-mers would be processed again, and LSH keys would be overwritten, possibly computed with a new random seed.
Hence, when a specific target is going to be built, `--from-library` must be given.
Initialization and library building can be done in a single command, only if all batches are going to be built together one by one.
This can be achieved by running the above command with `--target-batch 0` without the `--from-library` flag.


#### Parameters and library size
Two parameters define the shape of the final hash table: `-h` and `-b` (`--num-positions` and `--num-columns`, respectively).
KRANK uses $h$ many of positions (i.e., bases) of each *k*-mer to compute the hash value.
This also determines the number of rows of the table ($2^{2h}$).
As $b$ gives the number of columns, the total number of *k*-mers stored in a table is given by $2^{2h}b$.
So increasing $h$ and $b$ will directly translate into more memory usage, based on this.
If you are not sure about these values, just use `-b 16` and set `-h` to $k-16$.
The resulting reference library will use $2^(2h)16(4+2)+2^(2h)$ bytes. For $h=14$, $k=30$, and $b=16$, KRANK will construct a 26Gb library with high sensitivity.

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
                              Save macthing information to --output-dir for each query, this flag is given by default.
  --num-threads UINT          Number of threads to use for OpenMP based parallelism.
```
