# KRANK
:bangbang: This tutorial and examples are outdated, and will be updated soon. :bangbang:
## Quickstart
KRANK consists of two main subprograms: `build` and `query`.
You can extract *k*-mer sets from input sequences, and build a KRANK library using `krank build [OPTIONS]`.
Once you have a library (or multiple libraries), you can perform queries against it using `krank query [OPTIONS]`.
Running `krank --help` will display available subprograms together with brief descriptions regarding them and main arguments.

### Building a KRANK library
Please run `krank build --help` to see all options and defaults.
KRANK requires two input files: a taxonomy and a mapping from taxon IDs to paths of input sequences.
The path to the directory in which the krank library will be initialized, built, or updated must be stated with `--library-dir` or `-l`.

#### Input sequences
The input sequences could be assembled genomes (contigs or reads are also possible).
Alternatively, you can also use an external tool to extract *k*-mer sets, such as Jellyfish, and give *k*-mer sets directly as the input data.
Whatever input type you prefer, both FASTA and FASTQ formats can be used.
All you need is a mapping between taxon IDs and file paths, which will be a tab-separated file.
Each line should look similar to this.
```
562	/path/to/file/escherichia_coli-0001.fq
562	/path/to/file/escherichia_coli-0002.fq
54736 /path/to/file/salmonell_abongori-0001.fa
```
Note that, you can specify more than one file per taxon ID, but not vice versa.
Full paths are always preferred. File names are not relevant and can be anything.
This file must be specified with the `--input-file /path/to/file/mapping` argument.
If each line (2 lines actually, if FASTA) is a sequence (such as a scaffold) in the input data, then use the `--input-sequences` flag.
Otherwise, if each line corresponds to a *k*-mer extracted with an external tool, then use the `--input-kmers` flag (default).

Use `-k` or `--kmer-length` to specify the *k*-mer lengths, and `-w` or `--window-length` for minimizer window length.
The default value for $w$ is $k+3$, and $k=28$ by default.
You might want to increase $k$ to 31, $k=28$ would result in lighter-weight libraries with slightly worse precision.
Simply, set $k$ and $w$ to the same value if you do not want to use minimizers.
If `--input-kmers` is given, the length of each line must be equal to $w$.
Otherwise, it is undefined behavior.

#### Taxonomy
KRANK only expects a taxonomy tree, as a nodes file `nodes.dmp`, which you can download from the NCBI website.
The format of the file is given [here](https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/NCBI/sourcerepresentation.html#file3).
The keys (the first column) in the `--input-file` must appear in the taxonomy, i.e., the first field in the `nodes.dmp`.
Use `--taxonomy-dmp` or `-t` to give the path to the file, note that the file name is not important.

#### Library parameters
Two parameters define the shape of the final hash table: `-h` and `-b` (`--num-positions` and `--num-columns`, respectively).
KRANK uses $h$ many of positions (i.e., bases) of each *k*-mer to compute the hash value.
This also determines the number of rows of the table ($2^{2h}$).
As $b$ gives the number of columns, the total number of *k*-mers stored in a table is given by $2^{2h}b$.
So increasing $h$ and $b$ will directly translate into more memory usage, based on this.
If you are not sure about these values, just use `-b 16` and set `-h` to $k-16$.

The other two flags are related to *k*-mer selection.
Just use the defaults and do not change them if you do not know what you are doing.
Setting `--kmer-ranking 0` (default is 1) changes *k*-mer selection to random; and KRANK does not use our heuristic.
The other flag is `--free-size`, which disables our sampling strategy.
If you want to use multiple libraries to query against, using `--free-size` for one of them might improve the performance.

#### Logistics and usage details (with example commands)
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
