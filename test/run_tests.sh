##!/bin/bash
RANDOM_SEED=0
echo "========"
date

echo "KRANK is a taxonomic identification tool."
echo "This script will run multiple test to make sure that KRANK can run flawlessly on your system."
echo "You can configure KRANK by specifying each parameter - press ENTER to skip and use default test value."
echo "Note that these defaults are only for the test and might not perform sufficiently good."
sleep 5

if [ ! -d taxonomy ]; then
  mkdir -p taxonomy
  wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -xvf taxdump.tar.gz -C taxonomy
  rm -f taxdump.tar.gz
fi

echo "Please enter:"
read -p "Library directory (example: ./test-KRANK/library, default: testlib) -> " LIBDIR
read -p "Number of threads to use (example: 4, default: 16) -> " NTHREADS
read -p "Length of k-mers, k (default: 27) -> " MERLEN
read -p "Minimizer window size, w (default: 35) -> " WINSIZE
read -p "Number of LSH positions, h (default: 12) -> " NUMPOS
read -p "Number of columns of the table, b (default: 8) -> " NUMCOL
read -p "Number of bits to split batches, s (default: 2) -> " NBATCHB

LIBDIR=${LIBDIR:-./testlib}
NTHREADS=${NTHREADS:-16}
MERLEN=${MERLEN:-27}
WINSIZE=${WINSIZE:-35}
NUMPOS=${NUMPOS:-12}
NUMCOL=${NUMCOL:-8}
NBATCHB=${NBATCHB:-2}

echo "Initialzing the library..."
/usr/bin/time -v ../krank --seed ${RANDOM_SEED} build \
  -l ${LIBDIR} -t ./taxonomy/ \
  -i ./input_map.tsv  \
  -k ${MERLEN} -w ${WINSIZE} -h ${NUMPOS} -b ${NUMCOL} -s ${NBATCHB} \
  --from-scratch --input-sequences \
  --kmer-ranking 1 --adaptive-size \
  --num-threads ${NTHREADS} \
  && echo "The library has been successfully initialized and k-mers havebeen extracted."

echo "Building the library for each batch one by one..."
/usr/bin/time -v ../krank --seed ${RANDOM_SEED} build \
  -l ${LIBDIR} -t ./taxonomy/ \
  -i ./input_map.tsv  \
  -k ${MERLEN} -w ${WINSIZE} -h ${NUMPOS} -b ${NUMCOL} -s ${NBATCHB} \
  --target-batch 0 --fast-mode --from-library --input-sequences \
  --kmer-ranking 1 --adaptive-size \
  --num-threads ${NTHREADS} \
  && echo "All batches have been constucted, the library is ready to query against."

echo "Building the library in fast mode, for each batch one by one..."
/usr/bin/time -v ../krank --seed ${RANDOM_SEED} build \
  -l ${LIBDIR} -t ./taxonomy/ \
  -i ./input_map.tsv  \
  -k ${MERLEN} -w ${WINSIZE} -h ${NUMPOS} -b ${NUMCOL} -s ${NBATCHB} \
  --target-batch 0 --selection-mode --from-library --input-sequences \
  --kmer-ranking 1 --adaptive-size \
  --num-threads ${NTHREADS} \
  && echo "All batches have been constucted, the library is ready to query against."

/usr/bin/time -v ../krank --seed ${RANDOM_SEED} query \
  -l ${LIBDIR} -o ./ -q query.fq \
  --max-match-distance 5 \
  --num-threads ${NTHREADS} \
  && echo "Queries have been completed, results have been saved."

date
echo "========"
