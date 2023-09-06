#ifndef _TAXONOMY_H
#define _TAXONOMY_H

#include "common.h"
#include "io.h"
#include <unordered_map>

class TaxonomyNCBI
{
public:
  TaxonomyNCBI(const char* nodes_filepath);
  uint64_t getParent(uint64_t taxID);
  std::string getRank(uint64_t taxID);
  void printTaxonomyNCBI();

private:
  std::unordered_map<uint64_t, uint64_t> _parent_map;
  std::unordered_map<uint64_t, std::string> _rank_map;
};

template<typename T>
class TaxonomyRecord
{
private:
  T _num_input;
  T _num_nodes;
  std::unordered_map<T, std::string> _tID_to_rank;
  std::unordered_map<T, uint64_t> _tID_to_taxID;
  std::unordered_map<uint64_t, T> _taxID_to_tID;
  std::unordered_map<std::string, T> _input_to_tID;
  std::unordered_map<T, std::string> _tID_to_input;
  std::vector<T> _parent_vec;
  std::vector<uint8_t> _depth_vec;
  std::unordered_map<T, std::set<T>> _child_map;

public:
  TaxonomyRecord(const char* input_filepath, TaxonomyNCBI taxonomy);
  void printTaxonomyRecord();
  T getLowestCommonAncestor(T a, T b);
  uint64_t changeIDtax(T tID);
  T changeIDt(uint64_t taxID);
  bool isBasis(T tID);
  decltype(_tID_to_input)& tID_to_input() { return _tID_to_input; }
  decltype(_child_map)& child_map() { return _child_map; }
  bool saveTaxonomyRecord(const char* library_dirpath);
};

#endif
