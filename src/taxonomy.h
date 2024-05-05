#ifndef _TAXONOMY_H
#define _TAXONOMY_H

#include "common.h"
#include "io.h"

class TaxonomyInput
{
private:
  std::unordered_map<uint32_t, uint32_t> _parent_map;
  std::unordered_map<uint32_t, uint8_t> _depth_map;
  std::unordered_map<uint32_t, std::string> _rank_map;
  std::unordered_map<uint32_t, std::string> _name_map;

public:
  TaxonomyInput(const char *tax_dirpath);
  uint32_t getParent(uint32_t tiID);
  std::string getRank(uint32_t tiID);
  void printTaxonomyInput();
  decltype(_parent_map) &parent_map() { return _parent_map; }
  decltype(_rank_map) &rank_map() { return _rank_map; }
  decltype(_depth_map) &depth_map() { return _depth_map; }
  decltype(_name_map) &name_map() { return _name_map; }
};

template<typename T>
class TaxonomyRecord
{
private:
  uint32_t _num_input;
  T _num_nodes;
  uint32_t _full_size;
  std::unordered_map<T, std::string> _trID_to_rank;
  std::unordered_map<tT, tT> _trID_to_lsroot;
  std::unordered_map<T, uint32_t> _trID_to_tiID;
  std::unordered_map<uint32_t, T> _tiID_to_trID;
  std::unordered_map<std::string, T> _input_to_trID;
  std::unordered_map<T, std::vector<std::string>> _trID_to_input;
  std::vector<T> _parent_vec;
  std::vector<uint8_t> _depth_vec;
  std::unordered_map<T, std::set<T>> _child_map;
  std::unordered_map<uint32_t, uint32_t> _parent_inmap;
  std::unordered_map<uint32_t, uint8_t> _depth_inmap;
  std::unordered_map<uint32_t, std::string> _rank_inmap;
  std::unordered_map<uint32_t, std::string> _name_inmap;

public:
  TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy);
  void printTaxonomyRecord();
  T getLowestCommonAncestor(T a, T b);
  uint32_t tiID_from_trID(T trID);
  T trID_from_tiID(uint32_t tiID);
  bool isBasis(T trID);
  decltype(_trID_to_input) &trID_to_input() { return _trID_to_input; }
  decltype(_trID_to_lsroot) &trID_to_lsroot() { return _trID_to_lsroot; }
  decltype(_child_map) &child_map() { return _child_map; }
  decltype(_parent_vec) &parent_vec() { return _parent_vec; }
  decltype(_depth_vec) &depth_vec() { return _depth_vec; }
  decltype(_num_nodes) &get_num_nodes() { return _num_nodes; }
  bool saveTaxonomyRecord(const char *library_dirpath);
};

#endif
