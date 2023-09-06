#include "taxonomy.h"

TaxonomyNCBI::TaxonomyNCBI(const char* nodes_filepath)
{
  std::ifstream nodes_file(nodes_filepath);
  if (!nodes_file.good()) {
    std::cerr << "Error opening " << nodes_filepath << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  uint64_t node_id, parent_id;
  std::string name, rank, junk;

  const std::string delim = "\t|\t";
  while (getline(nodes_file, line)) {
    line.pop_back();
    line.pop_back();
    size_t pos1, pos2;
    pos1 = 0;
    int field_ct = 0;
    bool finished = false;
    while (field_ct++ < 10 && !finished) {
      pos2 = line.find(delim, pos1);
      std::string token;
      if (pos2 == std::string::npos) {
        token = line.substr(pos1);
        finished = true;
      } else {
        token = line.substr(pos1, pos2 - pos1);
        pos1 = pos2 + delim.size();
      }
      switch (field_ct) {
        case 1:
          node_id = (uint64_t)stoul(token);
          if (node_id == 0) {
            std::cerr << "Attempt to create taxonomy w/ node ID == 0" << std::endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 2:
          parent_id = (uint64_t)stoul(token);
          break;
        case 3:
          rank = token;
          finished = true;
          break;
      }
    }
    if (node_id == 1)
      parent_id = 0;
    _parent_map[node_id] = parent_id;
    _rank_map[node_id] = rank;
  }
  nodes_file.close();
}

uint64_t
TaxonomyNCBI::getParent(uint64_t taxID)
{
  return _parent_map[taxID];
}

std::string
TaxonomyNCBI::getRank(uint64_t taxID)
{
  return _rank_map[taxID];
}

template<typename T>
TaxonomyRecord<T>::TaxonomyRecord(const char* input_filepath, TaxonomyNCBI taxonomy)
{
  std::map<std::string, uint64_t> input_to_taxID;
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string genome, taxID;
    if (!(std::getline(iss, genome, '\t') && std::getline(iss, taxID, '\t'))) {
      std::cerr << "Failed to read file for genome to taxon ID map." << std::endl;
      exit(EXIT_FAILURE);
    }
    input_to_taxID[genome] = (uint64_t)stoul(taxID);
  }
  input_file.close();
  _num_input = input_to_taxID.size();

  _tID_to_taxID[0] = 0;
  _tID_to_taxID[1] = 1;
  _taxID_to_tID[0] = 0;
  _taxID_to_tID[1] = 1;
  T curr_tID = 2;
  uint64_t parent_taxID = 1;

  for (auto& kv : input_to_taxID) {
    _input_to_tID[kv.first] = curr_tID;
    _tID_to_input[curr_tID] = kv.first;
    _tID_to_taxID[curr_tID] = kv.second;
    _tID_to_rank[curr_tID] = taxonomy.getRank(kv.second);
    _taxID_to_tID[kv.second] = curr_tID;
    parent_taxID = taxonomy.getParent(kv.second);

    while ((_taxID_to_tID.find(parent_taxID) == _taxID_to_tID.end()) && parent_taxID != 1) {
      curr_tID++;
      _tID_to_taxID[curr_tID] = parent_taxID;
      _tID_to_rank[curr_tID] = taxonomy.getRank(parent_taxID);
      _taxID_to_tID[parent_taxID] = curr_tID;
      parent_taxID = taxonomy.getParent(parent_taxID);
    }

    curr_tID++;
  }

  _num_nodes = curr_tID - 1;

  _parent_vec.resize(curr_tID);
  for (auto& kv : _tID_to_taxID) {
    _parent_vec[kv.first] = _taxID_to_tID[taxonomy.getParent(kv.second)];
  }

  for (unsigned int i = 1; i < _parent_vec.size(); ++i) {
    _child_map[_parent_vec[i]].insert(i);
  }

  _depth_vec.resize(curr_tID);
  uint8_t depth;
  for (auto& kv : _tID_to_taxID) {
    depth = 0;
    curr_tID = kv.first;
    while (_parent_vec[curr_tID] != 0) {
      curr_tID = _parent_vec[curr_tID];
      depth++;
    }
    if (curr_tID == 1)
      depth++;
    _depth_vec[kv.first] = depth;
  }
}

void
TaxonomyNCBI::printTaxonomyNCBI()
{
  std::cout << "Child : Parent" << std::endl;
  for (auto& kv : _parent_map) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
  std::cout << "Taxon : Rank" << std::endl;
  for (auto& kv : _rank_map) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
}

template<typename T>
void
TaxonomyRecord<T>::printTaxonomyRecord()
{
  std::cout << "Taxonomy-record ID : Taxon ID" << std::endl;
  for (auto& kv : _tID_to_taxID) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
  std::cout << "Genomes name : Taxonomy-record ID : Taxon ID" << std::endl;
  for (auto& kv : _input_to_tID) {
    std::cout << kv.first << " " << kv.second << " " << _tID_to_taxID[kv.second] << std::endl;
  }
}

template<typename T>
T
TaxonomyRecord<T>::getLowestCommonAncestor(T a, T b)
{
  if (!a || !b) // LCA(x,0) = LCA(0,x) = x
    return a ? a : b;
  while (a != b) {
    if (_depth_vec[a] < _depth_vec[b])
      b = _parent_vec[b];
    else
      a = _parent_vec[a];
  }
  return a;
}

template<typename T>
T
TaxonomyRecord<T>::changeIDt(uint64_t taxID)
{
  return _taxID_to_tID[taxID];
}

template<typename T>
uint64_t
TaxonomyRecord<T>::changeIDtax(T tID)
{
  return _tID_to_taxID[tID];
}

template<typename T>
bool
TaxonomyRecord<T>::isBasis(T tID)
{
  return _tID_to_input.find(tID) != _tID_to_input.end();
}

template<typename T>
bool
TaxonomyRecord<T>::saveTaxonomyRecord(const char* library_dirpath)
{
  bool is_ok = true;
  std::string save_filepath(library_dirpath);
  std::vector<std::pair<T, uint64_t>> tIDs_taxIDs(_tID_to_taxID.begin(), _tID_to_taxID.end());

  FILE* taxonomyf = IO::open_file((save_filepath + "/taxonomy").c_str(), is_ok, "wb");
  std::fwrite(&_num_input, sizeof(T), 1, taxonomyf);
  std::fwrite(&_num_nodes, sizeof(T), 1, taxonomyf);
  std::fwrite(tIDs_taxIDs.data(), sizeof(std::pair<T, uint64_t>), _num_nodes, taxonomyf);
  std::fwrite(_parent_vec.data(), sizeof(T), _num_nodes, taxonomyf);
  std::fwrite(_depth_vec.data(), sizeof(uint8_t), _num_nodes, taxonomyf);

  if (std::ferror(taxonomyf)) {
    std::puts("I/O error when writing taxonomy-record file to the library.\n");
    is_ok = false;
  }
  std::fclose(taxonomyf);

  return is_ok;
}

template TaxonomyRecord<uint32_t>::TaxonomyRecord(const char* input_filepath, TaxonomyNCBI taxonomy);

template TaxonomyRecord<uint16_t>::TaxonomyRecord(const char* input_filepath, TaxonomyNCBI taxonomy);

template void
TaxonomyRecord<uint32_t>::printTaxonomyRecord();

template void
TaxonomyRecord<uint16_t>::printTaxonomyRecord();

template uint32_t
TaxonomyRecord<uint32_t>::getLowestCommonAncestor(uint32_t a, uint32_t b);

template uint16_t
TaxonomyRecord<uint16_t>::getLowestCommonAncestor(uint16_t a, uint16_t b);

template uint64_t
TaxonomyRecord<uint16_t>::changeIDtax(uint16_t tID);

template uint64_t
TaxonomyRecord<uint32_t>::changeIDtax(uint32_t tID);

template uint16_t
TaxonomyRecord<uint16_t>::changeIDt(uint64_t taxID);

template uint32_t
TaxonomyRecord<uint32_t>::changeIDt(uint64_t taxID);

template bool
TaxonomyRecord<uint16_t>::isBasis(uint16_t tID);

template bool
TaxonomyRecord<uint32_t>::isBasis(uint32_t tID);

template bool
TaxonomyRecord<uint64_t>::saveTaxonomyRecord(const char* library_dirpath);

template bool
TaxonomyRecord<uint32_t>::saveTaxonomyRecord(const char* library_dirpath);

template bool
TaxonomyRecord<uint16_t>::saveTaxonomyRecord(const char* library_dirpath);
