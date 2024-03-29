#include "taxonomy.h"

TaxonomyInput::TaxonomyInput(const char *taxonomy_dirpath)
{
  std::string main_dir = taxonomy_dirpath;
  std::ifstream nodes_file(main_dir + "/nodes.dmp");
  std::ifstream names_file(main_dir + "/names.dmp");

  if (!nodes_file.good() || !names_file.good()) {
    std::cerr << "Error opening taxonomy files in the given directory" << main_dir << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  uint64_t node_id, parent_id;
  std::string name, trank, junk;

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
          node_id = static_cast<uint64_t>(stoul(token));
          if (node_id == 0) {
            std::cerr << "Attempt to create taxonomy w/ node ID == 0" << std::endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 2: parent_id = static_cast<uint64_t>(stoul(token)); break;
        case 3:
          trank = token;
          finished = true;
          break;
      }
    }
    if (node_id == 1)
      parent_id = 0;
    _parent_map[node_id] = parent_id;
    _rank_map[node_id] = trank;
  }

  while (getline(names_file, line)) {
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
          node_id = static_cast<uint64_t>(stoul(token));
          if (node_id == 0) {
            std::cerr << "Attempt to create taxonomy w/ node ID == 0" << std::endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 2: // name
          name = token;
          break;
        case 4:
          if (token == "scientific name")
            _name_map[node_id] = name;
          finished = true;
          break;
      }
    }
  }

  uint64_t tmp_taxID;
  uint8_t depth;
  for (auto &kv : _parent_map) {
    depth = 0;
    tmp_taxID = kv.first;
    while (_parent_map[tmp_taxID] != 0) {
      tmp_taxID = _parent_map[tmp_taxID];
      depth++;
    }
    if (tmp_taxID == 1)
      depth++;
    _depth_map[kv.first] = depth;
  }

  nodes_file.close();
  names_file.close();
}

uint64_t TaxonomyInput::getParent(uint64_t taxID) { return _parent_map[taxID]; }

std::string TaxonomyInput::getRank(uint64_t taxID) { return _rank_map[taxID]; }

template<typename T>
TaxonomyRecord<T>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy)
{
  _parent_inmap = taxonomy.parent_map();
  _depth_inmap = taxonomy.depth_map();
  _rank_inmap = taxonomy.rank_map();
  _name_inmap = taxonomy.name_map();
  _full_size = _parent_inmap.size();

  std::map<std::string, uint64_t> input_to_taxID;
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string inputn, taxID;
    if (!(std::getline(iss, taxID, '\t') && std::getline(iss, inputn, '\t'))) {
      std::cerr << "Failed to read file for taxon ID to input map!" << std::endl;
      exit(EXIT_FAILURE);
    }
    input_to_taxID[inputn] = (uint64_t)stoul(taxID);
  }
  input_file.close();

  _num_input = input_to_taxID.size();

  _tID_to_taxID[0] = 0;
  _tID_to_taxID[1] = 1;
  _taxID_to_tID[0] = 0;
  _taxID_to_tID[1] = 1;
  T curr_tID = 2;
  uint64_t num_nodes = 1;
  uint64_t parent_taxID = 1;

  for (auto &kv : input_to_taxID) {
    if (_taxID_to_tID.find(kv.second) == _taxID_to_tID.end()) {
      _input_to_tID[kv.first] = curr_tID;
      _tID_to_input[curr_tID].push_back(kv.first);
      _tID_to_taxID[curr_tID] = kv.second;
      _tID_to_rank[curr_tID] = taxonomy.getRank(kv.second);
      _taxID_to_tID[kv.second] = curr_tID;
      parent_taxID = taxonomy.getParent(kv.second);
      while ((_taxID_to_tID.find(parent_taxID) == _taxID_to_tID.end()) && parent_taxID != 1) {
        curr_tID++;
        num_nodes++;
        _tID_to_taxID[curr_tID] = parent_taxID;
        _tID_to_rank[curr_tID] = taxonomy.getRank(parent_taxID);
        _taxID_to_tID[parent_taxID] = curr_tID;
        parent_taxID = taxonomy.getParent(parent_taxID);
      }
      curr_tID++;
      num_nodes++;
    } else {
      _input_to_tID[kv.first] = _taxID_to_tID[kv.second];
      _tID_to_input[_taxID_to_tID[kv.second]].push_back(kv.first);
    }
  }

#ifdef LARGE_TAXONOMY
  if (std::numeric_limits<uint32_t>::max() < num_nodes) {
    std::cerr << "The number of taxonomy nodes " << num_nodes
              << " exceeds supported limit (std::numeric_limits<uint32_t>::max)!" << std::endl;
    exit(EXIT_FAILURE);
  }
#else
  if (std::numeric_limits<uint16_t>::max() < num_nodes) {
    std::cerr << "The number of taxonomy nodes " << num_nodes
              << " exceeds supported limit (std::numeric_limits<uint16_t>::max)!" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  _num_nodes = curr_tID - 1;

  _parent_vec.resize(curr_tID);
  for (auto &kv : _tID_to_taxID) {
    _parent_vec[kv.first] = _taxID_to_tID[taxonomy.getParent(kv.second)];
  }

  for (unsigned int i = 1; i < _parent_vec.size(); ++i) {
    assert((_parent_vec[i] != 0) || i == 1);
    _child_map[_parent_vec[i]].insert(i);
  }

  _depth_vec.resize(curr_tID);
  uint8_t depth;
  for (auto &kv : _tID_to_taxID) {
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

  tT tmp_tID = 0;
  for (unsigned int tID_i = 1; tID_i <= _num_nodes; ++tID_i) {
    tmp_tID = tID_i;
    while (_depth_vec[tmp_tID] > LSR) {
      tmp_tID = _parent_vec[tmp_tID];
    }
    _tID_to_lsroot[tID_i] = tmp_tID;
  }
}

void TaxonomyInput::printTaxonomyInput()
{
  std::cout << "Child : Parent" << std::endl;
  for (auto &kv : _parent_map) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
  std::cout << "Taxon : Rank" << std::endl;
  for (auto &kv : _rank_map) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
}

template<typename T>
void TaxonomyRecord<T>::printTaxonomyRecord()
{
  std::cout << "Taxonomy-record ID : Taxon ID" << std::endl;
  for (auto &kv : _tID_to_taxID) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
  std::cout << "Genome path : Taxonomy-record ID : Taxon ID" << std::endl;
  for (auto &kv : _input_to_tID) {
    std::cout << kv.first << " " << kv.second << " " << _tID_to_taxID[kv.second] << std::endl;
  }
}

template<typename T>
T TaxonomyRecord<T>::getLowestCommonAncestor(T a, T b)
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
T TaxonomyRecord<T>::tID_from_taxID(uint64_t taxID)
{
  return _taxID_to_tID[taxID];
}

template<typename T>
uint64_t TaxonomyRecord<T>::taxID_from_tID(T tID)
{
  return _tID_to_taxID[tID];
}

template<typename T>
bool TaxonomyRecord<T>::isBasis(T tID)
{
  return _tID_to_input.find(tID) != _tID_to_input.end();
}

template<typename T>
bool TaxonomyRecord<T>::saveTaxonomyRecord(const char *library_dirpath)
{
  bool is_ok = true;
  std::string save_filepath(library_dirpath);
  std::vector<std::pair<T, uint64_t>> tIDs_taxIDs(_tID_to_taxID.begin(), _tID_to_taxID.end());
  std::vector<std::pair<uint64_t, uint64_t>> taxIDs_parents(_parent_inmap.begin(), _parent_inmap.end());
  std::vector<std::pair<uint64_t, uint8_t>> taxIDs_depths(_depth_inmap.begin(), _depth_inmap.end());

  FILE *taxonomy_f = IO::open_file((save_filepath + "/taxonomy").c_str(), is_ok, "wb");
  std::fwrite(&_num_nodes, sizeof(T), 1, taxonomy_f);
  std::fwrite(&_full_size, sizeof(uint64_t), 1, taxonomy_f);
  std::fwrite(&_num_input, sizeof(uint64_t), 1, taxonomy_f);
  std::fwrite(_parent_vec.data(), sizeof(T), _num_nodes, taxonomy_f);
  std::fwrite(_depth_vec.data(), sizeof(uint8_t), _num_nodes, taxonomy_f);
  std::fwrite(tIDs_taxIDs.data(), sizeof(std::pair<T, uint64_t>), _num_nodes, taxonomy_f);
  std::fwrite(taxIDs_parents.data(), sizeof(std::pair<uint64_t, uint64_t>), _full_size, taxonomy_f);
  std::fwrite(taxIDs_depths.data(), sizeof(std::pair<uint64_t, uint8_t>), _full_size, taxonomy_f);

  for (auto &kv : _rank_inmap) {
    size_t size_str = kv.second.size();
    std::fwrite(&kv.first, sizeof(uint64_t), 1, taxonomy_f);
    std::fwrite(&size_str, sizeof(size_t), 1, taxonomy_f);
    std::fwrite(&kv.second[0], sizeof(char), size_str, taxonomy_f);
  }
  for (auto &kv : _name_inmap) {
    size_t size_str = kv.second.size();
    std::fwrite(&kv.first, sizeof(uint64_t), 1, taxonomy_f);
    std::fwrite(&size_str, sizeof(size_t), 1, taxonomy_f);
    std::fwrite(&kv.second[0], sizeof(char), size_str, taxonomy_f);
  }

  if (std::ferror(taxonomy_f)) {
    std::puts("I/O error when writing taxonomy-record file to the library.\n");
    is_ok = false;
  }
  std::fclose(taxonomy_f);

  return is_ok;
}

template TaxonomyRecord<uint32_t>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy);

template TaxonomyRecord<uint16_t>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy);

template void TaxonomyRecord<uint32_t>::printTaxonomyRecord();

template void TaxonomyRecord<uint16_t>::printTaxonomyRecord();

template uint32_t TaxonomyRecord<uint32_t>::getLowestCommonAncestor(uint32_t a, uint32_t b);

template uint16_t TaxonomyRecord<uint16_t>::getLowestCommonAncestor(uint16_t a, uint16_t b);

template uint64_t TaxonomyRecord<uint16_t>::taxID_from_tID(uint16_t tID);

template uint64_t TaxonomyRecord<uint32_t>::taxID_from_tID(uint32_t tID);

template uint16_t TaxonomyRecord<uint16_t>::tID_from_taxID(uint64_t taxID);

template uint32_t TaxonomyRecord<uint32_t>::tID_from_taxID(uint64_t taxID);

template bool TaxonomyRecord<uint16_t>::isBasis(uint16_t tID);

template bool TaxonomyRecord<uint32_t>::isBasis(uint32_t tID);

template bool TaxonomyRecord<uint64_t>::saveTaxonomyRecord(const char *library_dirpath);

template bool TaxonomyRecord<uint32_t>::saveTaxonomyRecord(const char *library_dirpath);

template bool TaxonomyRecord<uint16_t>::saveTaxonomyRecord(const char *library_dirpath);
