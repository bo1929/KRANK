#include "taxonomy.h"

TaxonomyInput::TaxonomyInput(const char *tax_dirpath)
{
  std::string main_dir = tax_dirpath;
  std::ifstream nodes_file(main_dir + "/nodes.dmp");
  std::ifstream names_file(main_dir + "/names.dmp");

  if (!nodes_file.good() || !names_file.good()) {
    std::cerr << "Error opening taxonomy files in the given directory" << main_dir << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  uint32_t node_id, parent_id;
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
          node_id = static_cast<uint32_t>(stoul(token));
          if (node_id == 0) {
            std::cerr << "Attempt to create taxonomy w/ node ID == 0" << std::endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 2: parent_id = static_cast<uint32_t>(stoul(token)); break;
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
          node_id = static_cast<uint32_t>(stoul(token));
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

  uint32_t tmp_tiID;
  uint8_t depth;
  for (auto &kv : _parent_map) {
    depth = 0;
    tmp_tiID = kv.first;
    while (_parent_map[tmp_tiID] != 0) {
      tmp_tiID = _parent_map[tmp_tiID];
      depth++;
    }
    if (tmp_tiID == 1)
      depth++;
    _depth_map[kv.first] = depth;
  }

  nodes_file.close();
  names_file.close();
}

uint32_t TaxonomyInput::getParent(uint32_t tiID) { return _parent_map[tiID]; }

std::string TaxonomyInput::getRank(uint32_t tiID) { return _rank_map[tiID]; }

template<typename T>
TaxonomyRecord<T>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy)
{
  std::map<std::string, uint32_t> input_to_tiID;
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string inputn, tiID;
    if (!(std::getline(iss, tiID, '\t') && std::getline(iss, inputn, '\t'))) {
      std::cerr << "Failed to read file for taxon ID to input map!" << std::endl;
      exit(EXIT_FAILURE);
    }
    input_to_tiID[inputn] = (uint32_t)stoul(tiID);
  }
  input_file.close();

  _num_input = input_to_tiID.size();

  _trID_to_tiID[0] = 0;
  _trID_to_tiID[1] = 1;
  _tiID_to_trID[0] = 0;
  _tiID_to_trID[1] = 1;
  T curr_trID = 2;
  uint32_t num_nodes = 1;
  uint32_t parent_tiID = 1;

  for (auto &kv : input_to_tiID) {
    if (_tiID_to_trID.find(kv.second) == _tiID_to_trID.end()) {
      _input_to_trID[kv.first] = curr_trID;
      _trID_to_input[curr_trID].push_back(kv.first);
      _trID_to_tiID[curr_trID] = kv.second;
      _trID_to_rank[curr_trID] = taxonomy.getRank(kv.second);
      _tiID_to_trID[kv.second] = curr_trID;
      parent_tiID = taxonomy.getParent(kv.second);
      while ((_tiID_to_trID.find(parent_tiID) == _tiID_to_trID.end()) && parent_tiID != 1) {
        curr_trID++;
        num_nodes++;
        _trID_to_tiID[curr_trID] = parent_tiID;
        _trID_to_rank[curr_trID] = taxonomy.getRank(parent_tiID);
        _tiID_to_trID[parent_tiID] = curr_trID;
        parent_tiID = taxonomy.getParent(parent_tiID);
      }
      curr_trID++;
      num_nodes++;
    } else {
      _input_to_trID[kv.first] = _tiID_to_trID[kv.second];
      _trID_to_input[_tiID_to_trID[kv.second]].push_back(kv.first);
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

  _num_nodes = curr_trID - 1;

  for (auto &kv : _trID_to_tiID) {
    _parent_inmap[kv.second] = taxonomy.parent_map()[kv.second];
    _depth_inmap[kv.second] = taxonomy.depth_map()[kv.second];
    _rank_inmap[kv.second] = taxonomy.rank_map()[kv.second];
    _name_inmap[kv.second] = taxonomy.name_map()[kv.second];
  }

  _parent_vec.resize(curr_trID);
  for (auto &kv : _trID_to_tiID) {
    _parent_vec[kv.first] = _tiID_to_trID[taxonomy.getParent(kv.second)];
  }

  for (unsigned int i = 1; i < _parent_vec.size(); ++i) {
    assert((_parent_vec[i] != 0) || i == 1);
    _child_map[_parent_vec[i]].insert(i);
  }

  _depth_vec.resize(curr_trID);
  uint8_t depth;
  for (auto &kv : _trID_to_tiID) {
    depth = 0;
    curr_trID = kv.first;
    while (_parent_vec[curr_trID] != 0) {
      curr_trID = _parent_vec[curr_trID];
      depth++;
    }
    if (curr_trID == 1)
      depth++;
    _depth_vec[kv.first] = depth;
  }

  tT tmp_trID = 0;
  for (unsigned int trID_i = 1; trID_i <= _num_nodes; ++trID_i) {
    tmp_trID = trID_i;
    while (_depth_vec[tmp_trID] > LSR) {
      tmp_trID = _parent_vec[tmp_trID];
    }
    _trID_to_lsroot[trID_i] = tmp_trID;
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
  for (auto &kv : _trID_to_tiID) {
    std::cout << kv.first << " " << kv.second << std::endl;
  }
  std::cout << "Genome path : Taxonomy-record ID : Taxon ID" << std::endl;
  for (auto &kv : _input_to_trID) {
    std::cout << kv.first << " " << kv.second << " " << _trID_to_tiID[kv.second] << std::endl;
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
T TaxonomyRecord<T>::trID_from_tiID(uint32_t tiID)
{
  return _tiID_to_trID[tiID];
}

template<typename T>
uint32_t TaxonomyRecord<T>::tiID_from_trID(T trID)
{
  return _trID_to_tiID[trID];
}

template<typename T>
bool TaxonomyRecord<T>::isBasis(T trID)
{
  return _trID_to_input.find(trID) != _trID_to_input.end();
}

template<typename T>
bool TaxonomyRecord<T>::saveTaxonomyRecord(const char *library_dirpath)
{
  bool is_ok = true;
  T num_entry = _num_nodes + 1;
  std::string save_filepath(library_dirpath);
  std::vector<std::pair<T, uint32_t>> trIDs_tiIDs(_trID_to_tiID.begin(), _trID_to_tiID.end());
  std::vector<std::pair<uint32_t, uint32_t>> tiIDs_parents(_parent_inmap.begin(), _parent_inmap.end());
  std::vector<std::pair<uint32_t, uint8_t>> tiIDs_depths(_depth_inmap.begin(), _depth_inmap.end());

  FILE *tax_f = IO::open_file((save_filepath + "/taxonomy").c_str(), is_ok, "wb");
  std::fwrite(&_num_nodes, sizeof(T), 1, tax_f);
  std::fwrite(&_num_input, sizeof(uint32_t), 1, tax_f);
  std::fwrite(_parent_vec.data(), sizeof(T), num_entry, tax_f);
  std::fwrite(_depth_vec.data(), sizeof(uint8_t), num_entry, tax_f);
  std::fwrite(trIDs_tiIDs.data(), sizeof(std::pair<T, uint32_t>), num_entry, tax_f);
  std::fwrite(tiIDs_parents.data(), sizeof(std::pair<uint32_t, uint32_t>), num_entry, tax_f);
  std::fwrite(tiIDs_depths.data(), sizeof(std::pair<uint32_t, uint8_t>), num_entry, tax_f);

  for (auto &kv : _rank_inmap) {
    size_t size_str = kv.second.size();
    std::fwrite(&kv.first, sizeof(uint32_t), 1, tax_f);
    std::fwrite(&size_str, sizeof(size_t), 1, tax_f);
    std::fwrite(&kv.second[0], sizeof(char), size_str, tax_f);
  }
  for (auto &kv : _name_inmap) {
    size_t size_str = kv.second.size();
    std::fwrite(&kv.first, sizeof(uint32_t), 1, tax_f);
    std::fwrite(&size_str, sizeof(size_t), 1, tax_f);
    std::fwrite(&kv.second[0], sizeof(char), size_str, tax_f);
  }

  if (std::ferror(tax_f)) {
    std::puts("I/O error when writing taxonomy-record file to the library.\n");
    is_ok = false;
  }
  std::fclose(tax_f);

  return is_ok;
}

template TaxonomyRecord<uint32_t>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy);

template TaxonomyRecord<uint16_t>::TaxonomyRecord(const char *input_filepath, TaxonomyInput taxonomy);

template void TaxonomyRecord<uint32_t>::printTaxonomyRecord();

template void TaxonomyRecord<uint16_t>::printTaxonomyRecord();

template uint32_t TaxonomyRecord<uint32_t>::getLowestCommonAncestor(uint32_t a, uint32_t b);

template uint16_t TaxonomyRecord<uint16_t>::getLowestCommonAncestor(uint16_t a, uint16_t b);

template uint32_t TaxonomyRecord<uint16_t>::tiID_from_trID(uint16_t trID);

template uint32_t TaxonomyRecord<uint32_t>::tiID_from_trID(uint32_t trID);

template uint16_t TaxonomyRecord<uint16_t>::trID_from_tiID(uint32_t tiID);

template uint32_t TaxonomyRecord<uint32_t>::trID_from_tiID(uint32_t tiID);

template bool TaxonomyRecord<uint16_t>::isBasis(uint16_t trID);

template bool TaxonomyRecord<uint32_t>::isBasis(uint32_t trID);

template bool TaxonomyRecord<uint32_t>::saveTaxonomyRecord(const char *library_dirpath);

template bool TaxonomyRecord<uint16_t>::saveTaxonomyRecord(const char *library_dirpath);
