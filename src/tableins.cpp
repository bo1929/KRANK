#ifndef _TABLED_H
#define _TABLED_H

template uint64_t
StreamIM<uint32_t>::processInput(uint64_t rbatch_size);

template uint64_t
StreamIM<uint64_t>::processInput(uint64_t rbatch_size);

template bool
StreamIM<uint64_t>::save(const char* filepath);

template bool
StreamIM<uint32_t>::save(const char* filepath);

template bool
StreamIM<uint64_t>::load(const char* filepath);

template bool
StreamIM<uint32_t>::load(const char* filepath);

template void
StreamIM<uint32_t>::clear();

template void
StreamIM<uint64_t>::clear();

template uint64_t
StreamIM<uint64_t>::getBatch(vvec<uint64_t>& batch_table, uint32_t tbatch_size);

template uint64_t
StreamIM<uint32_t>::getBatch(vvec<uint32_t>& batch_table, uint32_t tbatch_size);

template void
StreamOD<uint64_t>::openStream();

template void
StreamOD<uint32_t>::openStream();

template uint64_t
StreamOD<uint64_t>::getBatch(vvec<uint64_t>& batch_table, uint32_t tbatch_size);

template uint64_t
StreamOD<uint32_t>::getBatch(vvec<uint32_t>& batch_table, uint32_t tbatch_size);

template std::unordered_map<uint8_t, uint64_t>
StreamIM<uint64_t>::histRowSizes();

template std::unordered_map<uint8_t, uint64_t>
StreamIM<uint32_t>::histRowSizes();

template void
HTd<uint32_t>::makeUnique(bool update_size);

template void
HTd<uint64_t>::makeUnique(bool update_size);

template void
HTs<uint32_t>::makeUnique(bool update_size);

template void
HTs<uint64_t>::makeUnique(bool update_size);

template void
HTd<uint32_t>::trimColumns(uint8_t b);

template void
HTd<uint64_t>::trimColumns(uint8_t b);

template void
HTd<uint32_t>::pruneColumns(uint8_t b_max);

template void
HTd<uint64_t>::pruneColumns(uint8_t b_max);

template void
HTd<uint64_t>::sortColumns();

template void
HTd<uint32_t>::sortColumns();

template void
HTs<uint64_t>::sortColumns();

template void
HTs<uint32_t>::sortColumns();

template bool
HTd<uint64_t>::areColumnsSorted();

template bool
HTd<uint32_t>::areColumnsSorted();

template bool
HTs<uint64_t>::areColumnsSorted();

template bool
HTs<uint32_t>::areColumnsSorted();

template void
HTd<uint32_t>::clearRows();

template void
HTd<uint64_t>::clearRows();

template void
HTs<uint32_t>::clearRows();

template void
HTs<uint64_t>::clearRows();

template void
HTd<uint64_t>::updateSize();

template void
HTd<uint32_t>::updateSize();

template void
HTs<uint64_t>::updateSize();

template void
HTs<uint32_t>::updateSize();

template void
HTd<uint64_t>::initBasis(tT tID);

template void
HTd<uint32_t>::initBasis(tT tID);

template void
HTd<uint64_t>::unionRows(HTd<uint64_t>& sibling, bool update_size);

template void
HTd<uint32_t>::unionRows(HTd<uint32_t>& sibling, bool update_size);

template void
HTs<uint64_t>::unionRows(HTs<uint64_t>& sibling, bool update_size);

template void
HTs<uint32_t>::unionRows(HTs<uint32_t>& sibling, bool update_size);

template void
HTd<uint64_t>::mergeRows(HTd<uint64_t>& sibling, bool update_size);

template void
HTd<uint32_t>::mergeRows(HTd<uint32_t>& sibling, bool update_size);

template void
HTs<uint64_t>::mergeRows(HTs<uint64_t>& sibling, bool update_size);

template void
HTs<uint32_t>::mergeRows(HTs<uint32_t>& sibling, bool update_size);

template void
HTd<uint64_t>::shrinkHT(uint64_t num_rm, uint8_t b_max);

template void
HTd<uint32_t>::shrinkHT(uint64_t num_rm, uint8_t b_max);

template void
HTs<uint64_t>::shrinkHT(uint64_t num_rm);

template void
HTs<uint32_t>::shrinkHT(uint64_t num_rm);

template std::unordered_map<uint8_t, uint64_t>
HTd<uint64_t>::histRowSizes();

template std::unordered_map<uint8_t, uint64_t>
HTd<uint32_t>::histRowSizes();

template std::unordered_map<uint8_t, uint64_t>
HTs<uint64_t>::histRowSizes();

template std::unordered_map<uint8_t, uint64_t>
HTs<uint32_t>::histRowSizes();

template void
HTd<uint32_t>::convertHTs(HTs<uint32_t>* table);

template void
HTd<uint64_t>::convertHTs(HTs<uint64_t>* table);

template void
HTd<uint32_t>::updateLCA();

template void
HTd<uint64_t>::updateLCA();

template void
HTs<uint64_t>::updateLCA();

template void
HTs<uint32_t>::updateLCA();

#endif
