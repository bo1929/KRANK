#ifndef _TABLED_H
#define _TABLED_H

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

template void
HTd<uint64_t>::rmCommon(std::vector<uint8_t>& depth_vec);

template void
HTd<uint32_t>::rmCommon(std::vector<uint8_t>& depth_vec);
#endif
