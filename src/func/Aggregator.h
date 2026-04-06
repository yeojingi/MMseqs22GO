#ifndef MMSEQS_AGGREGATOR_H
#define MMSEQS_AGGREGATOR_H
#include <map>
#include <vector>
#include <string>
#include "IndexReader.h"
#include "Scorer.h"


class Aggregator {
public:
  Aggregator() {}
  ~Aggregator() {}

  std::vector<std::pair<std::string, float>> aggregateOneQuery(
      const std::map<size_t, float>*,
      IndexReader*,
      unsigned int thread_idx);

  void aggregateAll(
      const EvidenceScore&,
      const IndexReader*,
      const IndexReader*,
      unsigned int thread_idx,
      FILE*
  );
};

#endif
