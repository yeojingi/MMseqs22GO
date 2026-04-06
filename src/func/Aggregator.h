#ifndef MMSEQS_AGGREGATOR_H
#define MMSEQS_AGGREGATOR_H
#include <map>
#include <string>
#include <utility>
#include "IndexReader.h"
#include "Scorer.h"


class Aggregator {
public:
  Aggregator() {}
  ~Aggregator() {}

  // returns {GO terms string, score}
  std::pair<std::string, float> aggregateOneQuery(const std::map<size_t, float>* , IndexReader* );
  void aggregateAll(
      const EvidenceScore& ,
      const IndexReader* ,
      const IndexReader* ,
      unsigned int ,
      FILE*
  );
private:
  // Private members and methods 

};



#endif