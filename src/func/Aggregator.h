#ifndef MMSEQS_AGGREGATOR_H
#define MMSEQS_AGGREGATOR_H
#include <map>
#include "IndexReader.h"
#include "Scorer.h"


class Aggregator {
public:
  Aggregator() {}
  ~Aggregator() {}

  size_t aggregateOneQuery(const std::map<size_t, float>* , IndexReader* );
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