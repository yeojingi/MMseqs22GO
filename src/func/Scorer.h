#ifndef MMSEQS_SCORER_H
#define MMSEQS_SCORER_H
#include "DBReader.h"
#include "IndexReader.h"
#include <map>

class EvidenceScore {
public:
  EvidenceScore() {}
  ~EvidenceScore() {}
  void addScore(unsigned int query_id, size_t target_id, float score) {
    aln_scores[query_id][target_id] = score;
  }

  float getScore(unsigned int query_id, size_t target_id) {
    return aln_scores[query_id][target_id];
  }
  
  auto begin() { return aln_scores.begin(); }
  auto end() { return aln_scores.end(); }
  auto begin() const { return aln_scores.begin(); }
  auto end() const { return aln_scores.end(); }

  void print() const {
    for (const auto& qpair : aln_scores) {
        unsigned int query_id = qpair.first;
        const auto& targetMap = qpair.second;

        std::cout << "Query " << query_id << ":\n";
        for (const auto& tpair : targetMap) {
            unsigned int target_id = tpair.first;
            float score = tpair.second;
            std::cout << "   Target " << target_id 
                      << " → Score: " << score 
                      << "\n";
        }
        std::cout << "\n";
    }
  }
private:
  std::map<unsigned int, std::map<size_t, float> > aln_scores;
};

void ScoreAlignments(DBReader<unsigned int>* , EvidenceScore* , IndexReader* );

#endif