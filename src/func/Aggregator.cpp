#include "Aggregator.h"
#include <map>
#include "Parameters.h"
#include "IndexReader.h"
#include "Scorer.h"
#include "Util.h"
#include "goparser.h"

std::pair<std::string, float> Aggregator::aggregateOneQuery(const std::map<size_t, float>* scores, IndexReader* tGoDbr) {
    Parameters &par = Parameters::getInstance();

    if (par.policy == 0) {
        // Policy 0: best e-value — score is always 1.0
        size_t bestTarget = UINT_MAX;
        float minEvalue = std::numeric_limits<float>::max();
        for (const auto& entry : *scores) {
            if (entry.second < minEvalue) {
                minEvalue = entry.second;
                bestTarget = entry.first;
            }
        }
        if (bestTarget == UINT_MAX) {
            return {"", 1.0f};
        }
        char* goData = tGoDbr->sequenceReader->getData(bestTarget, 0);
        size_t goLen  = tGoDbr->sequenceReader->getSeqLen(bestTarget);
        return {goParser(goData, goLen), 1.0f};

    } else {
        // Policy 1: voting — count votes per targetId, return winner GO + ratio
        std::map<size_t, size_t> voteCounts;
        for (const auto& entry : *scores) {
            voteCounts[entry.first]++;
        }
        size_t bestTarget = UINT_MAX;
        size_t maxVotes = 0;
        for (const auto& entry : voteCounts) {
            if (entry.second > maxVotes) {
                maxVotes = entry.second;
                bestTarget = entry.first;
            }
        }
        if (bestTarget == UINT_MAX) {
            return {"", 0.0f};
        }
        float ratio = static_cast<float>(maxVotes) / static_cast<float>(scores->size());
        char* goData = tGoDbr->sequenceReader->getData(bestTarget, 0);
        size_t goLen  = tGoDbr->sequenceReader->getSeqLen(bestTarget);
        return {goParser(goData, goLen), ratio};
    }
}

void Aggregator::aggregateAll(
    const EvidenceScore& evidenceScores,
    const IndexReader* qDbrHeader,
    const IndexReader* tGoDbr,
    unsigned int thread_idx,
    FILE* resultFP)
{
    for (const auto& [queryId, scores] : evidenceScores) {

        // 1. Query header 읽기
        size_t qHeaderId = qDbrHeader->sequenceReader->getId(queryId);
        const char* qHeader =
            qDbrHeader->sequenceReader->getData(qHeaderId, thread_idx);

        size_t qHeaderLen =
            qDbrHeader->sequenceReader->getSeqLen(qHeaderId);

        std::string queryIdStr =
            Util::parseFastaHeader(qHeader);

        // 2. aggregate target → {GO terms, score}
        auto [goTerms, score] =
            aggregateOneQuery(&scores, (IndexReader*)tGoDbr);

        // 3. target이 없으면 빈 줄 출력
        if (goTerms.empty()) {
            fprintf(resultFP, "%s\t\n", queryIdStr.c_str());
            continue;
        }

        // 4. 출력
        fprintf(resultFP, "%s\t%s\t%.4f\n", queryIdStr.c_str(), goTerms.c_str(), score);
    }
}