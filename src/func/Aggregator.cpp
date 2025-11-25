#include "Aggregator.h"
#include <map>
#include "IndexReader.h"
#include "Scorer.h"
#include "Util.h"
#include "goparser.h"

size_t Aggregator::aggregateOneQuery(const std::map<size_t, float>* scores, IndexReader* tGoDbr) {
    size_t aggregatedResult = UINT_MAX;

    // Choose the smallest evalue among all scores as an example
    float minEvalue = std::numeric_limits<float>::max();
    for (const auto& entry : *scores) {
        size_t targetId = entry.first;
        float evalue = entry.second;
        if (evalue < minEvalue) {
            minEvalue = evalue;
            aggregatedResult = targetId;
        }
    }

    return aggregatedResult;
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

        // 2. aggregate target
        size_t aggregatedTargetId =
            aggregateOneQuery(&scores, (IndexReader*)tGoDbr);

        // 3. target이 없으면 빈 줄 출력
        if (aggregatedTargetId == UINT_MAX) {
            fprintf(resultFP, "%s\t\n", queryIdStr.c_str());
            continue;
        }

        // 4. GO terms 읽기
        char* tGoData =
            tGoDbr->sequenceReader->getData(aggregatedTargetId, thread_idx);

        size_t tGoLen =
            tGoDbr->sequenceReader->getSeqLen(aggregatedTargetId);

        std::string goTerms = goParser(tGoData, tGoLen);

        // 5. 출력
        fprintf(resultFP, "%s\t%s\n", queryIdStr.c_str(), goTerms.c_str());
    }
}