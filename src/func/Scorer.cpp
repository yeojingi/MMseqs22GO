#include "Scorer.h"
#include "Debug.h"
#include "Matcher.h"
#include "IndexReader.h"

void ScoreAlignments(DBReader<unsigned int>* alnDbr, EvidenceScore* evidenceScore, IndexReader* tGoDbr) {
    const int thread_idx = 0;
    Debug::Progress progress(alnDbr->getSize());

    for (size_t i = 0; i < alnDbr->getSize(); i++) {
        progress.updateProgress();
        const unsigned int queryKey = alnDbr->getDbKey(i);

        char *data = alnDbr->getData(i, thread_idx);

        while (*data != '\0') {
            Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
            data = Util::skipLine(data);
            size_t tGoId = tGoDbr->sequenceReader->getId(res.dbKey);
            if (tGoId == UINT_MAX) continue; // GO annotation 없는 target 제외
            evidenceScore->addScore(queryKey, tGoId, res.eval);
        }
    }
}