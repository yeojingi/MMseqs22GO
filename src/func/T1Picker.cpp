#include "T1Picker.h"
#include "Util.h"
#include "Debug.h"

#include <cstring>
#include <cfloat>
#include <climits>
#include <set>
#include <map>
#include <algorithm>

T1Picker::T1Picker(const std::string &targetGoDB, unsigned int threads, bool preload, float voteThreshold,
                   DBWriter *t1Dbw, DBWriter *t1t2Dbw)
        : voteThreshold(voteThreshold), t1Dbw(t1Dbw), t1t2Dbw(t1t2Dbw) {
    tGoDbr = new IndexReader(targetGoDB, threads, IndexReader::GO,
                              preload ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                              DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
}

T1Picker::~T1Picker() {
    delete tGoDbr;
}

// _func entries are one GO id (first space-delimited token) per line.
std::vector<std::string> T1Picker::getGoTerms(unsigned int dbKey, int thread_idx) {
    std::vector<std::string> terms;
    size_t id = tGoDbr->sequenceReader->getId(dbKey);
    if (id == UINT_MAX) {
        return terms;
    }
    const char *data = tGoDbr->sequenceReader->getData(id, thread_idx);
    size_t len = tGoDbr->sequenceReader->getSeqLen(id);
    const char *ptr = data;
    const char *end = data + len;
    while (ptr < end) {
        const char *nl = (const char *) memchr(ptr, '\n', end - ptr);
        size_t lineLen = nl ? (size_t) (nl - ptr) : (size_t) (end - ptr);
        if (lineLen > 0) {
            const char *sp = (const char *) memchr(ptr, ' ', lineLen);
            size_t goLen = sp ? (size_t) (sp - ptr) : lineLen;
            terms.push_back(std::string(ptr, goLen));
        }
        if (!nl) {
            break;
        }
        ptr = nl + 1;
    }
    return terms;
}

void T1Picker::resolveQuery(
        unsigned int queryDbKey, size_t qId,
        char *origQuerySeqData, unsigned int origQueryLen,
        std::vector<Matcher::result_t> &swResults,
        char *origData,
        DBReader<unsigned int> *tdbr,
        Sequence &qSeq, Sequence &dbSeq, Matcher *realigner,
        int covMode, float realignCov, unsigned int realignSwMode, unsigned int lcaSwMode, int seqIdMode,
        bool includeIdentity, bool sameQTDB, unsigned int maxAccept, unsigned int maxReject,
        int thread_idx,
        std::vector<Matcher::result_t> &outResults
) {
    std::set<std::string> resolved;
    std::set<unsigned int> emitted;
    char buffer[1024];
    char fmtBuffer[1024 + 32768 * 4];
    std::string t1OutString;
    std::string t1t2OutString;

    for (size_t c = 0; c < swResults.size(); c++) {
        Matcher::result_t &cand = swResults[c];

        std::vector<std::string> candTermsAll = getGoTerms(cand.dbKey, thread_idx);
        std::vector<std::string> candTerms;
        for (size_t i = 0; i < candTermsAll.size(); i++) {
            if (resolved.find(candTermsAll[i]) == resolved.end()) {
                candTerms.push_back(candTermsAll[i]);
            }
        }
        if (candTerms.empty()) {
            // every GO term this target carries is already resolved -- skip
            // the expensive chop+research entirely
            continue;
        }

        // chop: realign cand against the ORIGINAL query to get an accurate T1_i
        qSeq.mapSequence(qId, queryDbKey, origQuerySeqData, origQueryLen);
        realigner->initQuery(&qSeq);

        size_t candDbId = tdbr->getId(cand.dbKey);
        char *candDbSeqData = tdbr->getData(candDbId, thread_idx);
        if (candDbSeqData == NULL) {
            Debug(Debug::ERROR) << "Sequence " << cand.dbKey << " is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
            EXIT(EXIT_FAILURE);
        }
        dbSeq.mapSequence(candDbId, cand.dbKey, candDbSeqData, tdbr->getSeqLen(candDbId));
        const bool candIsIdentity = (queryDbKey == cand.dbKey && (includeIdentity || sameQTDB)) ? true : false;
        Matcher::result_t topHit = realigner->getSWResult(&dbSeq, INT_MAX, false, covMode, realignCov, FLT_MAX, realignSwMode, seqIdMode, candIsIdentity);

        const unsigned int topHitKey = topHit.dbKey;
        size_t dbId = tdbr->getId(topHitKey);
        char *qSeqData = tdbr->getData(dbId, thread_idx);
        if (qSeqData == NULL) {
            Debug(Debug::ERROR) << "Sequence " << topHitKey << " is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
            EXIT(EXIT_FAILURE);
        }
        const unsigned int fragLen = topHit.dbEndPos - topHit.dbStartPos + 1;
        qSeq.mapSequence(dbId, topHitKey, qSeqData + topHit.dbStartPos, fragLen);
        realigner->initQuery(&qSeq);

        const int topHitScore = topHit.score;

        // query-T1 "result": accumulate this anchor's topHit as a normal
        // alignment row; all of this query's anchors get written together
        // under queryDbKey once the whole walk is done (same shape as a
        // regular alignment output: one key, one row per hit)
        if (t1Dbw != NULL) {
            size_t len = Matcher::resultToBuffer(fmtBuffer, topHit, false);
            t1OutString.append(fmtBuffer, len);
        }

        // 2nd search: T1_i fragment vs the ORIGINAL prefilter candidate list
        // (same maxAccept/maxReject cap as the normal pass-1 loop -- without
        // this, a weak topHitScore lets nearly the whole prefilter candidate
        // list through for every single anchor)
        std::vector<Matcher::result_t> roundHits;
        char *data = origData;
        unsigned int rejected = 0;
        unsigned int passedNum = 0;
        while (*data != '\0' && passedNum < maxAccept && rejected < maxReject) {
            Util::parseKey(data, buffer);
            const unsigned int dbKey = (unsigned int) strtoul(buffer, NULL, 10);
            data = Util::skipLine(data);

            dbId = tdbr->getId(dbKey);
            char *dbSeqData = tdbr->getData(dbId, thread_idx);
            if (dbSeqData == NULL) {
                Debug(Debug::ERROR) << "Sequence " << dbKey << " is required in the prefiltering, but is not contained in the target sequence database!\nPlease check your database.\n";
                EXIT(EXIT_FAILURE);
            }
            dbSeq.mapSequence(dbId, dbKey, dbSeqData, tdbr->getSeqLen(dbId));

            Matcher::result_t res = realigner->getSWResult(&dbSeq, INT_MAX, false, covMode, realignCov, FLT_MAX, lcaSwMode, seqIdMode, false);

            if (res.score >= topHitScore) {
                roundHits.emplace_back(res);
                rejected = 0;
                passedNum++;
            } else {
                rejected++;
            }
        }

        // persist this anchor's full 2nd-search hit list (no backtrace);
        // accumulated across anchors, written once under queryDbKey like t1Dbw
        if (t1t2Dbw != NULL) {
            for (size_t i = 0; i < roundHits.size(); i++) {
                size_t len = Matcher::resultToBuffer(fmtBuffer, roundHits[i], false);
                t1t2OutString.append(fmtBuffer, len);
            }
        }

        // voting (Aggregator policy-4 formula, restricted to this anchor's
        // still-unresolved candidate terms): S(t) = sum(bitscore where t in GO(h)) / sum(bitscore)
        double totalBs = 0.0;
        for (size_t i = 0; i < roundHits.size(); i++) {
            totalBs += roundHits[i].score;
        }
        if (totalBs <= 0.0) {
            continue;
        }

        std::map<std::string, double> goWeightSum;
        for (size_t i = 0; i < roundHits.size(); i++) {
            std::vector<std::string> hitTerms = getGoTerms(roundHits[i].dbKey, thread_idx);
            for (size_t t = 0; t < hitTerms.size(); t++) {
                if (resolved.find(hitTerms[t]) != resolved.end()) {
                    continue;
                }
                goWeightSum[hitTerms[t]] += roundHits[i].score;
            }
        }

        for (size_t t = 0; t < candTerms.size(); t++) {
            std::map<std::string, double>::const_iterator it = goWeightSum.find(candTerms[t]);
            if (it == goWeightSum.end()) {
                continue;
            }
            float support = (float) (it->second / totalBs);
            if (support < voteThreshold) {
                continue;
            }
            resolved.insert(candTerms[t]);

            for (size_t i = 0; i < roundHits.size(); i++) {
                if (emitted.find(roundHits[i].dbKey) != emitted.end()) {
                    continue;
                }
                std::vector<std::string> hitTerms = getGoTerms(roundHits[i].dbKey, thread_idx);
                if (std::find(hitTerms.begin(), hitTerms.end(), candTerms[t]) != hitTerms.end()) {
                    outResults.emplace_back(roundHits[i]);
                    emitted.insert(roundHits[i].dbKey);
                }
            }
        }
    }

    if (t1Dbw != NULL) {
        t1Dbw->writeData(t1OutString.c_str(), t1OutString.length(), queryDbKey, thread_idx);
    }
    if (t1t2Dbw != NULL) {
        t1t2Dbw->writeData(t1t2OutString.c_str(), t1t2OutString.length(), queryDbKey, thread_idx);
    }
}
