#ifndef MMSEQS_T1PICKER_H
#define MMSEQS_T1PICKER_H

#include <string>
#include <vector>

#include "Matcher.h"
#include "Sequence.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"

// Greedy per-GO-term resolution ("Greedy" annotation algorithm).
//
// Alignment::run() owns and configures all alignment primitives (realigner,
// query/target Sequence scratch, the target sequence DBReader, the original
// prefilter candidate list for the current query) and passes them in here by
// reference for every query. T1Picker does not do any alignment setup of its
// own and never touches the prefilter/target databases directly -- only the
// GO annotation DB, which it owns.
//
// For a single query, resolveQuery() walks swResults (pass-1, e-value sorted)
// rank by rank. A rank is skipped without doing any alignment work if every
// GO term carried by its target is already resolved -- this decision depends
// on the resolved-set as it stands AFTER earlier (higher-ranked) anchors in
// this same call have already been chopped, searched and voted on, so the
// whole walk has to stay sequential and in-process (no splitting the T1
// selection from the 2nd search into separate batch steps).
//
// Otherwise the target's matched region is chopped into a fragment (T1_i)
// and re-searched against the query's original prefilter candidate list
// (exactly like align2b's single-anchor step), and the resulting hits are
// used to bitscore-vote (DiamondScore-style, Aggregator policy 4 formula) on
// the still-unresolved candidate terms. Terms that clear voteThreshold are
// marked resolved, and the hits that supported them are appended to
// outResults -- this is the real "query-T2" output (plain Matcher::result_t,
// so the normal alignment output path can write it out unchanged).
//
// Two side DBs are additionally written for inspectability, since the final
// outResults alone doesn't say which T1 anchor found which T2. Both are keyed
// by queryDbKey (same key/shape as a normal alignment output -- plain
// Matcher::resultToBuffer rows, one query key, multiple rows), so both can be
// read back with a regular `convertalis <queryDB> <targetDB> <these>` -- per-
// anchor attribution is not kept, just "everything this query's T1s found":
//   - t1Dbw:   queryDbKey -> all of this query's accepted T1_i (query-T1)
//   - t1t2Dbw: queryDbKey -> the union of all anchors' 2nd-search hits (T1-T2)
// Callers that only care about the final query->T2 result never need to
// touch t1Dbw/t1t2Dbw; they exist purely for inspection.
class T1Picker {
public:
    T1Picker(const std::string &targetGoDB, unsigned int threads, bool preload, float voteThreshold,
              DBWriter *t1Dbw, DBWriter *t1t2Dbw);
    ~T1Picker();

    void resolveQuery(
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
    );

private:
    IndexReader *tGoDbr;
    float voteThreshold;

    DBWriter *t1Dbw;
    DBWriter *t1t2Dbw;

    // GO ids carried by a target, looked up by its dbKey. Empty if the target
    // has no GO annotation.
    std::vector<std::string> getGoTerms(unsigned int dbKey, int thread_idx);
};

#endif
