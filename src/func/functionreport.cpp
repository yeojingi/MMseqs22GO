#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "MappingReader.h"
#include "IndexReader.h"
#include "FuncReader.h"
#include "Matcher.h"
#include "goparser.h"
#include "GeneOntology.h"
#include "Aggregator.h"
#include "Scorer.h"
#include <unordered_set>
#include <unordered_map>

#include "function.html.h"

#ifdef OPENMP
#include <omp.h>
#endif

int functionreport(int argc, const char **argv, const Command &command) {

unsigned int thread_idx = 0;
#ifdef OPENMP
thread_idx = omp_get_thread_num();
#endif
    std::string queryHeaderBuffer;
    queryHeaderBuffer.reserve(1024);

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;

    GeneOntology go(par.db2 + "_func_gog");

    int dbaccessMode = DBReader<unsigned int>::USE_INDEX;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader= &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FILE *resultFP = FileUtil::openAndDelete(par.db4.c_str(), "w");

    IndexReader* tGoDbr = new IndexReader(par.db2, par.threads, IndexReader::GO , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    FuncReader* funcMapping = new FuncReader(db2NoIndexName);

    const int formatMode = par.funcFormatMode;

    EvidenceScore evidenceScore;
    ScoreAlignments(&alnDbr, &evidenceScore, tGoDbr);
    std::cout << "Scored" << std::endl;

    if (formatMode == 1) {
        fwrite(function_html, 1, function_html_len, resultFP);
    }

    Aggregator aggregator;
    aggregator.aggregateAll(
        evidenceScore,
        &qDbrHeader,
        tGoDbr,
        thread_idx,
        resultFP,
        formatMode
    );

    if (formatMode == 1) {
        fprintf(resultFP, "]);</script>");
    }

    std::cout << "Done" << std::endl;

    delete tGoDbr;
    delete funcMapping;
    if (!sameDB) {
        delete tDbr;
        delete tDbrHeader;
    }
    alnDbr.close();
    return EXIT_SUCCESS;
}
