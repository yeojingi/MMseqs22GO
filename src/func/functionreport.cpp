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

#define NO_ALIGNMNET 0

int functionreport(int argc, const char **argv, const Command &command) {

unsigned int thread_idx = 0;
#ifdef OPENMP
thread_idx = omp_get_thread_num();
#endif        
    std::string queryHeaderBuffer;
    queryHeaderBuffer.reserve(1024);

    Parameters &par = Parameters::getInstance();
    bool needFuncMapping = true;
    bool needBacktrace = false;

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;

    par.parseParameters(argc, argv, command, true, 0, 0);
    
    GeneOntology go(par.db2 + "_func_gog");
    // std::cout << go.getLineage(1) << std::endl;
    // const GoNode* node = go.getGo(1); // test

    // if (node) {
    //     std::cout << node->goName << std::endl;
    // } else {
    //     std::cout << "GO ID 1 not found\n";
    // }

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
    // fprintf(resultFP, "Hello\n");

    FuncReader* funcMapping = NULL;
    IndexReader* tGoDbr = NULL;
    
    if (needFuncMapping) {
        // IndexReader fDbr(par.db1, par.threads,  IndexReader::GO, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tGoDbr = new IndexReader(par.db2, par.threads, IndexReader::GO , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        funcMapping = new FuncReader(db2NoIndexName);
    }

    EvidenceScore evidenceScore;
    ScoreAlignments(&alnDbr, &evidenceScore, tGoDbr);
    std::cout << "Scored" << std::endl;

    // evidenceScore.print();

    // Iterate evidence Score and run aggregateOneQuery and store the GO results for each query
    Aggregator aggregator;
    std::cout << "Aggregator" << std::endl;
    aggregator.aggregateAll(
        evidenceScore,
        &qDbrHeader,
        tGoDbr,
        thread_idx,
        resultFP
    );

    // Debug::Progress progress(alnDbr.getSize());
    // if (needFuncMapping) {
    //     // IndexReader fDbr(par.db1, par.threads,  IndexReader::GO, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    //     tGoDbr = new IndexReader(par.db2, par.threads, IndexReader::GO , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    //     std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    //     funcMapping = new FuncReader(db2NoIndexName);
    // }
    // for (size_t i = 0; i < alnDbr.getSize(); i++) {
    //     progress.updateProgress();
    //     // fprintf(resultFP, "%d\n", i);
    //     const unsigned int queryKey = alnDbr.getDbKey(i);


    //     size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
    //     const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
    //     size_t qHeaderLen = qDbrHeader.sequenceReader->getSeqLen(qHeaderId);
    //     std::string queryId = Util::parseFastaHeader(qHeader);
    //     queryHeaderBuffer.assign(qHeader, qHeaderLen);
    //     qHeader = (char*) queryHeaderBuffer.c_str();

    //     char *data = alnDbr.getData(i, thread_idx);
    //     int tGoId_max = UINT_MAX;
    //     char *tGo_max;
    //     size_t tGoLen_max;
    //     double eval_min=100;
    //     while (*data != '\0') {
    //         Matcher::result_t res = Matcher::parseAlignmentRecord(data, true);
    //         data = Util::skipLine(data);

    //         if (res.backtrace.empty() && needBacktrace == true) {
    //             Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
    //                                     "Command: mmseqs align " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
    //             EXIT(EXIT_FAILURE);
    //         }

    //         if (res.eval < eval_min) {
    //             eval_min = res.eval;

    //             // Debug(Debug::ERROR) << res.dbKey << "\n";
    //             size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
    //             size_t tGoId = tGoDbr->sequenceReader->getId(res.dbKey);
    //             if (tGoId == UINT_MAX) {
    //                 continue;
    //             }
    //             // std::cout << tGoId << std::endl;
    //             tGo_max = tGoDbr->sequenceReader->getData(tGoId, thread_idx);
    //             tGoLen_max = tGoDbr->sequenceReader->getSeqLen(tGoId);
    //             tGoId_max = tGoId;
    //         }
    //     }
    //     if (eval_min != 100) {
    //         // fprintf(resultFP, "%s\t%s\n", queryId.c_str(), goParser(tGo_max, tGoLen_max).c_str());
    //         if (tGoId_max == UINT_MAX) {
    //             continue;
    //         }
    //         fprintf(resultFP, "%s\t%s\t%.4e\n", queryId.c_str(), goParser(tGo_max, tGoLen_max).c_str(), eval_min);
    //     }
    // }
    std::cout << "Done" << std::endl;
//     NcbiTaxonomy *taxDB = NcbiTaxonomy::openTaxonomy(par.db1);
//     // allow reading any kind of sequence database
//     const int readerDbType = FileUtil::parseDbType(par.db2.c_str());
//     const bool isSequenceDB = Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_HMM_PROFILE)
//                              || Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_AMINO_ACIDS)
//                              || Parameters::isEqualDbtype(readerDbType, Parameters::DBTYPE_NUCLEOTIDES);
//     int dataMode = DBReader<unsigned int>::USE_INDEX;
//     if (isSequenceDB == false) {
//         dataMode |= DBReader<unsigned int>::USE_DATA;
//     }
//     DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, dataMode);
//     reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

//     // support reading both LCA databases and result databases (e.g. alignment)
//     const bool isTaxonomyInput = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_TAXONOMICAL_RESULT);
//     MappingReader* mapping = NULL;
//     if (isTaxonomyInput == false) {
//         mapping = new MappingReader(par.db1);
//     }

//     FILE *resultFP = FileUtil::openAndDelete(par.db3.c_str(), "w");

//     std::unordered_map<TaxID, unsigned int> taxCounts;
//     Debug::Progress progress(reader.getSize());
// #pragma omp parallel
//     {
//         unsigned int thread_idx = 0;
// #ifdef OPENMP
//         thread_idx = (unsigned int) omp_get_thread_num();
// #endif

//         std::unordered_map<TaxID, unsigned int> localTaxCounts;
// #pragma omp for schedule(dynamic, 10)
//         for (size_t i = 0; i < reader.getSize(); ++i) {
//             progress.updateProgress();

//             if (isSequenceDB == true) {
//                 unsigned int taxon = mapping->lookup(reader.getDbKey(i));
//                 if (taxon != 0) {
//                     ++localTaxCounts[taxon];
//                 }
//                 continue;
//             }

//             char *data = reader.getData(i, thread_idx);
//             while (*data != '\0') {
//                 if (isTaxonomyInput) {
//                     TaxID taxon = Util::fast_atoi<int>(data);
//                     ++localTaxCounts[taxon];
//                 } else {
//                     // match dbKey to its taxon based on mapping
//                     unsigned int taxon = mapping->lookup(Util::fast_atoi<unsigned int>(data));
//                     if (taxon != 0) {
//                         ++localTaxCounts[taxon];
//                     }
//                 }
//                 data = Util::skipLine(data);
//             }
//         }

//         // merge maps again
// #pragma omp critical
//         for (std::unordered_map<TaxID, unsigned int>::const_iterator it = localTaxCounts.cbegin(); it != localTaxCounts.cend(); ++it) {
//             if (taxCounts[it->first]) {
//                 taxCounts[it->first] += it->second;
//             } else {
//                 taxCounts[it->first] = it->second;
//             }
//         }
//     }
//     Debug(Debug::INFO) << "Found " << taxCounts.size() << " different taxa for " << reader.getSize() << " different reads\n";
//     unsigned int unknownCnt = (taxCounts.find(0) != taxCounts.end()) ? taxCounts.at(0) : 0;
//     Debug(Debug::INFO) << unknownCnt << " reads are unclassified\n";
//     const size_t entryCount = reader.getSize();
//     reader.close();

//     std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxDB->getCladeCounts(taxCounts);
//     if (par.reportMode == 0) {
//         taxReport(resultFP, *taxDB, cladeCounts, entryCount);
//     } else {
//         fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), resultFP);
//         fprintf(resultFP, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", entryCount);
//         kronaReport(resultFP, *taxDB, cladeCounts, entryCount);
//         fprintf(resultFP, "</node></krona></div></body></html>");
//     }
//     delete taxDB;
//     if (fclose(resultFP) != 0) {
//         Debug(Debug::ERROR) << "Cannot close file " << par.db3 << "\n";
//         return EXIT_FAILURE;
//     }
    // delete evidenceScore;
    // delete(evidenceScore);
    delete tGoDbr;
    delete funcMapping;
    if (!sameDB) {
        delete tDbr;
        delete tDbrHeader;
    } 
    alnDbr.close();
    return EXIT_SUCCESS;
}

