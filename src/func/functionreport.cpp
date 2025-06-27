#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"
#include "MappingReader.h"

#include <unordered_map>

// #include "krona_prelude.html.h"

#ifdef OPENMP
#include <omp.h>
#endif


int functionreport(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

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
    return EXIT_SUCCESS;
}

