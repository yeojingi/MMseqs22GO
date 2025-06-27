#include "GeneOntology.h"
// #include "FileUtil.h"
// #include "MathUtil.h"
#include "Debug.h"
// #include "Util.h"
// #include "sys/mman.h"

#include <fstream>
// #include <algorithm>
// #include <cassert>

// const int NcbiTaxonomy::SERIALIZATION_VERSION = 2;

// int **makeMatrix(size_t maxNodes) {
//     size_t dimension = maxNodes * 2;
//     int **M = new int*[dimension];
//     int k = (int)(MathUtil::flog2(dimension)) + 1;
//     M[0] = new int[dimension * k]();
//     for(size_t i = 1; i < dimension; i++) {
//         M[i] = M[i-1] + k;
//     }

//     return M;
// }

// void deleteMatrix(int** M) {
//     delete[] M[0];
//     delete[] M;
// }

GeneOntology::GeneOntology(const std::string &fileName) {
    // block = new StringBlock<unsigned int>();
    std::vector<GoNode> tmpNodes;
    // loadNodes(tmpNodes, nodesFile);
    // loadMerged(mergedFile);
    // loadNames(tmpNodes, namesFile);
    loadGo(tmpNodes, fileName);

    // maxNodes = tmpNodes.size();
    // taxonNodes = new TaxonNode[maxNodes];
    // std::copy(tmpNodes.begin(), tmpNodes.end(), taxonNodes);

    // std::vector<int> tmpE;
    // tmpE.reserve(maxNodes * 2);

    // std::vector<int> tmpL;
    // tmpL.reserve(maxNodes * 2);

    // H = new int[maxNodes];
    // std::fill(H, H + maxNodes, 0);

    // std::vector<std::vector<TaxID>> children(tmpNodes.size());
    // for (std::vector<TaxonNode>::const_iterator it = tmpNodes.begin(); it != tmpNodes.end(); ++it) {
    //     if (it->parentTaxId != it->taxId) {
    //         children[nodeId(it->parentTaxId)].push_back(it->taxId);
    //     }
    // }

    // elh(children, 1, 0, tmpE, tmpL);
    // tmpE.resize(maxNodes * 2, 0);
    // tmpL.resize(maxNodes * 2, 0);

    // E = new int[maxNodes * 2];
    // std::copy(tmpE.begin(), tmpE.end(), E);
    // L = new int[maxNodes * 2];
    // std::copy(tmpL.begin(), tmpL.end(), L);

    // M = makeMatrix(maxNodes);
    // computeSparseTable();

    // mmapData = NULL;
    // mmapSize = 0;
}

GeneOntology::~GeneOntology() {
    // if (externalData) {
    //     delete[] M;
    // } else {
    //     delete[] taxonNodes;
    //     delete[] H;
    //     delete[] D;
    //     delete[] E;
    //     delete[] L;
    //     deleteMatrix(M);
    // }
    // delete block;
    // if (mmapData != NULL) {
    //     munmap(mmapData, mmapSize);
    // }
}

void GeneOntology::loadGo(std::vector<GoNode> tmpNodes, const std::string &goFile) {
    Debug(Debug::INFO) << "Loading a Gene Ontology file ...";
    std::ifstream ss(goFile);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << goFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
    //     if (line.find("scientific name") == std::string::npos) {
    //         continue;
    //     }

    //     std::pair<int, std::string> entry = parseName(line);
    //     if (!nodeExists(entry.first)) {
    //         Debug(Debug::ERROR) << "loadNames: Taxon " << entry.first << " not present in nodes file!\n";
    //         EXIT(EXIT_FAILURE);
    //     }
    //     tmpNodes[nodeId(entry.first)].nameIdx = block->append(entry.second.c_str(), entry.second.size());
    }
    Debug(Debug::INFO) << " Done\n";
}

GeneOntology * GeneOntology::openFuncDb(const std::string &database){
    std::string binFile = database + "_taxonomy";
    // if (FileUtil::fileExists(binFile.c_str())) {
    //     FILE* handle = fopen(binFile.c_str(), "r");
    //     struct stat sb;
    //     if (fstat(fileno(handle), &sb) < 0) {
    //         Debug(Debug::ERROR) << "Failed to fstat file " << binFile << "\n";
    //         EXIT(EXIT_FAILURE);
    //     }
    //     char* data = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fileno(handle), 0);
    //     if (data == MAP_FAILED){
    //         Debug(Debug::ERROR) << "Failed to mmap file " << binFile << " with error " << errno << "\n";
    //         EXIT(EXIT_FAILURE);
    //     }
    //     fclose(handle);
    //     NcbiTaxonomy* t = NcbiTaxonomy::unserialize(data);
    //     if (t != NULL) {
    //         t->mmapData = data;
    //         t->mmapSize = sb.st_size;
    //         return t;
    //     } else {
    //         Debug(Debug::WARNING) << "Outdated taxonomy information, please recreate with createtaxdb.\n";
    //     }
    // }
    // Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    // std::string nodesFile = database + "_nodes.dmp";
    // std::string namesFile = database + "_names.dmp";
    // std::string mergedFile = database + "_merged.dmp";
    // if (FileUtil::fileExists(nodesFile.c_str())
    //     && FileUtil::fileExists(namesFile.c_str())
    //     && FileUtil::fileExists(mergedFile.c_str())) {
    // } else if (FileUtil::fileExists("nodes.dmp")
    //            && FileUtil::fileExists("names.dmp")
    //            && FileUtil::fileExists("merged.dmp")) {
    //     nodesFile = "nodes.dmp";
    //     namesFile = "names.dmp";
    //     mergedFile = "merged.dmp";
    // } else {
    //     Debug(Debug::ERROR) << "names.dmp, nodes.dmp, merged.dmp from NCBI taxdump could not be found!\n";
    //     EXIT(EXIT_FAILURE);
    // }
    // return new NcbiTaxonomy(namesFile, nodesFile, mergedFile);
    return new GeneOntology(binFile);
}