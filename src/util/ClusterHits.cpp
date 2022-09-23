#include "Debug.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "itoa.h"
#include "Util.h"


#ifdef OPENMP
#include <omp.h>
#endif

struct hit{
    std::string alignment;
    double pval;
    unsigned int qPos;
    unsigned int tPos;
    bool qStrand;
    bool tStrand;
};

//std::sort by default sorts the hits by qPos in ascending order
bool operator<(const hit& a, const hit& b)
{
    return a.qPos < b.qPos;
}

//Log of binomial coefficient as defined in combinepvalperset
double LBinCoeff(double* lookup, int M, int k);


double hypergeoDensity(double* lookup, int k, int n, int N, int K){
    if(k < std::max(0,n+K-N) || k > std::min(K,n)){
        return 0.0;
    }
    else{
        return exp(LBinCoeff(lookup, n, k) + LBinCoeff(lookup, N - n, K - k) - LBinCoeff(lookup, N, K));
    }
}

double hypergeoDistribution(double* lookup, int k, int n, int N, int K){
    if(k < std::max(0,n+K-N)){
        return 0.0;
    }
    else if(k > std::min(K,n)){
        return 1.0;
    }
    else{
        double sum = 0;
        for (int i = 0; i < k + 1; i++){
            sum += hypergeoDensity(lookup, i, n, N, K);
        }
        return sum;
    }
}

double logClusterPval(double* lookup, int k, int m, int K, int Nq, int Nt) {
    size_t minKm = (K < m) ? K : m;
    double sum = 0;
    for (size_t Kp = k; Kp <= minKm; Kp++){
        double pHG = hypergeoDensity(lookup, (Kp - 1), (K - 1), (Nq - 1), (m - 1));
        double PHG = hypergeoDistribution(lookup, (k - 2), (Kp - 1), (Nt - 1), (m - 1));
        sum += pHG * (1 - pow(PHG, k));
    }
    if (sum < DBL_EPSILON){
        return log(K - k  + 1) + lookup[k+1] + (k-1) * log(1.0 * K / (Nq * Nt));
    } else {
        return log(1 - pow((1 - sum), (K - k + 1)));
    }
}

double clusterPval_fast(int m, int K, int Nq, int Nt) {
    double power = 2.0 * pow((m - 1), 2) * (K - 1);
    return 1 - pow((1 - 1.0 * K / (Nq * Nt)), power);
}


double logOrderingPval(double* lookup, int k, int m){
    return log(1 - 1.0 * m / k) - m * log(2) - lookup[m+1];
}

int findSpan(const std::vector<hit> &cluster){
    unsigned int iMax = 0;
    unsigned int iMin = INT_MAX;
    unsigned int jMax = 0;
    unsigned int jMin = INT_MAX;
    for (size_t l = 0; l < cluster.size(); l++){
        iMax = (cluster[l].qPos > iMax) ? cluster[l].qPos : iMax;
        iMin = (cluster[l].qPos < iMin) ? cluster[l].qPos : iMin;
        jMax = (cluster[l].tPos > jMax) ? cluster[l].tPos : jMax;
        jMin = (cluster[l].tPos < jMin) ? cluster[l].tPos : jMin;
    }
    int spanI =iMax - iMin + 1;
    int spanJ =jMax - jMin + 1;
    return (spanI > spanJ) ? spanI : spanJ;
}

int findConservedPairs(std::vector<hit> &cluster){
    //re-index hits in clusters with pos in query set in ascending order
    std::sort(cluster.begin(), cluster.end());
    int m = 0;
    for (size_t l = 0; l < cluster.size()-1; l++){
        bool isSameOrder = (cluster[l+1].tPos > cluster[l].tPos);
        bool isSameStrand = (cluster[l].qStrand == cluster[l].tStrand);
        bool isSameStrand2 = (cluster[l+1].qStrand == cluster[l+1].tStrand);
        if((isSameStrand == isSameOrder) && (isSameStrand2 == isSameOrder)){
            m++;
        }
    }
    return m;
}

double computeWeight(double* lookup, int k, int K, int Nq, int Nt){
    double factorNqNtK = 1.0 * Nq * Nt / K;
    double w = (-log(K-k+1) - lookup[k+1] + (k-1) * log(factorNqNtK)) / (-log(K-k+1) + (k-1) * log(2 * factorNqNtK));
    return w;
}

double clusterMatchScore(double* lookup, std::vector<hit> &cluster, int K, int Nq, int Nt, bool useWeight){
    if(cluster.size()== 0 ){
        return 0.0;
    }
    else {
        int span = findSpan(cluster);
        int k = cluster.size();
        double logpClu;
        double logpOrd;
        int m = findConservedPairs(cluster);
        if(k == 2){
            logpClu = log(clusterPval_fast(span, K, Nq, Nt));
            logpOrd = logOrderingPval(lookup, k, m);
        }
        else{
            logpClu = logClusterPval(lookup, k, span, K, Nq, Nt);
            logpOrd = logOrderingPval(lookup, k, m);
        }
        if (useWeight) {
            double w = computeWeight(lookup, k, K, Nq, Nt);
            //full score would be: -log(pClu) - log(pOrd) + log(1 -log(pClu) - log(pOrd))
            return - w * logpClu - (1 - w) * logpOrd;
        } else {
            return - logpClu - logpOrd;
        }
    }
}

//TODO: make this step more efficient
bool isCompatibleCluster(std::vector<hit> &cluster1, std::vector<hit> &cluster2, unsigned int d){
    unsigned int iMax1 = 0;
    unsigned int iMin1 = INT_MAX;
    unsigned int jMax1 = 0;
    unsigned int jMin1 = INT_MAX;
    for (size_t l = 0; l < cluster1.size(); l++){
        iMax1 = (cluster1[l].qPos > iMax1) ? cluster1[l].qPos : iMax1;
        iMin1 = (cluster1[l].qPos < iMin1) ? cluster1[l].qPos : iMin1;
        jMax1 = (cluster1[l].tPos > jMax1) ? cluster1[l].tPos : jMax1;
        jMin1 = (cluster1[l].tPos < jMin1) ? cluster1[l].tPos : jMin1;
    }
    unsigned int iMax2 = 0;
    unsigned int iMin2 = INT_MAX;
    unsigned int jMax2 = 0;
    unsigned int jMin2 = INT_MAX;
    for (size_t l = 0; l < cluster2.size(); l++){
        iMax2 = (cluster2[l].qPos > iMax2) ? cluster2[l].qPos : iMax2;
        iMin2 = (cluster2[l].qPos < iMin2) ? cluster2[l].qPos : iMin2;
        jMax2 = (cluster2[l].tPos > jMax2) ? cluster2[l].tPos : jMax2;
        jMin2 = (cluster2[l].tPos < jMin2) ? cluster2[l].tPos : jMin2;
    }
    return (std::min(jMin1-jMax2,jMin2-jMax1) <= d && std::min(iMin1-iMax2,iMin2-iMax1) <= d);
}


std::vector<hit> groupNodes(const std::vector<std::vector<int>> &nodeList, const std::vector<hit> &matchList, int i, int j, unsigned int d){
    std::vector<hit> cluster1;
    std::vector<hit> cluster2;
    std::vector<hit> cluster;
    //check if one node is empty or two nodes are incompatible, if so return an empty cluster
    if(nodeList[i].size() != 0 && nodeList[j].size() != 0 ){
        for(size_t m = 0; m < nodeList[i].size(); m++){
            cluster1.push_back(matchList[nodeList[i][m]]);
        }
        for(size_t n = 0; n < nodeList[j].size(); n++){
            cluster2.push_back(matchList[nodeList[j][n]]);
        }
        if(isCompatibleCluster(cluster1,cluster2,d)){
            cluster.insert(cluster.begin(), cluster1.begin(), cluster1.end());
            cluster.insert(cluster.end(), cluster2.begin(), cluster2.end());
        }
    }

    return cluster;
}


double multihitPval(double* lookup, std::vector<hit> &cluster, int Nq, double alpha){
    size_t k = 0;
    double r = 0;
    double pvalThreshold = alpha / (Nq + 1);
    double logPvalThr = log(pvalThreshold);
    for (size_t i = 0; i < cluster.size(); ++i) {
        double logPvalue = log(cluster[i].pval);
        if (logPvalue < logPvalThr) {
            k++;
            r -= logPvalue - logPvalThr;
        }
    }
    //multihit edge case
    if (r == 0) {
        return 1.0;
    }
    if (std::isinf(r)) {
        return 0.0;
    }        
    double expMinusR = exp(-r);
    //expMinusR can also underflow
    if (expMinusR == 0) {
        return 0.0;
    }
    double sum = 0;
    for (size_t i = 0; i < k - 1; ++i){
        sum += pow(r,i)/exp(lookup[i+1]);
    }
    return expMinusR * sum;
}

int clusterhits(int argc, const char **argv, const Command &command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> querySizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	querySizeReader.open(DBReader<unsigned int>::NOSORT);

    sizeDBName = std::string(par.db2) + "_set_size";
    sizeDBIndex = std::string(par.db2) + "_set_size.index";
    DBReader<unsigned int> targetSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	targetSizeReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> qlookupReader(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    qlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* qlookup = qlookupReader.getLookup();

    DBReader<unsigned int> tlookupReader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    tlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* tlookup = tlookupReader.getLookup();

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> headerReader(par.hdr3.c_str(), par.hdr3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? resultReader.getDbtype() : Parameters::DBTYPE_OMIT_FILE;
    const bool isDb = par.dbOut;
    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, shouldCompress, dbType);
    writer.open();

    DBWriter headerWriter(par.hdr4.c_str(), par.hdr4Index.c_str(), par.threads,  par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();
    //find the max number of genes/ORFs of Nq&Nt
    //TODO: can we do that faster by scanning the header file?
    unsigned int maxOrfCount = 0;
    for (size_t i = 0; i < querySizeReader.getSize(); ++i) { 
        unsigned int currentCount = Util::fast_atoi<unsigned int>(querySizeReader.getData(i, 0));
        if (currentCount > maxOrfCount) {
            maxOrfCount = currentCount;
        };
    }
    for (size_t i = 0; i < targetSizeReader.getSize(); ++i) { 
        unsigned int currentCount = Util::fast_atoi<unsigned int>(targetSizeReader.getData(i, 0));
        if (currentCount > maxOrfCount) {
            maxOrfCount = currentCount;
        };
    }

    //create a lookup table for all possible log gamma values (n-1)!, n! will be lookup[n + 1]
    double* lGammaLookup = new double[maxOrfCount + 2];
    for (size_t i = 0; i < maxOrfCount + 2; ++i) { 
        lGammaLookup[i] = lgamma(i);
    }


Debug::Progress progress(resultReader.getSize());
unsigned int cluster_idx = 0;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(100 * 1024);
        std::string headerBuffer;
        headerBuffer.reserve(1024);
        std::string header;
        header.reserve(1024);

        const char *entry[255];
        

        const unsigned int d = par.maxGeneGaps; //par.maxGeneGaps, d is the maximum number of genes allowed between two clusters to merge
        const unsigned int cls = par.clusterSize;
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            std::vector<hit> match;

            //from result header file read Nq & Nt
            header = headerReader.getData(i, thread_idx);
            std::vector<std::string>  hdrcolumns = Util::split(header, "\t");
            unsigned int qSet = Util::fast_atoi<size_t>(hdrcolumns[0].c_str());
            unsigned int tSet = Util::fast_atoi<size_t>(hdrcolumns[1].c_str());
            unsigned int Nq = Util::fast_atoi<size_t>(hdrcolumns[2].c_str()); //query set size
            unsigned int Nt = Util::fast_atoi<size_t>(hdrcolumns[3].c_str()); //target set size
            hdrcolumns.clear();
            header.clear();

            //read all hits from a match and extract pos and strand, append them to std::vector<hit> match
            char *data = resultReader.getData(i, thread_idx);
            while (*data != '\0'){
                //read first two column for qid and tid, get entryName from lookups
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int qid = Util::fast_atoi<unsigned int>(entry[0]);
                unsigned int tid = Util::fast_atoi<unsigned int>(entry[1]);

                //from seqid get entryName in lookup, split by "_" and retrieve info
                std::vector<std::string> qcolumns = Util::split(qlookup[qid].entryName, "_");
                std::vector<std::string> tcolumns = Util::split(tlookup[tid].entryName, "_");
                int qStart = Util::fast_atoi<size_t>(qcolumns[qcolumns.size()-2].c_str());
                int qEnd = Util::fast_atoi<size_t>(qcolumns.back().c_str());
                int tStart = Util::fast_atoi<size_t>(tcolumns[tcolumns.size()-2].c_str());
                int tEnd = Util::fast_atoi<size_t>(tcolumns.back().c_str());
                hit tmpHit;

                //pos is the protein index in the genome, strand is determined by start and end coordinates
                tmpHit.alignment = std::string(entry[0], data - entry[0]);
                tmpHit.pval = std::strtod(entry[2], NULL);
                tmpHit.qPos = Util::fast_atoi<size_t>(qcolumns[qcolumns.size()-3].c_str());
                tmpHit.tPos = Util::fast_atoi<size_t>(tcolumns[tcolumns.size()-3].c_str());
                tmpHit.qStrand = (qStart < qEnd) ? 1 : 0;
                tmpHit.tStrand = (tStart < tEnd) ? 1 : 0;

                match.push_back(tmpHit);
            }

            //TODO: total number of hits can be read from the header
            //int K = Util::fast_atoi<int>(hdrcolumns[3].c_str());//total number of hits
            size_t K = match.size(); 

            if(K == 1){
                continue;
            }

            //initiallize distance matrix to be [K][K]
            double** DistMat = new double*[K]; // Rows
            for (size_t i = 0; i < K; i++)
            {
                DistMat[i] = new double[K]; // Columns
            }

            std::vector<int> dmin(K); //index of closest cluster/highest score
            std::vector<std::vector<int>> nodes(K);
            //assign each node with the index of the singleton cluster
            for(size_t n = 0; n < K; n++){
                nodes[n].push_back(n);
            }

            for(size_t i = 0; i < K; i++){
                for(size_t j = 0; j < K; j++){
                    if(i == j){
                        DistMat[i][j] = 0.0; //set score = 0 to self similarities
                    }
                    else{
                        std::vector<hit> tmpCluster = groupNodes(nodes,match,i,j,d);
                        DistMat[i][j] = clusterMatchScore(lGammaLookup,tmpCluster,K,Nq,Nt,par.clusterUseWeight);//score(i,j)
                    }
                    dmin[i] = (DistMat[i][j] > DistMat[i][dmin[i]]) ? j : dmin[i];
                }
            }

            //Score is determined by the maxscore in the first iteration, due to varying Nq & Nt      
            double maxScore = DBL_MAX;
            bool isFirstIter = true;
            double sMin;
            if(par.clusterUseWeight){
                double w = computeWeight(lGammaLookup, 2, K, Nq, Nt);
                sMin  = -w * log(clusterPval_fast(d+1, K, Nq, Nt)) - (1-w)* logOrderingPval(lGammaLookup,2,1);
            } else{
                sMin  = -log(clusterPval_fast(d+1, K, Nq, Nt)) - logOrderingPval(lGammaLookup,2,1);
            }
            while(isFirstIter|| (maxScore >= sMin)){
                size_t i1 = 0;
                size_t i2;
                //find closest pair of clusters (i1,i2), i.e pair of clusters with highest score
                for(size_t s = 0; s < K-1; s++){
                    i1 = 0;
                    for (size_t i = 0; i < K; i++){
                        i1 = (DistMat[i][dmin[i]] > DistMat[i1][dmin[i1]]) ? i : i1;
                    }
                    i2 = dmin[i1];
                }

                maxScore = DistMat[i1][i2];

                if(maxScore != 0){
                    if(isFirstIter){
                        //sMin = maxScore;
                        isFirstIter = false;
                    }
                } else {
                    break;
                    }
                
                //delete node i2 and append all elements in i2 to i1
                for(size_t n = 0; n < nodes[i2].size() ;n++){
                    nodes[i1].push_back(nodes[i2][n]);
                }
                nodes[i2].clear();


                //overwrite row and column i1 with dist[i1,i2][j]
                for(size_t j = 0; j < K; j++){
                    if(i1 == j || i2 == j){
                        DistMat[i1][j] = 0.0; 
                        DistMat[j][i1] = 0.0; 
                    }
                    else{
                        std::vector<hit> tmpCluster = groupNodes(nodes,match,i1,j,d);
                        DistMat[i1][j] = clusterMatchScore(lGammaLookup,tmpCluster,K,Nq,Nt,par.clusterUseWeight);
                        DistMat[j][i1] = DistMat[i1][j];
                    }

                    //set scores in old row and column i2 to 0
                    DistMat[i2][j] = 0.0;
                    DistMat[j][i2] = 0.0;

                    //for row i1 find index dmin[i1] of closest cluster
                    if(j != 0){
                        dmin[i1] = (DistMat[i1][j] > DistMat[i1][dmin[i1]]) ? j : dmin[i1];
                    } else {
                        dmin[i1] = j;
                    }
                    
                    
                    //for every other row j != i1 && j != i2 compare the new score DistMat[j][i1] with DistMat[j][dmin[j]]
                    if(j != i1 && j != i2) {
                        dmin[j] = (DistMat[j][i1] > DistMat[j][dmin[j]]) ? i1 : dmin[j];
                    }

                }

            }

            //print qid & tid of clusters with >=2 hits and cluster + order P-values greater than thresholds
            for(size_t i = 0; i < nodes.size(); i++){
                if(nodes[i].size() >= cls){
                    std::vector<hit> cluster;
                    for(size_t j = 0; j < nodes[i].size(); j++){
                        cluster.push_back(match[nodes[i][j]]);
                    }
                    double pCO = exp(-clusterMatchScore(lGammaLookup, cluster, K, Nq, Nt, par.clusterUseWeight));
                    double pMH = multihitPval(lGammaLookup, cluster, Nq, par.alpha);
                    if(pCO < par.pCluThr && (pMH < par.pMHThr)){ // par.pCluThr,par.pMultiHitThr;
                        headerBuffer.append(SSTR(qSet));
                        headerBuffer.append("\t");
                        headerBuffer.append(SSTR(tSet));
                        headerBuffer.append("\t");
                        headerBuffer.append(SSTR(pCO));
                        headerBuffer.append("\t");
                        headerBuffer.append(SSTR(pMH));
                        headerBuffer.append("\t");
                        headerBuffer.append(SSTR(cluster.size()));
                        headerBuffer.append("\n");
                        for(size_t i = 0; i < cluster.size(); i++){
                            buffer.append(cluster[i].alignment);
                        }
                    unsigned int key = __sync_fetch_and_add(&(cluster_idx), 1);
                    writer.writeData(buffer.c_str(), buffer.length(), key, thread_idx,isDb);
                    headerWriter.writeData(headerBuffer.c_str(), headerBuffer.length(), key, thread_idx);
                    buffer.clear();
                    headerBuffer.clear();
                    }
                }
            }
            // deallocate memory using the delete operator
            for (size_t i = 0; i < K; i++) {
                delete[] DistMat[i];
            }
            delete[] DistMat;
            match.clear();
        }
    }
    writer.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    headerWriter.close(true);
    resultReader.close();
    headerReader.close();


    qlookupReader.close();
    tlookupReader.close();
    querySizeReader.close();
    targetSizeReader.close();
    delete[] lGammaLookup;
    return EXIT_SUCCESS;
}

