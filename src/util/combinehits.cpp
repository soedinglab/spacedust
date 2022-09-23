#include "Debug.h"
#include "LocalParameters.h"
#include "Aggregation.h"
#include "itoa.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif


double LBinCoeff(double* lookup, int M, int k);

// Precompute coefficients logB[i] = log(B[i])
void precomputeLogB(const unsigned int orfCount, const double pvalThreshold, double* lGammaLookup, double *logB);


class PvalueAggregator : public Aggregation {
public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, float alpha, unsigned int threads, unsigned int compressed, int aggregationMode, bool filterSelfMatch) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), alpha(alpha), aggregationMode(aggregationMode), filterSelfMatch(filterSelfMatch) {

        std::string sizeDBName = queryDbName + "_set_size";
        std::string sizeDBIndex = queryDbName + "_set_size.index";
        querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        sizeDBName = targetDbName + "_set_size";
        sizeDBIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);

        unsigned int maxOrfCount = 0;
        for (size_t i = 0; i < querySizeReader->getSize(); ++i) { 
            unsigned int currentCount = Util::fast_atoi<unsigned int>(querySizeReader->getData(i, 0));
            if (currentCount > maxOrfCount) {
                maxOrfCount = currentCount;
            };
        }

        lGammaLookup = new double[maxOrfCount + 2];
        for (size_t i = 0; i < maxOrfCount + 2; ++i) { 
            lGammaLookup[i] = lgamma(i);
        }

        logBiLookup = new double*[threads];
        for (size_t i = 0; i < threads; ++i) {
            logBiLookup[i] = new double[maxOrfCount];
        }
    }

    ~PvalueAggregator() {
        for (size_t i = 0; i < threads; ++i) {
            delete[] logBiLookup[i];
        }
        delete[] logBiLookup;

        delete[] lGammaLookup;

        targetSizeReader->close();
        delete targetSizeReader;

        querySizeReader->close();
        delete querySizeReader;
    }

    void prepareInput(unsigned int querySetKey, unsigned int thread_idx) {
        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        precomputeLogB(orfCount, alpha/(orfCount + 1), lGammaLookup, logBiLookup[thread_idx]);
    }

    //Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey,
                               unsigned int targetSetKey, unsigned int thread_idx) {
        
        const size_t numTargetSets = targetSizeReader->getSize();  

        std::string buffer;
        buffer.reserve(10 * 1024);
        std::vector<std::vector<std::string>*> entries;

        if(filterSelfMatch && querySetKey == targetSetKey) {
            return "";
        }

        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx)); 
        unsigned int targetOrfCount = Util::fast_atoi<unsigned int>(targetSizeReader->getDataByDBKey(targetSetKey, thread_idx));

        buffer.append(SSTR(querySetKey));
        buffer.append("\t");
        buffer.append(SSTR(targetSetKey));
        buffer.append("\t");
        buffer.append(SSTR(orfCount));
        buffer.append("\t");
        buffer.append(SSTR(targetOrfCount));
        buffer.append("\t");

        //0) multihit P-values
        if(aggregationMode == Parameters::AGGREGATION_MODE_MULTIHIT){
            double pvalThreshold = alpha / (orfCount + 1);

            //multihit edge case p0 = 0
            if (pvalThreshold == 0.0) {
                return "";
            }

            size_t k = 0;
            double r = 0;
            const double logPvalThr = log(pvalThreshold);
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][2].c_str(), NULL);
                if (logPvalue < logPvalThr) {
                    k++;
                    r -= logPvalue - logPvalThr;
                    entries.push_back(&dataToAggregate[i]);
                }
            }
            //multihit edge case r = 0
            if (r == 0 || k == 0) {
                return "";
            }

            buffer.append(SSTR(k));
            buffer.append("\t");

            if (std::isinf(r)) {
                buffer.append(SSTR(0.0));
                goto next;
            }

            const double expMinusR = exp(-r);
            //expMinusR can also underflow
            if (expMinusR == 0) {
                buffer.append(SSTR(0.0));
                goto next;
            }

            //multihit edge case p0 = 1
            if (pvalThreshold == 1.0) {
                //TODO: add parameter to cutoff by cEval?
                buffer.append(SSTR(expMinusR * numTargetSets)); 
                goto next;
            }


            double truncatedFisherPval = 0;
            const double logR = log(r);
            for (size_t i = 0; i < orfCount; ++i) { 
                truncatedFisherPval += exp(i*logR - lGammaLookup[i+1] + logBiLookup[thread_idx][i]);
            }
            double updatedPval = expMinusR * truncatedFisherPval;       
            double updatedEval = updatedPval * numTargetSets;
            buffer.append(SSTR(updatedEval));
        }
        // //1) the minimum of all P-values(as a baseline) (deprecated)
        // else if(aggregationMode == Parameters::AGGREGATION_MODE_MIN_PVAL){
        //     unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        //     double minLogPval = 0;
        //     for (size_t i = 0; i < dataToAggregate.size(); ++i) { 
        //         double currentLogPval = std::strtod(dataToAggregate[i][2].c_str(), NULL);
        //         if (currentLogPval < minLogPval) {
        //             minLogPval = currentLogPval;
        //         };
        //     }
        //     updatedPval = 1 -  exp( - exp(minLogPval) * orfCount);    
        // }

        //2) the P-value for the product-of-P-values
        else if (aggregationMode == Parameters::AGGREGATION_MODE_PRODUCT){
            if(dataToAggregate.size() == 0){
                return "";
            }
            double  sumLogPval= 0;
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][2].c_str(), NULL);
                sumLogPval += logPvalue;
                entries.push_back(&dataToAggregate[i]);
            }
            buffer.append(SSTR(dataToAggregate.size()));
            buffer.append("\t");
            buffer.append(SSTR(exp(sumLogPval) * numTargetSets));
        }

        //3) the P-values of the truncated product method 
        else if(aggregationMode == Parameters::AGGREGATION_MODE_TRUNCATED_PRODUCT){
            double logPvalThreshold = log(alpha / (orfCount + 1));
            double sumLogPval = 0; 
            size_t k = 0;
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][2].c_str(), NULL);
                if (logPvalue < logPvalThreshold) {
                    sumLogPval += logPvalue;
                    k++;
                }
            }
            if(k == 0){
                return "";
            }
            else {
                buffer.append(SSTR(k));
                buffer.append("\t");
                buffer.append(SSTR(exp(sumLogPval)));
            }
        }
        else {
            Debug(Debug::ERROR) << "Invalid aggregation function!\n";
            EXIT(EXIT_FAILURE);
        }
        
        next:
        buffer.append(",");

        for (size_t j = 0; j < entries.size(); j++) {
            // Aggregate the full line into string
            for (size_t i = 0; i < entries[j]->size(); ++i) {
                if (i == 2) {
                    double logPval = std::strtod(entries[j]->at(i).c_str(), NULL);
                    char tmpBuf[15];
                    sprintf(tmpBuf, "%.3E", exp(logPval));
                    buffer.append(tmpBuf);
                } else {
                    buffer.append(entries[j]->at(i));
                }
                if (i != (entries[j]->size() - 1)) {
                    buffer.append(1, '\t');
                }
            }
            if (j != entries.size() -1 ){
                buffer.append(1, '\n');
            }
        }
        return buffer;
    }

private:
    double alpha;
    int aggregationMode;
    bool filterSelfMatch;
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
    double* lGammaLookup;
    double** logBiLookup;
};

int combinehits(int argc, const char **argv, const Command &command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, par.alpha, (unsigned int) par.threads, par.compressed, par.aggregationMode, par.filterSelfMatch);
    return aggregation.runWithHeader();
}
