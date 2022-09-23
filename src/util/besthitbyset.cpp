#include "Debug.h"
#include "LocalParameters.h"
#include "Aggregation.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

double ComputelogPval(double eval, double logCalibration){
    if(eval == 0) {
        return log(DBL_MIN)-logCalibration;
    }
    else if (eval > 0 && eval < 10e-4){
        return log(eval)-logCalibration;
    }
    else {
        return log(1 - exp(-eval))-logCalibration;
    }
}

class BestHitBySetFilter : public Aggregation {
public :
    BestHitBySetFilter(const std::string &targetDbName, const std::string &resultDbName,
                       const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads, unsigned int compressed) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), simpleBestHitMode(simpleBestHitMode) {
        std::string sizeDbName = targetDbName + "_set_size";
        std::string sizeDbIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~BestHitBySetFilter() {
        targetSizeReader->close();
        delete targetSizeReader;
    }


    void prepareInput(unsigned int, unsigned int) {}

    std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int, unsigned int targetSetKey, unsigned int thread_idx)  {
        double bestScore = -DBL_MAX;
        double secondBestScore = -DBL_MAX;
        double bestEval = DBL_MAX;
        double logBestHitCalibration = log(1);

        // Look for the lowest p-value and retain only this line
        // dataToAggregate = [nbrTargetGene][Field of result]
        size_t targetId = targetSizeReader->getId(targetSetKey);
        if (targetId == UINT_MAX) {
            Debug(Debug::ERROR) << "Invalid target size database key " << targetSetKey << ".\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = targetSizeReader->getData(targetId, thread_idx);
        unsigned int nbrGenes = Util::fast_atoi<unsigned int>(data);

        std::vector<std::string> *bestEntry = NULL;
        std::vector<std::vector<std::string>*> allBestEntries;
        std::vector<double> logCorrectedPvalList;
        std::vector<double> bestEvalList;
        for (size_t i = 0; i < dataToAggregate.size(); i++) {
            double eval = strtod(dataToAggregate[i][4].c_str(), NULL);
            double pval = eval/nbrGenes;
            //prevent log(0)
            if (pval == 0) {
                pval = DBL_MIN;
            }
            double score = std::min(DBL_MAX, -log(eval));

            //if only one hit use simple best hit
            if(simpleBestHitMode ||dataToAggregate.size() < 2) {
                if( eval < bestEval){
                    bestEval = eval;
                    bestEntry = &dataToAggregate[i];
                }
            }
            else {
                if (score >= bestScore) {
                    secondBestScore = bestScore;
                    bestScore = score;
                    bestEntry = &dataToAggregate[i];
                } 
                else if (score > secondBestScore) {
                    secondBestScore = score;
                }
            }
        }
        //TODO: includeParalog with a factor e.g. 100. Go through the list again if this option is turned on for those with more than 1 hits
        if(false && simpleBestHitMode && dataToAggregate.size() > 1){
            double evalThr = bestEval * 100; //TODO: hardcoded
            for (size_t i = 0; i < dataToAggregate.size(); i++) {
                double eval = strtod(dataToAggregate[i][4].c_str(), NULL);
                if( eval <= evalThr){
                    allBestEntries.push_back(&dataToAggregate[i]);
                    bestEvalList.push_back(eval);
                }
            }
        } else{
            allBestEntries.push_back(bestEntry);
        }


        if(allBestEntries.size() > 1){
            for (size_t i = 0; i < allBestEntries.size(); i++) {
                logCorrectedPvalList.push_back(ComputelogPval(bestEvalList[i],logBestHitCalibration));
            }
        } else {
            if (simpleBestHitMode ||dataToAggregate.size() < 2) {
                logCorrectedPvalList.push_back(ComputelogPval(bestEval,logBestHitCalibration));
            } 
            else {
                logCorrectedPvalList.push_back(secondBestScore - bestScore);
            }
        }


        if (bestEntry == NULL) {
            return "";
        }

        std::string buffer;
        buffer.reserve(1024);

        for (size_t j = 0; j < allBestEntries.size(); j++) {
            // Aggregate the full line into string
            for (size_t i = 0; i < allBestEntries[j]->size(); ++i) {
                if (i == 2) {
                    char tmpBuf[15];
                    sprintf(tmpBuf, "%.3E", logCorrectedPvalList[j]);
                    buffer.append(tmpBuf);
                } else {
                    buffer.append(allBestEntries[j]->at(i));
                }
                if (i != (allBestEntries[j]->size() - 1)) {
                    buffer.append(1, '\t');
                }
            }
            if (j != allBestEntries.size() -1 ){
                buffer.append(1, '\n');
            }
        }

        return buffer;
    }

private:
    DBReader<unsigned int> *targetSizeReader;
    bool simpleBestHitMode;
};


int besthitbyset(int argc, const char **argv, const Command &command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    BestHitBySetFilter aggregation(par.db2, par.db3, par.db4, par.simpleBestHit, (unsigned int) par.threads, par.compressed);
    return aggregation.run();
}
