#include "Aggregation.h"
#include "Util.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif

Aggregation::Aggregation(const std::string &targetDbName, const std::string &resultDbName,
                         const std::string &outputDbName, unsigned int threads, unsigned int compressed)
        : resultDbName(resultDbName), outputDbName(outputDbName), threads(threads), compressed(compressed) {
    std::string sizeDbName = targetDbName + "_member_to_set";
    std::string sizeDbIndex = targetDbName + "_member_to_set.index";
    targetSetReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    targetSetReader->open(DBReader<unsigned int>::NOSORT);
}

Aggregation::~Aggregation() {
    targetSetReader->close();
    delete targetSetReader;
}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
void Aggregation::buildMap(char *data, int thread_idx, std::map<unsigned int, std::vector<std::vector<std::string>>> &dataToAggregate) {
    while (*data != '\0') {
        char *current = data;
        data = Util::skipLine(data);
        size_t length = data - current;
        std::string line(current, length - 1);
        if (line.empty() == true) {
            continue;
        }

        std::vector<std::string> columns = Util::split(line, "\t");
        unsigned int targetKey = Util::fast_atoi<unsigned int>(columns[1].c_str());
        size_t setId = targetSetReader->getId(targetKey);
        if (setId == UINT_MAX) {
            Debug(Debug::ERROR) << "Invalid target database key " << columns[1] << ".\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = targetSetReader->getData(setId, thread_idx);
        unsigned int setKey = Util::fast_atoi<unsigned int>(data);
        dataToAggregate[setKey].push_back(columns);
    }
}

int Aggregation::run() {
    std::string inputDBIndex = resultDbName + ".index";
    DBReader<unsigned int> reader(resultDbName.c_str(), inputDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string outputDBIndex = outputDbName + ".index";
    DBWriter writer(outputDbName.c_str(), outputDBIndex.c_str(), threads, compressed, Parameters::DBTYPE_ALIGNMENT_RES); //TODO: check if this is true
    writer.open();
    Debug::Progress progress(reader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);

        std::map<unsigned int, std::vector<std::vector<std::string>>> dataToMerge;
#pragma omp for
        for (size_t i = 0; i < reader.getSize(); i++) {
            progress.updateProgress();
            dataToMerge.clear();

            unsigned int key = reader.getDbKey(i);
            buildMap(reader.getData(i, thread_idx), thread_idx, dataToMerge);
            prepareInput(key, thread_idx);
            
            for (std::map<unsigned int, std::vector<std::vector<std::string>>>::const_iterator it = dataToMerge.begin();
                 it != dataToMerge.end(); ++it) {
                unsigned int targetKey = it->first;
                std::vector<std::vector<std::string>> columns = it->second;
                buffer.append(aggregateEntry(columns, key, targetKey, thread_idx));
                buffer.append("\n");
            }
            writer.writeData(buffer.c_str(), buffer.length(), key, thread_idx);
            buffer.clear();
        }
    };
    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

int Aggregation::runWithHeader() {
    std::string inputDBIndex = resultDbName + ".index";
    DBReader<unsigned int> reader(resultDbName.c_str(), inputDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string outputDBIndex = outputDbName + ".index";
    DBWriter writer(outputDbName.c_str(), outputDBIndex.c_str(), threads, compressed, Parameters::DBTYPE_ALIGNMENT_RES); //TODO: check if this is true
    writer.open();

    std::string outputHdrName = outputDbName + "_h";
    std::string outputHdrIndex = outputDbName + "_h.index";
    DBWriter hdrwriter(outputHdrName.c_str(), outputHdrIndex.c_str(), threads, compressed, Parameters::DBTYPE_GENERIC_DB);
    hdrwriter.open();
    Debug::Progress progress(reader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);
        std::string header;
        header.reserve(1024);
        unsigned int match_idx = 0;

        std::map<unsigned int, std::vector<std::vector<std::string>>> dataToMerge;
#pragma omp for
        for (size_t i = 0; i < reader.getSize(); i++) {
            progress.updateProgress();
            dataToMerge.clear();

            unsigned int querykey = reader.getDbKey(i);
            buildMap(reader.getData(i, thread_idx), thread_idx, dataToMerge);
            prepareInput(querykey, thread_idx);

            for (std::map<unsigned int, std::vector<std::vector<std::string>>>::const_iterator it = dataToMerge.begin();
                 it != dataToMerge.end(); ++it) {
                unsigned int targetKey = it->first;
                std::vector<std::vector<std::string>> columns = it->second;
                std::string entryLine = aggregateEntry(columns, querykey, targetKey, thread_idx);
                if (entryLine ==""){
                    continue;
                }
                //header info is printed in the same string separated by comma
                std::vector<std::string> column = Util::split(entryLine, ",");
                header.append(column[0]);
                header.append("\n");
                buffer.append(column[1]);
                buffer.append("\n");
                unsigned int key = __sync_fetch_and_add(&(match_idx), 1);
                hdrwriter.writeData(header.c_str(), header.length(), key, thread_idx);
                writer.writeData(buffer.c_str(), buffer.length(), key, thread_idx);
                buffer.clear();
                header.clear();
            }
        }
    };
    hdrwriter.close(true);
    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}
