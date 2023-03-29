#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

std::map<unsigned int, std::string> readSetToSource(const std::string& file);

int summarizeresults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> qlookupReader(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    qlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* qlookup = qlookupReader.getLookup();

    DBReader<unsigned int> tlookupReader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_LOOKUP);
    tlookupReader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* tlookup = tlookupReader.getLookup();

    std::map<unsigned int, std::string> qSetToSource;
    std::map<unsigned int, std::string> tSetToSource;
    std::string file1 = par.db1 + ".source";
    std::string file2 = par.db2 + ".source";
    qSetToSource = readSetToSource(file1);
    tSetToSource = readSetToSource(file2);

    DBReader<unsigned int> hdrReader(par.hdr3.c_str(), par.hdr3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    hdrReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> alnReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads, shouldCompress, dbType);
    dbw.open();

    Debug::Progress progress(hdrReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char *entry[255];

        std::string buffer;
        buffer.reserve(1024 * 1024);

        std::string tmpBuffer;
        tmpBuffer.reserve(1024 * 1024);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < hdrReader.getSize(); ++id) {
            progress.updateProgress();

            unsigned int matchKey = hdrReader.getDbKey(id);
            char *alnData = alnReader.getData(id, thread_idx);
            char *matchData = hdrReader.getData(id, thread_idx);
            while (*matchData != '\0') {
                size_t columns = Util::getWordsOfLine(matchData, entry, 255);
                matchData = Util::skipLine(matchData);
                if (columns < 5) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                char* alnCurrent = alnData;
                size_t qsetid =  Util::fast_atoi<unsigned int>(entry[0]);
                size_t tsetid =  Util::fast_atoi<unsigned int>(entry[1]);
                buffer.append("#");
                buffer.append(SSTR(id));
                buffer.append("\t");
                buffer.append(SSTR(qSetToSource[qsetid]));
                buffer.append("\t");
                buffer.append(SSTR(tSetToSource[tsetid]));
                buffer.append("\t");
                buffer.append(entry[2], entry[3] - entry[2] - 1);
                buffer.append("\t");
                buffer.append(entry[3], entry[4] - entry[3] - 1);
                buffer.append("\t");
                buffer.append(entry[4], entry[5] - entry[4]);
                buffer.append("\n");

                while (*alnCurrent != '\0') {
                    // // read only key
                    // unsigned int dbKey = Util::fast_atoi<unsigned int>(alnCurrent);

                    //entry is overwritten by getWordsOfLine
                    columns = Util::getWordsOfLine(alnCurrent, entry, 255);
                    if (columns < 10) {
                        Debug(Debug::ERROR) << "Invalid alignment result record\n";
                        EXIT(EXIT_FAILURE);
                    }
                    alnCurrent = Util::skipLine(alnCurrent);
                    size_t qid = Util::fast_atoi<unsigned int>(entry[0]);
                    size_t tid = Util::fast_atoi<unsigned int>(entry[1]);

                    tmpBuffer.append(">");
                    tmpBuffer.append(SSTR(qlookup[qid].entryName)); //queryname
                    tmpBuffer.append("\t");
                    tmpBuffer.append(SSTR(tlookup[tid].entryName)); //targetname
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[2], alnCurrent - entry[2]); //rest of the line
                }
                buffer.append(tmpBuffer);
                tmpBuffer.clear();
            }
            dbw.writeData(buffer.c_str(), buffer.length(), matchKey, thread_idx, par.dbOut);
            buffer.clear();
        }
    }
    dbw.close(par.dbOut == false);
    if (par.dbOut == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    hdrReader.close();
    alnReader.close();
    qlookupReader.close();
    tlookupReader.close();

    return EXIT_SUCCESS;
}

