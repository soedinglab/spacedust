#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "aa2foldseek.sh.h"

void setaa2FoldseekDefaults(Parameters *p) {
    p->exactKmerMatching = 1;
    p->maxResListLen = 10;
    //seqid = 0.5 && cov = 0.9 for more sensitivity?
    p->seqIdThr = 0.9;
    p->covThr=0.9;
}

int aa2foldseek(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setaa2FoldseekDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (FileUtil::directoryExists(par.db3.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db3 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db3.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db3 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db3 << "\n";
        }
    }
    size_t hash = par.hashParameter(command.databases, par.filenames, par.aa2foldseek);
    std::string tmpDir = par.db3 + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    FileUtil::symlinkAlias(tmpDir, "latest");

    CommandCaller cmd;
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    par.stat = "linecount";
    cmd.addVariable("RESULT2STATS_PAR", par.createParameterString(par.result2stats).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/aa2foldseek.sh", aa2foldseek_sh, aa2foldseek_sh_len);
    std::string program(tmpDir + "/aa2foldseek.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
