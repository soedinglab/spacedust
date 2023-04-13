#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "clusterdb.sh.h"

void setclusterDbDefaults(LocalParameters *p) {
    if(p->clusterSearchMode == 1){
        p->seqIdThr = 0.5;
    }else{
        p->seqIdThr = 0.7;
    }
    p->covThr=0.8;
    p->covMode = Parameters::COV_MODE_BIDIRECTIONAL;
}

int clusterdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setclusterDbDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (FileUtil::directoryExists(par.db2.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db2 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db2.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db2 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db2 << "\n";
        }
    }
    size_t hash = par.hashParameter(command.databases, par.filenames, par.clusterdb);
    std::string tmpDir = par.db2 + "/" + SSTR(hash);
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
    cmd.addVariable("USE_FOLDSEEK", par.clusterSearchMode == 1 ? "TRUE" : NULL);
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clusterworkflow).c_str());
    cmd.addVariable("CONSENSUS_PAR", par.createParameterString(par.profile2seq).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());

    FileUtil::writeFile(tmpDir + "/clusterdb.sh", clusterdb_sh, clusterdb_sh_len);
    std::string program(tmpDir + "/clusterdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
