#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "clusterdb.sh.h"

void setclusterDbDefaults(LocalParameters *p) {
    p->seqIdThr = 0.7;
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

    bool useFoldseek = false;
    if (par.clusterSearchMode == 1) {
        useFoldseek = true;
        struct stat st;
        if (stat(par.foldseekPath.c_str(), &st) != 0) {
            Debug(Debug::ERROR) << "Cannot find foldseek binary " << par.foldseekPath << ".\n";
            EXIT(EXIT_FAILURE);
        }
        bool isExecutable = (st.st_mode & S_IXUSR) || (st.st_mode & S_IXGRP) || (st.st_mode & S_IXOTH);
        if (isExecutable == false) {
            Debug(Debug::ERROR) << "Cannot execute foldseek binary " << par.foldseekPath << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    cmd.addVariable("FOLDSEEK", par.foldseekPath.c_str());
    cmd.addVariable("USE_FOLDSEEK", useFoldseek ? "TRUE" : NULL);
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clusterworkflow).c_str());
    cmd.addVariable("CONSENSUS_PAR", par.createParameterString(par.profile2seq).c_str());
    cmd.addVariable("PROFILE_PAR", par.createParameterString(par.result2profile).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("MERGECLU_PAR", par.createParameterString(par.threadsandcompression).c_str());
    par.pca = 1.4;
    par.pcb = 1.5;
    par.scoringMatrixFile = "3di.out";
    par.seedScoringMatrixFile = "3di.out";
    par.maskProfile = 0;
    par.compBiasCorrection = 0;
    if(par.PARAM_E_PROFILE.wasSet == false){
        par.evalProfile = 0.1;
        par.evalThr = 0.1;
    }
    std::vector<MMseqsParameter*> result2profile_ss;
    for (size_t i = 0; i < par.result2profile.size(); i++) {
        if (par.result2profile[i]->uniqid != par.PARAM_GAP_PSEUDOCOUNT.uniqid) {
            result2profile_ss.push_back(par.result2profile[i]);
        }
    }
    cmd.addVariable("PROFILE_SS_PAR", par.createParameterString(result2profile_ss).c_str());

    FileUtil::writeFile(tmpDir + "/clusterdb.sh", clusterdb_sh, clusterdb_sh_len);
    std::string program(tmpDir + "/clusterdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
