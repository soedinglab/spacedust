#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "iterativeclustersearch.sh.h"

void setIterativeClusterSearchWorkflowDefaults(Parameters *p) {
    p->sensitivity = 5.7;
    // TODO: Check query cov maybe?
    p->covThr = 0.8;
    p->evalProfile = 0.1;
    // TODO: Needs to be more than the count of target sets (10x?)
    //p->maxSequences = 1500;
    // TODO: Why??
    //p->scoreBias = 0.3;
    p->simpleBestHit = false;
    p->alpha = 1;
    // TODO: add a minimum alignment length cutoff, 4 residue alignments dont seem useful
    p->alnLenThr = 30;
    // Set alignment mode
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV; 
    p->addBacktrace = 1;
    p->dbOut = 0;
}

int iterativeclustersearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setIterativeClusterSearchWorkflowDefaults(&par);

    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.splitsequence.size(); i++) {
        par.splitsequence[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);

    if (FileUtil::directoryExists(par.db4.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << par.db4 << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(par.db4.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << par.db4 << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << par.db4 << "\n";
        }
    }
    size_t hash = par.hashParameter(command.databases, par.filenames, par.iterativeclusearchworkflow);
    std::string tmpDir = par.db4 + "/" + SSTR(hash);
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
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    // if (par.removeTmpFiles) {
    //     cmd.addVariable("REMOVE_TMP", "TRUE");
    // }
    cmd.addVariable("NUM_IT", SSTR(par.numIterations).c_str());
    cmd.addVariable("SUBSTRACT_PAR", par.createParameterString(par.subtractdbs).c_str());
    double originalEval = par.evalThr;
    par.evalThr = (par.evalThr < par.evalProfile) ? par.evalThr  : par.evalProfile;
    for (int i = 0; i < par.numIterations; i++) {
        if (i == 0) {
            par.realign = true;
        }
        if (i > 0) {
            //par.queryProfile = true;
            par.realign = false;
        }
        if (i == (par.numIterations - 1)) {
            par.evalThr = originalEval;
        }
        cmd.addVariable(std::string("PREFILTER_PAR_" + SSTR(i)).c_str(),
                        par.createParameterString(par.prefilter).c_str());
        cmd.addVariable(std::string("ALIGNMENT_PAR_" + SSTR(i)).c_str(),
                        par.createParameterString(par.align).c_str());
        cmd.addVariable(std::string("PROFILE_PAR_" + SSTR(i)).c_str(),
                        par.createParameterString(par.result2profile).c_str());
    }
    cmd.addVariable("MERGE_PAR", par.createParameterString(par.mergedbs).c_str());
    //cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    cmd.addVariable("BESTHITBYSET_PAR", par.createParameterString(par.besthitbyset).c_str());
    cmd.addVariable("COMBINEHITS_PAR", par.createParameterString(par.combinehits).c_str());
    cmd.addVariable("CLUSTERHITS_PAR", par.createParameterString(par.clusterhits).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("ALIGN_MODULE", "align");

    FileUtil::writeFile(tmpDir + "/iterativeclustersearch.sh", iterativeclustersearch_sh, iterativeclustersearch_sh_len);
    std::string program(tmpDir + "/iterativeclustersearch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
