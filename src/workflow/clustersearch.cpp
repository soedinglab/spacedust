#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "CommandCaller.h"

#include "clustersearch.sh.h"

void setClusterSearchWorkflowDefaults(LocalParameters *p) {
    p->sensitivity = 5.7;
    // TODO: Check query cov maybe?
    p->covMode = Parameters::COV_MODE_BIDIRECTIONAL;
    p->covThr = 0.8;
    p->evalThr = 10;
    
    // TODO: Needs to be more than the count of target sets (10x?)
    //p->maxSequences = 1500;

    // TODO: Why??
    //p->scoreBias = 0.3;
    p->simpleBestHit = true;
    p->alpha = 1;
    // TODO: add a minimum alignment length cutoff, 4 residue alignments dont seem useful
    p->alnLenThr = 30;
    // Set alignment mode
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV; 
    p->addBacktrace = 1;
    p->dbOut = 1;
    //profile-sequence clustersearch
    if(p->profileClusterSearch == 1){
        p->exhaustiveSearch = 1;
    }
}

int clustersearch(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setClusterSearchWorkflowDefaults(&par);

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
    size_t hash = par.hashParameter(command.databases, par.filenames, par.clustersearchworkflow);
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
    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }
    cmd.addVariable("USE_PROFILE", par.profileClusterSearch == 1 ? "TRUE" : NULL);
    cmd.addVariable("USE_FOLDSEEK", par.clusterSearchMode == 1 ? "TRUE" : NULL);
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clusterworkflow).c_str());
    std::vector<MMseqsParameter*> searchwithoutnumiter;
    for (size_t i = 0; i < par.clustersearchworkflow.size(); i++) {
        if (par.clustersearchworkflow[i]->uniqid != par.PARAM_NUM_ITERATIONS.uniqid) {
            searchwithoutnumiter.push_back(par.clustersearchworkflow[i]);
        }
    }
    cmd.addVariable("SEARCH_PAR", par.createParameterString(searchwithoutnumiter).c_str());
    cmd.addVariable("FOLDSEEKSEARCH_PAR", par.createParameterString(par.foldseeksearch).c_str());
    cmd.addVariable("BESTHITBYSET_PAR", par.createParameterString(par.besthitbyset).c_str());
    cmd.addVariable("COMBINEHITS_PAR", par.createParameterString(par.combinehits).c_str());
    cmd.addVariable("CLUSTERHITS_PAR", par.createParameterString(par.clusterhits).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/clustersearch.sh", clustersearch_sh, clustersearch_sh_len);
    std::string program(tmpDir + "/clustersearch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
