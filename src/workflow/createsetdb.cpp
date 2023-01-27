#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "DBReader.h"
#include "CommandCaller.h"

#include "createsetdb.sh.h"

void setCreatesetdbWorkflowDefaults(Parameters *p) {
    p->orfMinLength = 30;
    //TODO: for protein seq we have to force --shuffle 0 
    p->shuffleDatabase = 0;
}

int createsetdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setCreatesetdbWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.filenames.back();
    par.filenames.pop_back();

    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        Debug(Debug::INFO) << "Tmp " << tmpDir << " folder does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created dir " << tmpDir << "\n";
        }
    }
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.createsetdb));
    if(par.reuseLatest == true){
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest" );
    }
    tmpDir = tmpDir + "/" + hash;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Can not create sub tmp folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    FileUtil::symlinkAlias(tmpDir, "latest");

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();

    std::string gffDir = par.gffDir;

    CommandCaller cmd;
    cmd.addVariable("OUTDB", outDb.c_str());
    cmd.addVariable("GFFDIR", gffDir.c_str());
    cmd.addVariable("TMP_PATH", tmpDir.c_str());

    if (par.removeTmpFiles) {
        cmd.addVariable("REMOVE_TMP", "TRUE");
    }

    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("GFF2DB_PAR", par.createParameterString(par.gff2db).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    par.stat = "linecount";
    cmd.addVariable("RESULT2STATS_PAR", par.createParameterString(par.result2stats).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    FileUtil::writeFile(tmpDir + "/createsetdb.sh", createsetdb_sh, createsetdb_sh_len);
    std::string program(tmpDir + "/createsetdb.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
