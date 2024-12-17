#include "LocalParameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "DBReader.h"
#include "CommandCaller.h"
#include "PatternCompiler.h"

#include "createsetdb.sh.h"

#include <dirent.h>
#include <sys/stat.h>

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

    PatternCompiler include(par.fileInclude.c_str());
    PatternCompiler exclude(par.fileExclude.c_str());
    
    for (size_t i = 1; i < par.filenames.size(); ++i) {
        if (FileUtil::directoryExists(par.filenames[i].c_str()) || Util::endsWith(".tsv", par.filenames[i].c_str())) {
            Debug(Debug::ERROR) << "Only one directory or tsv file (" << par.filenames[i] << ") or a list of files can be given\n";
            EXIT(EXIT_FAILURE);
        }
    }
    
    if (FileUtil::directoryExists(par.filenames[0].c_str())) {
        if (par.filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one directory can be given\n";
            EXIT(EXIT_FAILURE);
        }
        std::vector<std::string> dirs;
        dirs.push_back(par.filenames.back());
        par.filenames.pop_back();
        while (dirs.size() != 0) {
            std::string dir = dirs.back();
            dirs.pop_back();
            DIR* handle = opendir(dir.c_str());
            if (handle == NULL) {
                continue;
            }
            while (dirent* entry = readdir(handle)) {
                std::string filename(entry->d_name);
                if (filename != "." && filename != "..") {
                    std::string fullpath = dir + "/" + filename;
                    struct stat info;
                    stat(fullpath.c_str(), &info);
                    if (info.st_mode & S_IFDIR) {
                        dirs.push_back(fullpath);
                    } else if (include.isMatch(filename.c_str()) == true && exclude.isMatch(filename.c_str()) == false) {
                        par.filenames.push_back(fullpath);
                    }
                }
            }
            closedir(handle);
        }
    } else if (Util::endsWith(".tsv", par.filenames[0])) {
        if (par.filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one tsv file can be given\n";
            EXIT(EXIT_FAILURE);
        }
        std::string tsv = par.filenames.back();
        par.filenames.pop_back();
        FILE* file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
        char* line = NULL;
        size_t len = 0;
        ssize_t read;
        while ((read = getline(&line, &len, file)) != -1) {
            if (line[read - 1] == '\n') {
                line[read - 1] = '\0';
                read--;
            }
            par.filenames.push_back(line);
        }
        free(line);
        fclose(file);
    }

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
