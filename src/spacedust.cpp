#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "Prefiltering.h"
#include "DownloadDatabase.h"

const char* binary_name = "spacedust";
const char* tool_name = "spacedust";
const char* tool_introduction = "Spacedust is a tool to discover conserved gene clusters between any pairs of contig/genomes";
const char* main_author = "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de>";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
extern const char* MMSEQS_CURRENT_INDEX_VERSION;
const char* index_version_compatible = MMSEQS_CURRENT_INDEX_VERSION;
bool hide_base_commands = true;
void (*validatorUpdate)(void) = 0;
std::vector<DatabaseDownload> externalDownloads = {};
std::vector<KmerThreshold> externalThreshold = {};
bool hide_base_downloads = false;

std::vector<int> FlatfileAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_DIRECTORY};

LocalParameters& localPar = LocalParameters::getLocalInstance();

extern std::vector<Command> baseCommands;
std::vector<Command> spacedustCommands = {
        // {"easy-search",             easysearch,           &localPar.easyclustersearch, COMMAND_EASY,
        //         "Find clusters of colocalized hits between any query-target set pairs from FASTA input",
        //         NULL,
        //         "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de>",
        //         "<i:fastaFile1[.txt]> ... <i:fastaFileN[.txt]> <i:targetDB> <o:output[.tsv]> <tmpDir>",
        //         CITATION_MMSEQS2, {{"fast[a|q]File[.gz|bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
        //                                {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
        //                                {"result", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
        //                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"createsetdb",           createsetdb,           &localPar.createsetdb,           COMMAND_MAIN,
                "Create sequence set database from FASTA (and GFF3) input of contigs/genomes",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:fastaFile1[.gz|bz2]> ... <i:fastaFileN[.gz|bz2]>|<directory>|<listOfFastaFiles.tsv> <o:setDB> <tmpDir>",
                CITATION_MMSEQS2, {{"fast[a|q]File[.gz|bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &FlatfileAndFolder },
                                                           {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"gff2db",                gff2db,                &localPar.gff2db,                COMMAND_SPECIAL,
                "Extract regions from a sequence database based on a GFF3 file",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de>",
                "<i:gff3File1> ... <i:gff3FileN> <i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"gff3File", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                                                           {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},
        {"aa2foldseek",           aa2foldseek,           &localPar.aa2foldseek,           COMMAND_MAIN,
                "Map a sequence DB to reference foldseek DB",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:inputDB> <i:targetDB> <tmpDir>",
                CITATION_MMSEQS2, {{"inputDB", DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"clusterdb",           clusterdb,           &localPar.clusterdb,           COMMAND_MAIN,
                "Build a searchable cluster database from sequence DB or foldseek structure DB",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:inputDB> <tmpDir>",
                CITATION_MMSEQS2, {{"inputDB", DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"clustersearch",       clustersearch,       &localPar.clustersearchworkflow,       COMMAND_MAIN,
                "Find clusters of colocalized hits between any query-target sequence/profile set database",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <o:output[.tsv]> <tmpDir>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"output", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"besthitbyset",        besthitbyset,        &localPar.besthitbyset,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "For each set of sequences compute the best element and update p-value",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                " <i:targetSetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"combinehits",    combinehits,    &localPar.combinehits,     COMMAND_SPECIAL | COMMAND_EXPERT,
                "Group hits and compute a combined E-value for each query-target set pair",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <i:resultDB> <o:pvalDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"pvalDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"summarizeresults",    summarizeresults,    &localPar.clusterhits,     COMMAND_SPECIAL | COMMAND_EXPERT,
                "Summarize results on clustered hits",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de>",
                "<i:queryDB> <i:targetDB> <i:alnDB> <o:output[.tsv]>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"alnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}, 
                                       {"output", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"clusterhits",    clusterhits,    &localPar.clusterhits,      COMMAND_SPECIAL | COMMAND_EXPERT,
                "Find clusters of hits by agglomerative hierarchical clustering and compute their clustering and ordering P-values",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <i:resultDB> <o:outDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_LOOKUP, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_LOOKUP, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::resultDb },
                                                            {"outDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
};

void init() {
    registerCommands(&baseCommands);
    registerCommands(&spacedustCommands);
}
void (*initCommands)(void) = init;

void initParameterSingleton() { 
    LocalParameters::initInstance(); 
}
