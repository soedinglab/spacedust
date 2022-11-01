#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "Prefiltering.h"
#include "DownloadDatabase.h"

const char* binary_name = "clustersearch";
const char* tool_name = "clustersearch";
const char* tool_introduction = "clustersearch is a tool to discover conserved gene clusters between any pairs of contig/genomes";
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

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
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
                "<i:fastaFile1[.gz|bz2]> ... <i:fastaFileN[.gz|bz2]> <o:setDB> <tmpDir>",
                CITATION_MMSEQS2, {{"fast[a|q]File[.gz|bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                                           {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"aa2foldseek",           aa2foldseek,           &localPar.aa2foldseek,           COMMAND_MAIN,
                "Map an sequence DB to reference foldseek DB",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:inputDB> <i:targetDB> <o:outDB> <tmpDir>",
                CITATION_MMSEQS2, {{"inputDB", DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"outDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"clustersearch",       clustersearch,       &localPar.clustersearchworkflow,       COMMAND_MAIN,
                "Find clusters of colocalized hits between any query-target sequence/profile set database",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <o:resultDB> <tmpDir>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"iterativeclustersearch",       iterativeclustersearch,       &localPar.iterativeclusearchworkflow,       COMMAND_SPECIAL | COMMAND_EXPERT,
                "Iterative profile cluster search",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <o:resultDB> <tmpDir>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
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
        {"counthits",    counthits,    &localPar.counthits,     COMMAND_SPECIAL | COMMAND_EXPERT,
                "Count hits by query and filter",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpinat.mpg.de> & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:resultDB> <o:outDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"outDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
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

