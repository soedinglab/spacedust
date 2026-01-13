#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H
 
#include "Parameters.h"
#include "FileUtil.h"

const int CITATION_SPACEDUST = CITATION_END;

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> createsetdb; 
    std::vector<MMseqsParameter*> aa2foldseek;    
    std::vector<MMseqsParameter*> clustersearchworkflow;
    std::vector<MMseqsParameter*> besthitbyset;
    std::vector<MMseqsParameter*> combinehits;
    std::vector<MMseqsParameter*> clusterhits;
    std::vector<MMseqsParameter*> foldseeksearch;
    std::vector<MMseqsParameter*> counthits;
    std::vector<MMseqsParameter*> clusterdb;

    PARAMETER(PARAM_CLUSTERSEARCH_MODE)
    PARAMETER(PARAM_SUBOPTIMAL_HITS)
    PARAMETER(PARAM_FILTER_SELF_MATCH)
    PARAMETER(PARAM_MULTIHIT_PVAL)
    PARAMETER(PARAM_CLUSTER_PVAL)
    PARAMETER(PARAM_MAX_GENE_GAP)
    PARAMETER(PARAM_CLUSTER_SIZE)
    PARAMETER(PARAM_CLUSTER_USE_WEIGHT)
    PARAMETER(PARAM_PROFILE_CLUSTER_SEARCH)
    PARAMETER(PARAM_FILE_INCLUDE)
    PARAMETER(PARAM_FILE_EXCLUDE)
    PARAMETER(PARAM_GFF_DIR)
    PARAMETER(PARAM_FOLDSEEK_PATH)

    int clusterSearchMode;
    float pMHThr;
    float pCluThr;
    bool clusterUseWeight;
    bool profileClusterSearch;
    int maxGeneGaps;
    int clusterSize;
    std::string gffDir;
    bool filterSelfMatch;
    int suboptHitsFactor;
    std::string fileInclude;
    std::string fileExclude;
    std::string foldseekPath;

    static std::string getAbsExePath();

    LocalParameters() : 
        Parameters(),
        PARAM_CLUSTERSEARCH_MODE(PARAM_CLUSTERSEARCH_MODE_ID, "--search-mode", "Cluster Search Mode", "0: sequence search with MMseqs2, 1: structure comparison with Foldseek, 2: Foldseek + ProstT5", typeid(int), (void *) &clusterSearchMode, "^[0-2]{1}"),
        PARAM_SUBOPTIMAL_HITS(PARAM_SUBOPTIMAL_HITS_ID, "--suboptimal-hits", "Include sub-optimal hits with factor", "Include sub-optimal hits of query sequence up to a factor of its E-value. 0: only include one best hit", typeid(int), (void *) &suboptHitsFactor, "^(0|[1-9]{1}[0-9]*)$"),
        PARAM_FILTER_SELF_MATCH(PARAM_FILTER_SELF_MATCH_ID, "--filter-self-match", "Filter self match", "Remove hits between the same set. 0: do not filter, 1: filter", typeid(bool), (void *) &filterSelfMatch, ""),
        PARAM_MULTIHIT_PVAL(PARAM_MULTIHIT_PVAL_ID, "--multihit-pval", "Multihit P-value cutoff", "Multihit P-value threshold for cluster match output", typeid(float), (void *) &pMHThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_CLUSTER_PVAL(PARAM_CLUSTER_PVAL_ID, "--cluster-pval", "Clustering and Ordering P-value cutoff","Clustering and Ordering P-value threshold for cluster match output", typeid(float), (void *) &pCluThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_MAX_GENE_GAP(PARAM_MAX_GENE_GAP_ID, "--max-gene-gap", "Maximum gene gaps", "Maximum number of genes allowed between two clusters to merge", typeid(int), (void *) &maxGeneGaps, "^[1-9]{1}[0-9]*$"),
        PARAM_CLUSTER_SIZE(PARAM_CLUSTER_SIZE_ID, "--cluster-size", "Minimal cluster size", "Minimum number of genes to define cluster", typeid(int), (void *) &clusterSize, "^[1-9]{1}[0-9]*$"),
        PARAM_CLUSTER_USE_WEIGHT(PARAM_CLUSTER_USE_WEIGHT_ID, "--cluster-use-weight", "Cluster weighting factor ", "Use weighting factor to calculate cluster match score", typeid(bool), (void *) &clusterUseWeight, ""),
        PARAM_PROFILE_CLUSTER_SEARCH(PARAM_PROFILE_CLUSTER_SEARCH_ID, "--profile-cluster-search", "Cluster search against profiles", "Perform profile(target)-sequence searches in clustersearch", typeid(bool), (void *) &profileClusterSearch, ""),
        PARAM_FILE_INCLUDE(PARAM_FILE_INCLUDE_ID, "--file-include", "File Inclusion Regex", "Include file names based on this regex", typeid(std::string), (void *) &fileInclude, "^.*$"),
        PARAM_FILE_EXCLUDE(PARAM_FILE_EXCLUDE_ID, "--file-exclude", "File Exclusion Regex", "Exclude file names based on this regex", typeid(std::string), (void *) &fileExclude, "^.*$"),
        PARAM_GFF_DIR(PARAM_GFF_DIR_ID, "--gff-dir", "gff dir file", "Path to gff dir file", typeid(std::string), (void *) &gffDir, ""),
        PARAM_FOLDSEEK_PATH(PARAM_FOLDSEEK_PATH_ID, "--foldseek-path", "Path to Foldseek", "Path to Foldseek binary", typeid(std::string), (void *) &foldseekPath, "")
    {

        // clusterhits
        clusterhits.push_back(&PARAM_MULTIHIT_PVAL);
        clusterhits.push_back(&PARAM_CLUSTER_PVAL);
        clusterhits.push_back(&PARAM_MAX_GENE_GAP);
        clusterhits.push_back(&PARAM_CLUSTER_SIZE);
        clusterhits.push_back(&PARAM_CLUSTER_USE_WEIGHT);
        clusterhits.push_back(&PARAM_DB_OUTPUT);
        clusterhits.push_back(&PARAM_ALPHA);
        clusterhits.push_back(&PARAM_THREADS);
        clusterhits.push_back(&PARAM_COMPRESSED);
        clusterhits.push_back(&PARAM_V);


        // besthitperset
        besthitbyset.push_back(&PARAM_SIMPLE_BEST_HIT);
        besthitbyset.push_back(&PARAM_SUBOPTIMAL_HITS);
        besthitbyset.push_back(&PARAM_THREADS);
        besthitbyset.push_back(&PARAM_COMPRESSED);
        besthitbyset.push_back(&PARAM_V);


        // combinehits
        combinehits.push_back(&PARAM_ALPHA);
        combinehits.push_back(&PARAM_AGGREGATION_MODE);
        combinehits.push_back(&PARAM_FILTER_SELF_MATCH);
        combinehits.push_back(&PARAM_THREADS);
        combinehits.push_back(&PARAM_COMPRESSED);
        combinehits.push_back(&PARAM_V);

        // counthits
        counthits.push_back(&PARAM_DB_OUTPUT);
        counthits.push_back(&PARAM_THREADS);
        counthits.push_back(&PARAM_COMPRESSED);
        counthits.push_back(&PARAM_V);

        // foldseeksearch
        foldseeksearch.push_back(&PARAM_E);
        foldseeksearch.push_back(&PARAM_ADD_BACKTRACE);
        foldseeksearch.push_back(&PARAM_C);
        foldseeksearch.push_back(&PARAM_COV_MODE);
        foldseeksearch.push_back(&PARAM_NUM_ITERATIONS);
        foldseeksearch.push_back(&PARAM_MAX_SEQS);
        foldseeksearch.push_back(&PARAM_THREADS);
        foldseeksearch.push_back(&PARAM_COMPRESSED);
        foldseeksearch.push_back(&PARAM_V);

        // multi hit db
        createsetdb = combineList(createdb, extractorfs);
        createsetdb = combineList(createsetdb, translatenucs);
        createsetdb = combineList(createsetdb, gff2db);
        createsetdb = combineList(createsetdb, result2stats);
        createsetdb.push_back(&PARAM_FILE_INCLUDE);
        createsetdb.push_back(&PARAM_FILE_EXCLUDE);
        createsetdb.push_back(&PARAM_GFF_DIR);

        // multi hit search
        clustersearchworkflow = combineList(searchworkflow, besthitbyset);
        clustersearchworkflow = combineList(clustersearchworkflow, combinehits);
        clustersearchworkflow = combineList(clustersearchworkflow, clusterhits);
        clustersearchworkflow.push_back(&PARAM_PROFILE_CLUSTER_SEARCH);
        clustersearchworkflow.push_back(&PARAM_CLUSTERSEARCH_MODE);
        clustersearchworkflow.push_back(&PARAM_FOLDSEEK_PATH);

        //aa2foldseek
        aa2foldseek = combineList(prefilter, align);
        aa2foldseek = combineList(aa2foldseek,result2stats);
        aa2foldseek.push_back(&PARAM_REMOVE_TMP_FILES);

        //clusterdb
        clusterdb = combineList(clusterworkflow, profile2seq);
        clusterdb.push_back(&PARAM_CLUSTERSEARCH_MODE);
        clusterdb.push_back(&PARAM_FOLDSEEK_PATH);

        clusterSearchMode = 0;
        suboptHitsFactor = 0;
        filterSelfMatch = 0;
        maxGeneGaps = 3;
        clusterSize = 2;
        pMHThr = 0.01;
        pCluThr = 0.01;
        clusterUseWeight = 0;
        profileClusterSearch = 0;
        fileInclude = ".*";
        fileExclude = "^$";
        gffDir = "";
        std::string binaryDir = FileUtil::dirName(getAbsExePath());
        foldseekPath =  binaryDir.empty() ? "foldseek" : binaryDir + "/foldseek";

        //TODO: add citations (foldseek & mmseqs & clustersearch)
        citations.emplace(CITATION_SPACEDUST, "");
    }
private:
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif
