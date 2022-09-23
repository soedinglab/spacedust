#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H
 
#include <Parameters.h>

const int CITATION_CLUSTERSEARCH = CITATION_END;

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
    std::vector<MMseqsParameter*> iterativeclusearchworkflow;
    std::vector<MMseqsParameter*> besthitbyset;
    std::vector<MMseqsParameter*> combinehits;
    std::vector<MMseqsParameter*> clusterhits;
    std::vector<MMseqsParameter*> foldseeksearch;

    PARAMETER(PARAM_FILTER_SELF_MATCH)
    PARAMETER(PARAM_MULTIHIT_PVAL)
    PARAMETER(PARAM_CLUSTER_PVAL)
    PARAMETER(PARAM_MAX_GENE_GAP)
    PARAMETER(PARAM_CLUSTER_SIZE)
    PARAMETER(PARAM_CLUSTER_USE_WEIGHT)
    PARAMETER(PARAM_PROFILE_CLUSTER_SEARCH)
    PARAMETER(PARAM_GFF_DIR)


    float pMHThr;
    float pCluThr;
    bool clusterUseWeight;
    bool profileClusterSearch;
    int maxGeneGaps;
    int clusterSize;
    std::string gffDir;
    bool filterSelfMatch;
    
private:
    LocalParameters() : 
        Parameters(),
        PARAM_FILTER_SELF_MATCH(PARAM_FILTER_SELF_MATCH_ID, "--filter-self-match", "Filter self match", "Remove hits between the same set. 0: do not filter, 1: filter", typeid(bool), (void *) &filterSelfMatch, ""),
        PARAM_MULTIHIT_PVAL(PARAM_MULTIHIT_PVAL_ID, "--multihit-pval", "Multihit P-value cutoff", "Multihit P-value threshold for cluster match output", typeid(float), (void *) &pMHThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_CLUSTER_PVAL(PARAM_CLUSTER_PVAL_ID, "--cluster-pval", "Clustering and Ordering P-value cutoff","Clustering and Ordering P-value threshold for cluster match output", typeid(float), (void *) &pCluThr, "^0(\\.[0-9]+)?|^1(\\.0+)?$"),
        PARAM_MAX_GENE_GAP(PARAM_MAX_GENE_GAP_ID, "--max-gene-gap", "Maximum gene gaps", "Maximum number of genes allowed between two clusters to merge", typeid(int), (void *) &maxGeneGaps, "^[1-9]{1}[0-9]*$"),
        PARAM_CLUSTER_SIZE(PARAM_CLUSTER_SIZE_ID, "--cluster-size", "Minimal cluster size", "Minimum number of genes to define cluster", typeid(int), (void *) &clusterSize, "^[1-9]{1}[0-9]*$"),
        PARAM_CLUSTER_USE_WEIGHT(PARAM_CLUSTER_USE_WEIGHT_ID, "--cluster-use-weight", "Cluster weighting factor ", "Use weighting factor to calculate cluster match score", typeid(bool), (void *) &clusterUseWeight, ""),
        PARAM_PROFILE_CLUSTER_SEARCH(PARAM_PROFILE_CLUSTER_SEARCH_ID, "--profile-cluster-search", "Cluster search against profiles", "Perform profile(target)-sequence searches in clustersearch", typeid(bool), (void *) &profileClusterSearch, ""),
        PARAM_GFF_DIR(PARAM_GFF_DIR_ID, "--gff-dir", "gff dir file", "Path to gff dir file", typeid(std::string), (void *) &gffDir, "")
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

        // foldseeksearch
        foldseeksearch.push_back(&PARAM_E);
        foldseeksearch.push_back(&PARAM_ADD_BACKTRACE);
        foldseeksearch.push_back(&PARAM_C);
        foldseeksearch.push_back(&PARAM_COV_MODE);
        foldseeksearch.push_back(&PARAM_THREADS);
        foldseeksearch.push_back(&PARAM_COMPRESSED);
        foldseeksearch.push_back(&PARAM_V);

        // multi hit db
        createsetdb = combineList(createdb, extractorfs);
        createsetdb = combineList(createsetdb, translatenucs);
        createsetdb = combineList(createsetdb, gff2db);
        createsetdb = combineList(createsetdb, aa2foldseek);
        createsetdb.push_back(&PARAM_GFF_DIR);

        // multi hit search
        clustersearchworkflow = combineList(searchworkflow, besthitbyset);
        clustersearchworkflow = combineList(clustersearchworkflow, combinehits);
        clustersearchworkflow = combineList(clustersearchworkflow, clusterhits);
        clustersearchworkflow.push_back(&PARAM_PROFILE_CLUSTER_SEARCH);

        //aa2foldseek
        aa2foldseek = combineList(prefilter, align);
        aa2foldseek.push_back(&PARAM_COMPRESSED);
        aa2foldseek.push_back(&PARAM_THREADS);
        aa2foldseek.push_back(&PARAM_V);

        // iterative cluster search
        iterativeclusearchworkflow = combineList(searchworkflow, besthitbyset);
        iterativeclusearchworkflow = combineList(iterativeclusearchworkflow, mergedbs);
        iterativeclusearchworkflow = combineList(iterativeclusearchworkflow, subtractdbs);
        iterativeclusearchworkflow = combineList(iterativeclusearchworkflow, combinehits);
        iterativeclusearchworkflow = combineList(iterativeclusearchworkflow, clusterhits);


        filterSelfMatch = 0;
        maxGeneGaps = 3;
        clusterSize = 2;
        pMHThr = 0.01;
        pCluThr = 0.01;
        clusterUseWeight = 1;
        profileClusterSearch = 0;
        gffDir = "";

        
        citations.emplace(CITATION_CLUSTERSEARCH, "");
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif
