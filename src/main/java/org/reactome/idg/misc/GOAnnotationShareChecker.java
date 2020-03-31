package org.reactome.idg.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.reactome.data.GODataAnalyzerV2;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to check if two genes have shared go terms.
 * @author wug
 *
 */
public class GOAnnotationShareChecker extends FINetworkBuildConfig {
    
    private static final Logger logger = Logger.getLogger(GOAnnotationShareChecker.class);
    private GODataAnalyzerV2 goAnalyzer;
    
    public GOAnnotationShareChecker() {
        goAnalyzer = new GODataAnalyzerV2();
    }
    
    /**
     * Load all gene pairs that have domain-domain interactions
     * @return about 1.7 million pairs should be expected from this method.
     * @throws Exception
     */
    public Set<String> loadGenePairsViaGOBPShare() throws IOException {
        Map<String, String> geneToUniprot = ApplicationConfig.getConfig().getGeneToUniProMap();
        List<String> geneList = new ArrayList<>(geneToUniprot.keySet());
        Collections.sort(geneList);
        logger.info("Total genes to be checked: " + geneList.size());
        Map<String, Set<String>> proteinToGO = goAnalyzer.loadProteinToGOBPTerms();
        logger.info("Total proteins with GO BP annotated: " + proteinToGO.size());
        Set<String> rtn = new HashSet<>();
        for (int i = 0; i < geneList.size() - 1; i++) {
            String gene1 = geneList.get(i);
            String uniprot1 = geneToUniprot.get(gene1);
            for (int j = i + 1; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                String uniprot2 = geneToUniprot.get(gene2);
                // The order doesn't matter for this check
                boolean isInteracting = goAnalyzer.isTermShared(uniprot1 + "\t" + uniprot2,
                                                                proteinToGO);
                if (isInteracting)
                    rtn.add(gene1 + "\t" + gene2);
            }
        }
        logger.info("Total gene pairs having GO BP shared: " + rtn.size());
        return rtn;
    } 

}
