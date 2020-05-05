package org.reactome.idg.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.reactome.data.PfamAnalyzer;
import org.reactome.idg.annotations.FeatureDesc;
import org.reactome.idg.model.FeatureSource;
import org.reactome.idg.model.FeatureType;
import org.reactome.idg.util.ApplicationConfig;
import org.reactome.idg.util.FINetworkBuildConfig;

/**
 * In order to use this class, a pre-generated Uni2Pfam.txt should exist.
 * See the documentation for building the FI network for details on how to
 * generate this file.
 * @author wug
 *
 */
public class ProteinDDIChecker extends FINetworkBuildConfig {
    private static final Logger logger = Logger.getLogger(ProteinDDIChecker.class);
    private PfamAnalyzer pfamAnayzer;

    public ProteinDDIChecker() {
        pfamAnayzer = new PfamAnalyzer();
        pfamAnayzer.setpFamDirName(ApplicationConfig.getConfig().getAppConfig("pfam.dir"));
        pfamAnayzer.setUniprotDirName(ApplicationConfig.getConfig().getAppConfig("uniprot.dir"));
    }

    /**
     * Load all gene pairs that have domain-domain interactions
     * @return about 1.7 million pairs should be expected from this method.
     * @throws Exception
     */
    @FeatureDesc(sources = {FeatureSource.pFam},
                 type = FeatureType.Domain_Interaction)
    public Set<String> loadGenePairsViaDDIs() throws IOException {
        Map<String, String> geneToUniprot = ApplicationConfig.getConfig().getGeneToUniProMap();
        List<String> geneList = new ArrayList<>(geneToUniprot.keySet());
        Collections.sort(geneList);
        Set<String> rtn = new HashSet<>();
        logger.info("Total genes to be checked: " + geneList.size());
        for (int i = 0; i < geneList.size() - 1; i++) {
            String gene1 = geneList.get(i);
            String uniprot1 = geneToUniprot.get(gene1);
            for (int j = i + 1; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                String uniprot2 = geneToUniprot.get(gene2);
                boolean isInteracting = pfamAnayzer.checkIfInteracting(uniprot1, uniprot2);
                if (isInteracting)
                    rtn.add(gene1 + "\t" + gene2);
            }
        }
        logger.info("Total gene pairs having domain-domain interactions: " + rtn.size());
        return rtn;
    } 

}
