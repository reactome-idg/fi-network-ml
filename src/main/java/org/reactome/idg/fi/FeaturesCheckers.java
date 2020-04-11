package org.reactome.idg.fi;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.idg.coexpression.CoExpressionLoader;
import org.reactome.idg.harmonizome.HarmonizomePairwiseLoader;
import org.reactome.idg.misc.GOAnnotationShareChecker;
import org.reactome.idg.misc.ProteinDDIChecker;
import org.reactome.idg.ppi.MappedPPIDataHandler;
import org.reactome.idg.util.ApplicationConfig;
import org.reactome.idg.util.FINetworkBuildConfig;

/**
 * This class is used to handle all features collected for FI ML prediction.
 * @author wug
 *
 */
public class FeaturesCheckers extends FINetworkBuildConfig {
    private final Logger logger = Logger.getLogger(FeaturesCheckers.class);
    private FeatureChecker checker; // The actual object used to check features.
    
    public FeaturesCheckers() {
    }
    
    public static void main(String[] args) {
        FeaturesCheckers handler = new FeaturesCheckers();
        try {
            handler.checkPPIFeatures();
            handler.checkMiscFeatures();
            handler.checkHarmonizomeFeatures();
            handler.checkGeneExpressionFeatures();
        }
        catch(Exception e) {
            handler.logger.error(e.getMessage(), e);
        }
    }
    
    private void setUpChecker() throws IOException {
        if (checker == null) {
            checker = new FeatureChecker();
            Set<String> reactomeFIs = ApplicationConfig.getConfig().loadReactomeFIsInGenes();
            checker.setInteractionSet(reactomeFIs);
        }
    }
    
    private void checkFeatureOddsRatio(Set<String> pairwiseRels) throws Exception {
        setUpChecker();
        checker.checkFeatureOddsRatio(pairwiseRels);
    }
    
    /**
     * Check PPI features provided via package org.reactome.idg.ppi.
     * @throws Exception
     */
    public void checkPPIFeatures() throws Exception {
        MappedPPIDataHandler ppiHandler = new MappedPPIDataHandler();
        Set<String> humanPPIs = ppiHandler.loadHumanPPIs();
        logger.info("Check human PPIs...");
        checkFeatureOddsRatio(humanPPIs);
        Set<String> mousePPIs = ppiHandler.loadMousePPIs();
        logger.info("Check PPIs mapped from mouse...");
        checkFeatureOddsRatio(mousePPIs);
        Set<String> flyPPIs = ppiHandler.loadFlyPPIs();
        logger.info("Check PPIs mapped from fly...");
        checkFeatureOddsRatio(flyPPIs);
        Set<String> wormPPIs = ppiHandler.loadWormPPIs();
        logger.info("Check PPIs mapped from worm...");
        checkFeatureOddsRatio(wormPPIs);
        Set<String> yeastPPIs = ppiHandler.loadYeastPPIs();
        logger.info("Check PPIs mapped from yeast...");
        checkFeatureOddsRatio(yeastPPIs);
    }
    
    /**
     * Check other features provided via package org.reactome.idg.misc
     * @throws Exception
     */
    public void checkMiscFeatures() throws Exception {
        ProteinDDIChecker ddiChecker = new ProteinDDIChecker();
        Set<String> ddis = ddiChecker.loadGenePairsViaDDIs();
        logger.info("Check domain-domain interactions...");
        checkFeatureOddsRatio(ddis);
        GOAnnotationShareChecker goBPChecker = new GOAnnotationShareChecker();
        Set<String> bpShared = goBPChecker.loadGenePairsViaGOBPShare();
        logger.info("Check BP sharing...");
        checkFeatureOddsRatio(bpShared);
    }
    
    /**
     * Check gene coexpression features provided via package org.reactome.idg.coexpression.
     * @throws Exception
     */
    public void checkGeneExpressionFeatures() throws Exception {
        String cutoff = ApplicationConfig.getConfig().getAppConfig("coexpression.cutoff");
        if (cutoff == null || cutoff.length() == 0)
            cutoff = "0.5d";
        logger.info("Coexpression cutoff: " + cutoff);
        CoExpressionLoader loader = new CoExpressionLoader();
        List<File> files = loader.getGTExCoExpressionFiles();
        logger.info("Total GTEx files: " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            Set<String> rels = loader.loadCoExpression(file, new Double(cutoff));
            checkFeatureOddsRatio(rels);
        }
    }
    
    /**
     * Check gene similarity features provided via package org.reactome.idg.harmonizome
     * @throws Exception
     */
    public void checkHarmonizomeFeatures() throws Exception {
        HarmonizomePairwiseLoader loader = new HarmonizomePairwiseLoader();
        List<File> files = loader.getPairwiseFiles();
        logger.info("Total Harmonizome files: " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            Set<String> rels = loader.loadPairwises(file);
            checkFeatureOddsRatio(rels);
        }
    }

}
