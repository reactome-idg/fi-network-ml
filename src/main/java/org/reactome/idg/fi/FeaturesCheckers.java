package org.reactome.idg.fi;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;
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
//            handler.checkPPIFeatures();
//            handler.checkMiscFeatures();
//            handler.checkHarmonizomeFeatures();
//            handler.checkGeneExpressionFeatures();
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
    public void checkGeneExpressionFeaturesViaCutoff() throws Exception {
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
        files = loader.getTCGACoExpressionFiles();
        logger.info("Total TCGA files: " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            Set<String> rels = loader.loadCoExpression(file, new Double(cutoff));
            checkFeatureOddsRatio(rels);
        }
    }
    
    /**
     * Check gene coexpression features provided via package org.reactome.idg.coexpression.
     * @throws Exception
     */
    public void checkGeneExpressionFeatures() throws Exception {
        String percentile = ApplicationConfig.getConfig().getAppConfig("coexpression.percentile");
        if (percentile == null || percentile.length() == 0)
            percentile = "0.001";
        logger.info("Coexpression precentile: " + percentile);
        CoExpressionLoader loader = new CoExpressionLoader();
        List<File> files = loader.getGTExCoExpressionFiles();
        logger.info("Total GTEx files: " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            double cutoff = loader.getCutoffValueForRatio(file, new Double(percentile));
            logger.info("Found cutoff: " + cutoff);
            Set<String> rels = loader.loadCoExpression(file, cutoff);
            checkFeatureOddsRatio(rels);
        }
        files = loader.getTCGACoExpressionFiles();
//        files = loader.getTCGAFilesFromList();
        logger.info("Total TCGA files: " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            double cutoff = loader.getCutoffValueForRatio(file, new Double(percentile));
            logger.info("Found cutoff: " + cutoff);
            Set<String> rels = loader.loadCoExpression(file, cutoff);
            checkFeatureOddsRatio(rels);
        }
    }
    
    /**
     * Check gene similarity features provided via package org.reactome.idg.harmonizome
     * @throws Exception
     */
    public void checkHarmonizomeFeatures() throws Exception {
        String percentile = ApplicationConfig.getConfig().getAppConfig("harmonizome.percentile");
        if (percentile == null || percentile.length() == 0)
            percentile = "0.01"; // Default for harmonizome is 0.01
        logger.info("Chosen percentile: " + percentile);
        HarmonizomePairwiseLoader loader = new HarmonizomePairwiseLoader();
//        List<File> files = loader.getPairwiseFiles();
//        List<File> files = loader.getProcessedFiles();
        List<File> files = loader.getDownloadedFiles();
//        List<File> files = loader.getSelectedDownloadFiles();
        logger.info("Total Harmonizome files (selected downloaded): " + files.size());
        for (File file : files) {
            logger.info("Check " + file.getName() + "...");
            Set<String> rels = loader.loadPairwisesFromDownload(file, new Double(percentile));
            checkFeatureOddsRatio(rels);
        }
    }
    
    @Test
    public void collectResults() throws IOException {
        String dir = "results/features_check/";
        String[] files = {"out_041120.txt", "out_041220.txt"};
        files = new String[] {
//                "out_041320.txt",
//                "out_041420_tcga.txt",
//                "out_041420_harmonizome.txt",
                "out_041420_harmonizome_1.txt"
        };
        FileUtility fu = new FileUtility();
        StringBuilder builder = new StringBuilder();
        String line = null;
        System.out.println("Feature\tTotalPairs\tCutoff\tFilteredPairs\tMappedPairs\tRatio\tOddsRatio\tOR_SD");
        Double cutoff = null;
        for (String file : files) {
            fu.setInput(dir + file);
            while ((line = fu.readLine()) != null) {
//                System.out.println(line);
                if (line.contains("org.reactome.idg.fi.FeaturesCheckers  - Check")) {
                    if (builder.length() > 0) {
                        System.out.println(builder.toString());
                        builder.setLength(0);
                        cutoff = null;
                    }
                    int index = line.lastIndexOf("Check");
                    int index1 = line.lastIndexOf("...");
                    String feature = line.substring(index + "Check".length(), index1).trim();
                    builder.append(feature);
                }
                else if (line.contains("Total values:") ||
                         line.contains("All values have been loaded")) {
                    int index = line.lastIndexOf(":");
                    String totalValues = line.substring(index + 1).trim();
                    builder.append("\t").append(totalValues);
                }
                else if (line.contains("Cutoff value: ") ||
                         line.contains("Found cutoff:")) {
                    int index = line.lastIndexOf(":");
                    cutoff = new Double(line.substring(index + 1).trim()); // Hold on for the time being
                }
                else if (line.contains("Cutoff adjusted to:")) {
                    int index = line.lastIndexOf(":");
                    cutoff = new Double(line.substring(index + 1).trim());
                }
                else if (line.startsWith("Total checked pairs:")) {
                    int index = line.indexOf(":");
                    String totalPairs = line.substring(index + 1).trim();
                    // Don't forget to push this first
                    builder.append("\t").append(cutoff);
                    builder.append("\t").append(totalPairs);
                }
                else if (line.startsWith("Mapped to ppi:")) {
                    String[] tokens = line.split(": | \\(|\\)");
                    builder.append("\t").append(tokens[1]);
                    builder.append("\t").append(tokens[2].substring(1, tokens[2].length() - 1));
                }
                else if (line.startsWith("Average odds ratio:")) {
                    String[] tokens = line.split(": | \\(| \\+- ");
                    builder.append("\t").append(tokens[1]);
                    builder.append("\t").append(tokens[2]);
                }
            }
            System.out.println(builder.toString());
            fu.close();
        }
    }
    
    @Test
    public void test() {
        String line = "Mapped to ppi: 31612 (0.222111)";
        String[] tokens = line.split(": | \\(|\\)");
        Arrays.asList(tokens).stream().forEach(System.out::println);
    }

}
