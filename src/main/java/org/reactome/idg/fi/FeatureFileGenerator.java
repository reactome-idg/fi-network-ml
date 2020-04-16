package org.reactome.idg.fi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.coexpression.CoExpressionLoader;
import org.reactome.idg.harmonizome.HarmonizomePairwiseLoader;
import org.reactome.idg.misc.GOAnnotationShareChecker;
import org.reactome.idg.misc.ProteinDDIChecker;
import org.reactome.idg.ppi.MappedPPIDataHandler;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to generate a feature file subject to ML.
 * @author wug
 *
 */
public class FeatureFileGenerator {
    private static final Logger logger = Logger.getLogger(FeatureFileGenerator.class);
    
    public FeatureFileGenerator() {
    }
    
    public static void main(String[] args) {
        if (args.length == 0) {
            System.err.println("java -Xmx16G -jar XXX {check_features|generate_matrix} {out_file}");
            System.exit(1);
        }
        if (args[0].equals("check_features")) {
            new FeatureFileGenerator().checkFeatures();
            return;
        }
        if (args[0].equals("generate_matrix")) {
            if (args.length == 1) {
                System.err.println("Please provide a file name for output.");
                System.exit(1);
            }
            try {
                new FeatureFileGenerator().buildFeatureMatrix(args[1]);
            }
            catch(Exception e) {
                logger.error(e.getMessage(), e);
            }
        }
    }
    
    /**
     * Dump the feature file into a tab-delimited matrix file.
     * Note: The output is "," delimited to avoid replace FI's tab!
     * @param outFileName
     * @throws Exception
     */
    public void buildFeatureMatrix(String outFileName) throws Exception {
        logger.info("Loading all features...");
        Map<String, Set<String>> featureToPairs = loadAllFeatures();
        logger.info("Feature loading is done. Total features: " + featureToPairs.size());
        // Check the features
        featureToPairs.forEach((feature, pair) -> {
            logger.info(feature + ": " + pair.size());
        });
        // Positive training data set
        Set<String> reactomeFIs = ApplicationConfig.getConfig().loadReactomeFIsInGenes();
        logger.info("Total Reactome FIs: " + reactomeFIs.size()); 
        // Filter FIs that don't have any positive feature since these FIs will not contribute
        // anything to the training
        boolean isValid = false;
        for (Iterator<String> it = reactomeFIs.iterator(); it.hasNext();) {
            String fi = it.next();
            isValid = false;
            for (String feature : featureToPairs.keySet()) {
                Set<String> pairs = featureToPairs.get(feature);
                if (pairs.contains(fi)) {
                    isValid = true;
                    break;
                }
            }
            if (isValid)
                continue;
            it.remove();
        }
        logger.info("Total Reactome FIs after filtering FIs having no feature: " + reactomeFIs.size());
        logger.info("Start dumping...");
        // Let's start dump
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        // Generate the header
        // Make sure we have a fixed order
        List<String> features = new ArrayList<>(featureToPairs.keySet());
        StringBuilder builder = new StringBuilder();
        builder.append("GenePair,FI");
        features.forEach(feature -> builder.append(",").append(feature));
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String fi : reactomeFIs) {
            builder.append(fi).append(",1"); // 1 for true, 0 for false
            for (String feature : features) {
                Set<String> pairs = featureToPairs.get(feature);
                builder.append(",").append(pairs.contains(fi) ? "1" : "0");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        // Random pairs as negative data set
        Set<String> geneIds = InteractionUtilities.grepIDsFromInteractions(reactomeFIs);
        String ratio = ApplicationConfig.getConfig().getAppConfig("negative.to.positive.ratio");
        if (ratio == null || ratio.length() == 0)
            ratio = "100";
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(geneIds,
                                                                           (int)(reactomeFIs.size() * Double.parseDouble(ratio)), 
                                                                           reactomeFIs);
        logger.info("Total random pairs as the negative dataset: " + randomPairs.size());
        // As noted in the original FINetworkContruction project, we will not do filtering for the negative
        // data set. (see in class NBCAnalyzer.java).
        for (String fi : randomPairs) {
            builder.append(fi).append(",0"); // 1 for true, 0 for false
            for (String feature : features) {
                Set<String> pairs = featureToPairs.get(feature);
                builder.append(",").append(pairs.contains(fi) ? "1" : "0");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
        logger.info("All done!");
    }
    
    public void checkFeatures() {
        try {
            logger.info("Check features...");
            long memory = Runtime.getRuntime().totalMemory();
            logger.info("Total memory before loading: " + memory / (1024 * 1024.0d) + " MB.");
            Map<String, Set<String>> feature2pairs = loadAllFeatures();
            logger.info("Total memory after loading all features: " + memory / (1024 * 1024.0d) + " MB.");
            // Want to print out all features
            for (String feature : feature2pairs.keySet()) {
                Set<String> pairs = feature2pairs.get(feature);
                logger.info(feature + ": " + pairs.size());
            }
        }
        catch(Exception e) {
            logger.error(e.getMessage(), e);
        }
    }
    
    /**
     * Load all used features into a key (feature name) to a set of pairwise relationships.
     * As of April 16, 2020, 106 features have been collected.
     * @throws IOException
     */
    public Map<String, Set<String>> loadAllFeatures() throws Exception {
        // We want to control the order of the insertion. Therefore, 
        // a LinkedHashMap, instead of a usual HashMap, is used here.
        Map<String, Set<String>> feature2pairs = new LinkedHashMap<>();
        loadPPIFeatures(feature2pairs);
        loadMiscFeatures(feature2pairs);
        // Used to sort files
        Comparator<File> fileSorter = (file1, file2) -> file1.getName().compareTo(file2.getName());
        loadHarmonizomeFeatures(feature2pairs, fileSorter);
        // Get the percentile for coexpression data
        String coexpPercentile = ApplicationConfig.getConfig().getAppConfig("coexpression.percentile");
        if (coexpPercentile == null || coexpPercentile.length() == 0)
            coexpPercentile = "0.001";
        logger.info("Coexpression precentile: " + coexpPercentile);
        double coexpPercentValue = new Double(coexpPercentile);
        // GTEx
        CoExpressionLoader coexpressionHandler = new CoExpressionLoader();
        List<File> gteFiles = coexpressionHandler.getGTExCoExpressionFiles();
        gteFiles.sort(fileSorter);
        logger.info("Loading GTEx features...");
        loadCoExpFeatures(coexpressionHandler,
                          gteFiles,
                          "GTEx",
                          coexpPercentValue,
                          feature2pairs);
        logger.info("GTEx features loading is done.");
        // TCGA
        List<File> tcgaFiles = coexpressionHandler.getTCGACoExpressionFiles();
        tcgaFiles.sort(fileSorter);
        logger.info("Loading TCGA features...");
        loadCoExpFeatures(coexpressionHandler,
                          gteFiles,
                          null, // Provide by the file name directly
                          coexpPercentValue,
                          feature2pairs);
        logger.info("TCGA features loading is done.");
        return feature2pairs;
    }
    
    private void loadCoExpFeatures(CoExpressionLoader loader,
                                   List<File> files,
                                   String featureType,
                                   double percentile,
                                   Map<String, Set<String>> feature2pairs) throws IOException {
        for (File file : files) {
            logger.info("Loading " + file.getName() + "...");
            String feature = file.getName();
            feature = feature.split("_")[0];
            if (featureType != null)
                feature = featureType + "-" + feature;
            Set<String> pairs = loader.loadCoExpressionViaPercentile(file, percentile);
            feature2pairs.put(feature, pairs);
            logger.info("Done.");
        }
    }

    private Comparator<File> loadHarmonizomeFeatures(Map<String, Set<String>> feature2pairs,
                                                     Comparator<File> fileSorter) throws Exception {
        logger.info("Loading harmonizome features...");
        HarmonizomePairwiseLoader harmonizomeHandler = new HarmonizomePairwiseLoader();
        Map<File, Double> file2percentile = harmonizomeHandler.getSelectedDownloadFiles();
        List<File> files = file2percentile.keySet()
                .stream()
                .sorted(fileSorter)
                .collect(Collectors.toList());
        for (File file : files) {
            logger.info("Loading " + file.getName() + "...");
            Double percentile = file2percentile.get(file);
            // Get the feature name from the file name
            String feature = file.getName();
            feature = feature.split("\\.")[0]; // We only need the first part as our feature name
            Set<String> pairs = harmonizomeHandler.loadPairwisesFromDownload(file, percentile);
            feature2pairs.put(feature, pairs);
            logger.info("Done.");
        }
        logger.info("Harmonizome features loading is done.");
        return fileSorter;
    }

    private void loadMiscFeatures(Map<String, Set<String>> feature2pairs) throws IOException {
        // Domain interaction
        logger.info("Loading domain-domain interactions...");
        ProteinDDIChecker ddiHandler = new ProteinDDIChecker();
        Set<String> ddis = ddiHandler.loadGenePairsViaDDIs();
        logger.info("Done.");
        feature2pairs.put("DomainInteractions", ddis);
        // GO BP sharing
        logger.info("Loading GO BO sharing...");
        GOAnnotationShareChecker goHandler = new GOAnnotationShareChecker();
        Set<String> goPairs = goHandler.loadGenePairsViaGOBPShare();
        feature2pairs.put("GOBPSharing", goPairs);
        logger.info("Done.");
    }

    private void loadPPIFeatures(Map<String, Set<String>> feature2pairs) throws IOException {
        // PPI first
        MappedPPIDataHandler ppiHandler = new MappedPPIDataHandler();
        logger.info("Loading HumanPPIs...");
        Set<String> humanPPIs = ppiHandler.loadHumanPPIs();
        feature2pairs.put("HumanPPI", humanPPIs);
        logger.info("Done.");
        logger.info("Loading MousePPIs...");
        Set<String> mousePPIs = ppiHandler.loadMousePPIs();
        feature2pairs.put("MousePPI", mousePPIs);
        logger.info("Done.");
        logger.info("Loading FlyPPIs...");
        Set<String> flyPPIs = ppiHandler.loadFlyPPIs();
        feature2pairs.put("FlyPPI", flyPPIs);
        logger.info("Done.");
        logger.info("Loading WormPPIs...");
        Set<String> wormPPIs = ppiHandler.loadWormPPIs();
        feature2pairs.put("WormPPI", wormPPIs);
        logger.info("Done.");
        logger.info("Loading YeastPPIs...");
        Set<String> yeastPPIs = ppiHandler.loadYeastPPIs();
        feature2pairs.put("YeatPPI", yeastPPIs);
        logger.info("Done.");
    }
    
    
}
