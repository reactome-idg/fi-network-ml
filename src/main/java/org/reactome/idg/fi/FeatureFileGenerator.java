package org.reactome.idg.fi;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.annotations.FeatureLoader;
import org.reactome.idg.coexpression.CoExpressionLoader;
import org.reactome.idg.harmonizome.HarmonizomePairwiseLoader;
import org.reactome.idg.misc.GOAnnotationShareChecker;
import org.reactome.idg.misc.ProteinDDIChecker;
import org.reactome.idg.model.FeatureSource;
import org.reactome.idg.ppi.MappedPPIDataHandler;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to generate a feature file subject to ML.
 * @author wug
 *
 */
public class FeatureFileGenerator {
    private static final Logger logger = Logger.getLogger(FeatureFileGenerator.class);
    // Control if some features should be generated according postive and negative
    private boolean needNegative = false;
    // Originally Harmonizome- is not added to features for ML. However, for the 
    // idg.reactome.org, we need to add this. Therefore, we have this flag
    private boolean prefixHarmonizomeInFeature = false;
    
    public FeatureFileGenerator() {
    }
    
    /**
     * This method is used to generate a file having one feature is positive only.
     * @throws IOException
     */
    @Test
    public void generateSingleFeatureFile() throws IOException {
        String featureFile = "results/feature_files/test/feature_test_matrix_051120.csv";
        String outFile = "results/feature_files/prediction/single_feature_file_060320.csv";
        FileUtility fu = new FileUtility();
        fu.setInput(featureFile);
        fu.setOutput(outFile);
        String line = fu.readLine();
        String[] tokens = line.split(",");
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < tokens.length; i++) {
            if (i == 1)
                continue;
            builder.append(tokens[i]).append(",");
        }
        builder.deleteCharAt(builder.length() - 1);
        fu.printLine(builder.toString());
        builder.setLength(0);
        int length = tokens.length - 2; // Don't need the first and second columns
        for (int i = 0; i < length; i++) {
            builder.append("FI" + i);
            for (int j = 0; j < length; j++) {
                builder.append(",");
                if (i == j)
                    builder.append("1");
                else
                    builder.append("0");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    public static void main(String[] args) {
        if (args.length == 0) {
            System.err.println("java -Xmx48G -jar XXX {check_features|generate_matrix|generate_test_matrix|generate_prediction_file} {out_file} {training_file}");
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
            return;
        }
        if (args[0].equals("generate_test_matrix")) {
            if (args.length < 3) {
                System.err.println("Please provide a file name for output and the file name for"
                        + " the training data generated from the generate_matrix running.");
                System.exit(1);
            }
            try {
                new FeatureFileGenerator().buildFeatureMatrixForTest(args[2], args[1]);
            }
            catch(Exception e) {
                logger.error(e.getMessage(), e);
            }
            return;
        }
        // Generate a file used for prediction. FIs extracted from Reactome and other
        // pathway databases are not excluded in the final output so that they can be used
        // as an internal control.
        if (args[0].equals("generate_prediction_file")) {
            if (args.length < 2) {
                System.err.println("Please provide a file name for output.");
                System.exit(1);
            }
            try {
                new FeatureFileGenerator().generatePredictionFile(args[1]);
            }
            catch(Exception e) {
                logger.error(e.getMessage(), e);
            }
        }
    }
    
    public void generatePredictionFile(String outFileName) throws Exception {
        logger.info("Loading all features...");
        Map<String, Set<String>> feature2pairs = loadAllFeatures();
        logger.info("All features have been loaded.");
        // Need to collect all genes in the features
        List<String> allGenes = feature2pairs.values()
                                             .stream()
                                             .flatMap(v -> v.stream())
                                             .map(pair -> pair.split("\t"))
                                             .flatMap(genes -> Stream.of(genes))
                                             .distinct() // Avoid duplication
                                             .sorted()
                                             .collect(Collectors.toList());
        logger.info("Total genes collected from all features: " + allGenes.size());
        logger.info("Starting generating the predict file...");
        // Pair-wise generations
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("GenePair");
        List<String> features = new ArrayList<>(feature2pairs.keySet());
        features.stream().forEach(feature -> builder.append(",").append(feature));
        fu.printLine(builder.toString());
        builder.setLength(0);
        // All genes have been sorted
        boolean isNeeded = false;
        int count = 0;
        for (int i = 0; i < allGenes.size() - 1; i++) {
            String gene1 = allGenes.get(i);
            for (int j = i + 1; j < allGenes.size(); j++) {
                String gene2 = allGenes.get(j);
                String pair = gene1 + "\t" + gene2;
                builder.append(pair);
                isNeeded = false;
                for (String feature : features) {
                    Set<String> pairs = feature2pairs.get(feature);
                    builder.append(",");
                    if (pairs.contains(pair)) {
                        isNeeded = true;
                        builder.append("1");
                    }
                    else
                        builder.append("0");
                }
                // We will collect genes having at least one feature
                if (isNeeded) {// This check should be reliable
                    fu.printLine(builder.toString());
                    count ++;
                }
                builder.setLength(0);
            }
        }
        fu.close();
        logger.info("Done. Output in file " + outFileName);
        logger.info("Total pairs having at least one feature: " + count);
    }
    
    public boolean isNeedNegative() {
        return needNegative;
    }

    public void setNeedNegative(boolean needNegative) {
        this.needNegative = needNegative;
    }

    public boolean isPrefixHarmonizomeInFeature() {
        return prefixHarmonizomeInFeature;
    }

    public void setPrefixHarmonizomeInFeature(boolean prefixHarmonizomeInFeature) {
        this.prefixHarmonizomeInFeature = prefixHarmonizomeInFeature;
    }

    /**
     * Use this method to create an independent test data set based on non-Reactome FIs.
     * @param trainingFileName: the training data file name used to exclude pairs there.
     * @param outFileName
     * @throws Exception
     */
    public void buildFeatureMatrixForTest(String trainingFileName,
                                          String outFileName) throws Exception {
        // Load all pairs in the training data set. Since random pairs are used as negative,
        // different training data set may be different.
        Set<String> excludedPairs = null;
        try (Stream<String> lines = Files.lines(Paths.get(trainingFileName))) {
            excludedPairs = lines.skip(1)
                                 .map(line -> line.split(",")[0])
                                 .collect(Collectors.toSet());
        }
        logger.info("Total pairs that loaded from the training data and will be excluded: " + excludedPairs.size());
        // The following steps are very similar to ones used to generate the training dataset.
        // The only difference is that pairs from the training data sets will be removed.
        Map<String, Set<String>> feature2pairs = loadAllFeatures();
        Set<String> nonReactomeFIs = ApplicationConfig.getConfig().loadNonReactomeFIsInGenes();
        logger.info("Total non-Reactome FIs: " + nonReactomeFIs.size());
        buildFeatureMatrix(nonReactomeFIs,
                           excludedPairs, 
                           feature2pairs,
                           outFileName);
        logger.info("All done. The output is in: " + outFileName);
    }
    
    /**
     * Dump the feature file into a tab-delimited matrix file. This method is used to 
     * create the training data set.
     * Note: The output is "," delimited to avoid replace FI's tab!
     * @param outFileName
     * @throws Exception
     */
    public void buildFeatureMatrix(String outFileName) throws Exception {
        Map<String, Set<String>> featureToPairs = loadAllFeatures();
        // Positive training data set
        Set<String> reactomeFIs = ApplicationConfig.getConfig().loadReactomeFIsInGenes();
        logger.info("Total Reactome FIs: " + reactomeFIs.size()); 
        buildFeatureMatrix(reactomeFIs, 
                           new HashSet<>(),
                           featureToPairs, 
                           outFileName);
        logger.info("All done!");
    }

    private void buildFeatureMatrix(Set<String> fis, 
                                   Set<String> toBeExcluded, // FIs in this set should not be used for both positive and negative sets
                                   Map<String, Set<String>> featureToPairs, 
                                   String outFileName) throws IOException {
        logger.info("Total FIs passed into the method: " + fis.size());
        boolean isChanged = fis.removeAll(toBeExcluded);
        logger.info("Filter by removing in toBeExcluded: " + fis.size());
        if (isChanged) {
            // Check left gene ids
            Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
            logger.info("Total ids left after removing excluded pairs: " + genes.size());
        }
        filterFIsToOneFeatureMinimum(featureToPairs, fis);
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
        for (String fi : fis) {
            builder.append(fi).append(",1"); // 1 for true, 0 for false
            for (String feature : features) {
                Set<String> pairs = featureToPairs.get(feature);
                builder.append(",").append(pairs.contains(fi) ? "1" : "0");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        // Random pairs as negative data set
        Set<String> geneIds = InteractionUtilities.grepIDsFromInteractions(fis);
        String ratio = ApplicationConfig.getConfig().getAppConfig("negative.to.positive.ratio");
        if (ratio == null || ratio.length() == 0)
            ratio = "100";
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(geneIds,
                                                                           (int)(fis.size() * Double.parseDouble(ratio)), 
                                                                           fis);
        logger.info("Total random pairs as the negative dataset: " + randomPairs.size());
        randomPairs.removeAll(toBeExcluded);
        logger.info("Total random pairs after removing pairs in tobeRemoved: " + randomPairs.size());
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
    }

    private void filterFIsToOneFeatureMinimum(Map<String, Set<String>> featureToPairs, Set<String> reactomeFIs) {
        logger.info("Total Fis before filtering: " + reactomeFIs.size());
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
        logger.info("Total FIs after filtering FIs having no feature: " + reactomeFIs.size());
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
        logger.debug("Loading all features...");
        long memory = Runtime.getRuntime().totalMemory();
        logger.debug("Total memory before loading all features: " + memory / (1024 * 1024.0d) + " MB.");
        // We want to control the order of the insertion. Therefore, 
        // a LinkedHashMap, instead of a usual HashMap, is used here.
        Map<String, Set<String>> feature2pairs = new LinkedHashMap<>();
        loadPPIFeatures(feature2pairs);
        loadMiscFeatures(feature2pairs);
        // Used to sort files
        Comparator<File> fileSorter = getFileSorter();
        loadHarmonizomeFeatures(feature2pairs, fileSorter);
        Double coexpPercentValue = getCoExpressionPercentile();
        // GTEx
        loadGTExCoExpressions(feature2pairs, fileSorter, coexpPercentValue);
        // TCGA
        loadTCGACoExpressions(feature2pairs, fileSorter, coexpPercentValue);
        logger.info("Feature loading is done. Total features: " + feature2pairs.size());
        memory = Runtime.getRuntime().totalMemory();
        logger.debug("Total memory after loading all features: " + memory / (1024 * 1024.0d) + " MB.");
        // Check the features
        feature2pairs.forEach((feature, pair) -> {
            logger.info(feature + ": " + pair.size());
        });
        return feature2pairs;
    }

    private Comparator<File> getFileSorter() {
        return (file1, file2) -> file1.getName().compareTo(file2.getName());
    }

    public Double getCoExpressionPercentile() {
        // Get the percentile for coexpression data
        String coexpPercentile = ApplicationConfig.getConfig().getAppConfig("coexpression.percentile");
        if (coexpPercentile == null || coexpPercentile.length() == 0)
            coexpPercentile = "0.001";
        logger.info("Coexpression precentile: " + coexpPercentile);
        return new Double(coexpPercentile);
    }
    
    public void loadTCGACoExpressions(Map<String, Set<String>> feature2pairs) throws IOException {
        loadTCGACoExpressions(feature2pairs, null, null);
    }

    @FeatureLoader(methods = {"org.reactome.idg.coexpression.CoExpressionLoader.loadCoExpressionViaPercentile"},
                   source = FeatureSource.TCGA)
    private void loadTCGACoExpressions(Map<String, Set<String>> feature2pairs,
                                      Comparator<File> fileSorter,
                                      Double coexpPercentValue) throws IOException {
        if (coexpPercentValue == null)
            coexpPercentValue = getCoExpressionPercentile();
        fileSorter = fileSorter == null ? getFileSorter() : fileSorter;
        CoExpressionLoader coexpressionHandler = new CoExpressionLoader();
        coexpressionHandler.setNeedNegative(needNegative);
        List<File> tcgaFiles = coexpressionHandler.getTCGACoExpressionFiles();
        tcgaFiles.sort(fileSorter);
        logger.info("Loading TCGA features...");
        loadCoExpFeatures(coexpressionHandler,
                          tcgaFiles,
                          null, // Provide by the file name directly
                          coexpPercentValue,
                          feature2pairs);
        logger.info("TCGA features loading is done.");
    }
    
    public void loadGTExCoExpressions(Map<String, Set<String>> feature2pairs) throws IOException {
        loadGTExCoExpressions(feature2pairs, null, null);
    }
    
    @FeatureLoader(methods = {"org.reactome.idg.coexpression.CoExpressionLoader.loadCoExpressionViaPercentile"},
                   source = FeatureSource.GTEx)
    private void loadGTExCoExpressions(Map<String, Set<String>> feature2pairs,
                                      Comparator<File> fileSorter,
                                      Double coexpPercentValue) throws IOException {
        if (coexpPercentValue == null)
            coexpPercentValue = getCoExpressionPercentile();
        fileSorter = fileSorter == null ? getFileSorter() : fileSorter;
        // GTEx
        CoExpressionLoader coexpressionHandler = new CoExpressionLoader();
        coexpressionHandler.setNeedNegative(needNegative);
        List<File> gteFiles = coexpressionHandler.getGTExCoExpressionFiles();
        gteFiles.sort(fileSorter);
        logger.info("Loading GTEx features...");
        loadCoExpFeatures(coexpressionHandler,
                          gteFiles,
                          "GTEx",
                          coexpPercentValue,
                          feature2pairs);
        logger.info("GTEx features loading is done.");
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
    
    public void loadHarmonizomeFeatures(Map<String, Set<String>> feature2pairs) throws Exception {
        loadHarmonizomeFeatures(feature2pairs, null);
    }

    @FeatureLoader(methods = {"org.reactome.idg.harmonizome.HarmonizomePairwiseLoader.loadPairwisesFromDownload"},
                   source = FeatureSource.Harmonizome)
    private void loadHarmonizomeFeatures(Map<String, Set<String>> feature2pairs,
                                         Comparator<File> fileSorter) throws Exception {
        logger.info("Loading harmonizome features...");
        fileSorter = fileSorter == null ? getFileSorter() : fileSorter;
        HarmonizomePairwiseLoader harmonizomeHandler = new HarmonizomePairwiseLoader();
        harmonizomeHandler.setNeedNegative(needNegative);
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
            // Make sure the feature name starting with Harmonizome to downstream analysis
            if (prefixHarmonizomeInFeature)
                feature = "Harmonizome-" + feature;
            feature2pairs.put(feature, pairs);
            logger.info("Done.");
        }
        logger.info("Harmonizome features loading is done.");
    }

    @FeatureLoader(methods= {"DomainInteractions,org.reactome.idg.misc.ProteinDDIChecker.loadGenePairsViaDDIs",
                             "GOBPSharing,org.reactome.idg.misc.GOAnnotationShareChecker.loadGenePairsViaGOBPShare"})
    public void loadMiscFeatures(Map<String, Set<String>> feature2pairs) throws IOException {
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

    @FeatureLoader(methods= {"HumanPPI,org.reactome.idg.ppi.MappedPPIDataHandler.loadHumanPPIs",
                             "MousePPI,org.reactome.idg.ppi.MappedPPIDataHandler.loadMousePPIs",
                             "FlyPPI,org.reactome.idg.ppi.MappedPPIDataHandler.loadFlyPPIs",
                             "WormPPI,org.reactome.idg.ppi.MappedPPIDataHandler.loadWormPPIs",
                             "YeastPPI,org.reactome.idg.ppi.MappedPPIDataHandler.loadYeastPPIs"})
    public void loadPPIFeatures(Map<String, Set<String>> feature2pairs) throws IOException {
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
        feature2pairs.put("YeastPPI", yeastPPIs);
        logger.info("Done.");
    }
    
}
