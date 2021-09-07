package org.reactome.idg.fi;

import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.annotations.FeatureDesc;
import org.reactome.idg.annotations.FeatureLoader;
import org.reactome.idg.util.ApplicationConfig;

public class FeatureFileGeneratorTests {
    
    public FeatureFileGeneratorTests() {
    }
    
    @Test
    public void checkPredictedPartners() throws Exception {
        double scoreCutoff = 0.80d;
        String dirName = "results/feature_files/prediction/";
        String scoreFile = dirName + "prd_probabilities_prediction_061820_w_ppi_names.csv";
        String gene = "PRKY";
        
        FileUtility fu = new FileUtility();
        fu.setInput(scoreFile);
        String line = fu.readLine();
        Set<String> partners = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            if (!line.contains(gene))
                continue;
            String[] tokens = line.split(",");
            if (Double.parseDouble(tokens[4]) < scoreCutoff)
                continue;
            String[] genes = tokens[0].split("_");
            partners.addAll(Arrays.asList(genes));
        }
        partners.remove(gene);
        System.out.println("Total partners: " + partners.size());
        partners.stream().sorted().forEach(System.out::println);
    }
    
    /**
     * Print out a feature list for one gene or a pair of genes.
     * @throws Exception
     */
    @Test
    public void checkFeatures() throws Exception {
        double scoreCutoff = 0.50d;
        String dirName = "results/feature_files/prediction/";
        String featureFileName = dirName + "prediction_061820.csv";
        String scoreFile = dirName + "prd_probabilities_prediction_061820_w_ppi_names.csv";
        String gene = "PRKY";
        gene = "SBK2";
        gene = "CSNK2A3";
        gene = "ABL1";
        
//        gene = "CACNG6";
//        scoreCutoff = 0.65d;
        
        String outName = dirName + gene + "_prediction_061820.csv";
        // Get a list of pairs having scores > scoreCutoff
        Map<String, Double> pairToScore = new HashMap<>();
        FileUtility fu = new FileUtility();
        fu.setInput(scoreFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            if (!line.contains(gene))
                continue;
            String[] tokens = line.split(",");
            if (Double.parseDouble(tokens[4]) < scoreCutoff)
                continue;
            pairToScore.put(tokens[0].replace('_', '\t'), new Double(tokens[4]));
        }
        fu.close();
        System.out.println("Done loading scores: " + pairToScore.size() + " pairs.");
        Set<String> querySet = pairToScore.keySet();
        fu.setInput(featureFileName);
        fu.setOutput(outName);
        line = fu.readLine();
        String[] tokens = line.split(",");
        StringBuilder builder = new StringBuilder();
        builder.append(tokens[0]);
        builder.append(",RF_Score");
        for (int i = 1; i < tokens.length; i++)
            builder.append(",").append(tokens[i]);
        fu.printLine(builder.toString());
        builder.setLength(0);
        while ((line = fu.readLine()) != null) {
            tokens = line.split(",");
            if (!querySet.contains(tokens[0])) {
                continue;
            }
            builder.append(tokens[0]).append(",").append(pairToScore.get(tokens[0]));
            for (int i = 1; i < tokens.length; i++)
                builder.append(",").append(tokens[i]);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkPredictionResults() throws Exception {
        // This file contains predicted results for pairs based on trained RF
        String fileName = "results/feature_files/prediction/prd_probabilities_prediction_061820_w_ppi_names.csv";
        int totalPairs = 0;
        int selectedPairs = 0;
        double threshold = 0.67d; // For at least two features
//        threshold = 0.90d; // To control precision at 0.50
        threshold = 0.80; // The point the recall curve and the precision curve cross.
        threshold = 0.85;
        threshold = 0.875;
        threshold = 0.86124659652243;
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            Double score = new Double(tokens[4]);
            if (score >= threshold)
                selectedPairs ++;
            totalPairs ++;
            max = Math.max(max, score);
            min = Math.min(min, score);
        }
        fu.close();
        System.out.println("Max score: " + max);
        System.out.println("Min score: " + min);
        System.out.println("Threshold: " + threshold);
        System.out.println("Total pairs: " + totalPairs);
        System.out.println("Selected pairs: " + selectedPairs);
        System.out.println("Ratio: " + (double)selectedPairs / totalPairs);
    }
   
    @Test
    public void checkNonReactomeFIs() throws Exception {
        Set<String> nonReactomeFIs = ApplicationConfig.getConfig().loadNonReactomeFIsInGenes();
        System.out.println("Non Reactome FIs: " + nonReactomeFIs.size());
    }
    
    @Test
    public void checkPairsInTrainingData() throws Exception {
        String fileName = "results/feature_files/feature_matrix_041720.csv";
        Set<String> excludedPairs = null;
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            excludedPairs = lines.skip(1)
                                 .map(line -> line.split(",")[0])
                                 .collect(Collectors.toSet());
        }
        System.out.println("Total excluded pairs: " + excludedPairs.size());
        Set<String> nonReactomeFIs = ApplicationConfig.getConfig().loadNonReactomeFIsInGenes();
        System.out.println("Total non-Reactome FIs: " + nonReactomeFIs.size());
        nonReactomeFIs.removeAll(excludedPairs);
        System.out.println("Total non-reactome FIs after removing pairs in the training: " + nonReactomeFIs.size());
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(nonReactomeFIs);
        System.out.println("Total genes in non-Reactome FIs: " + genes.size());
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(genes, nonReactomeFIs.size() * 100, nonReactomeFIs);
        System.out.println("Random pairs generated: " + randomPairs.size());
        randomPairs.removeAll(excludedPairs);
        System.out.println("Random pairs after removing pairs in the training: " + randomPairs.size());
    }
    
    @Test
    public void checkAnnotations() throws Exception {
        Method[] methods = FeatureFileGenerator.class.getDeclaredMethods();
        for (Method method : methods) {
//            System.out.println(method.getName());
            FeatureLoader featureLoader = method.getAnnotation(FeatureLoader.class);
            if (featureLoader != null) {
                System.out.println(method.getName());
                String[] loaderMethods = featureLoader.methods();
                for (String methodName : loaderMethods) {
                    String[] tokens = methodName.split(",");
                    String loaderMethod = null;
                    if (tokens.length == 2)
                        loaderMethod = tokens[1];
                    else
                        loaderMethod = tokens[0];
                    // Need to extra parameters
                    int index = loaderMethod.lastIndexOf(".");
                    Class<?> cls = Class.forName(loaderMethod.substring(0, index));
                    Method clsMethod = searchMethod(loaderMethod.substring(index + 1), cls);
                    System.out.println("\t" + clsMethod.getName());
                    // Expect a FeatureDesc annotation
                    FeatureDesc desc = clsMethod.getAnnotation(FeatureDesc.class);
                    String sourceText = Arrays.asList(desc.sources())
                            .stream()
                            .map(s -> s.toString())
                            .collect(Collectors.joining(", "));
                    System.out.println("\t\t" + desc.type() + ": " + sourceText);
                }
            }
        }
    }
    
    private <T> Method searchMethod(String methodName, Class<T> cls) {
        Method[] methods = cls.getDeclaredMethods();
        for (Method method : methods) {
            if (method.getName().equals(methodName) && method.getAnnotation(FeatureDesc.class) != null)
                return method;
        }
        return null;
    }

}
