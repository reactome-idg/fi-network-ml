package org.reactome.idg.fi;

import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.annotations.FeatureDesc;
import org.reactome.idg.annotations.FeatureLoader;
import org.reactome.idg.util.ApplicationConfig;

public class FeatureFileGeneratorTests {
    
    public FeatureFileGeneratorTests() {
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
