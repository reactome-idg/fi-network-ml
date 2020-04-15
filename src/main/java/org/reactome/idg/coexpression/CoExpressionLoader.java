package org.reactome.idg.coexpression;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to load coexpression data.
 * @author wug
 *
 */
public class CoExpressionLoader {
    private static final Logger logger = Logger.getLogger(CoExpressionLoader.class);
    
    public CoExpressionLoader() {
    }

    private List<File> getGeneCoExpressionFiles(File dir) {
        return Arrays.asList(dir.listFiles())
                     .stream()
                     .filter(file -> file.getName().endsWith("_Spearman_Adj.csv"))
                     .collect(Collectors.toList());
    }
    
    public List<File> getGTExCoExpressionFiles() {
        String dir = ApplicationConfig.getConfig().getAppConfig("gtex.coexpression.dir");
        if (dir == null || dir.length() == 0)
            return new ArrayList<>();
        return getGeneCoExpressionFiles(new File(dir));
    }
    
    public List<File> getTCGACoExpressionFiles() {
        String dir = ApplicationConfig.getConfig().getAppConfig("tcga.coexpression.dir");
        if (dir == null || dir.length() == 0)
            return new ArrayList<>();
        return getGeneCoExpressionFiles(new File(dir));
    }
    
    public List<File> getTCGAFilesFromList() throws IOException {
        String dir = ApplicationConfig.getConfig().getAppConfig("tcga.coexpression.dir");
        if (dir == null || dir.length() == 0)
            return new ArrayList<>();
        // Used to handle some files that have not be executed
        File file = new File("TCGA_CoExpression_Files_List.txt");
        if (!file.exists()) return new ArrayList<>();
        try (Stream<String> lines = Files.lines(Paths.get(file.getAbsolutePath()))) {
            return lines.filter(line -> !line.startsWith("#"))
                        .map(line -> new File(dir, line))
                        .collect(Collectors.toList());
        }
    }
    
    /**
     * Get the absolute coexpression threshold value for the suggested top percentile (e.g.
     * 0.001).
     * @param ratio
     * @return
     * @throws IOException
     */
    public double getCutoffValueForRatio(File file, double percentile) throws IOException {
        // Try to load all values into a list
        List<Float> values = new ArrayList<>(); // Use Float to control memory usage.
        FileUtility fu = new FileUtility();
        fu.setInput(file.getAbsolutePath());
        String line = fu.readLine();
        int c = 1;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            // We want to look at the top-right triangle values
            for (int i = c + 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA") || tokens[i].equals("TRUE")) {
                    //                        System.err.println(line);
                    continue; // Just ignore NA
                }
                Float value = Float.parseFloat(tokens[i]);
                values.add(Math.abs(value));
            }
            c ++;
        }
        fu.close();
        logger.info("All values have been loaded for " + file.getName() + ": " + values.size());
        values.sort(Comparator.reverseOrder());
        logger.info("Sorting is done.");
        int index = (int) (values.size() * percentile);
        return values.get(index);
    }
    
    /**
     * Some of the code below was copied from org.reactome.idg.pairwise.main.GTExDataProcessor.java.
     * @param file
     * @return
     * @throws IOException
     */
    public Set<String> loadCoExpression(File file, double cutoff) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(file.getAbsolutePath());
        String line = fu.readLine();
        Set<String> rels = new HashSet<>();
        String[] genes = line.split(",");
        // The first token will be empty
        List<String> geneList = Arrays.asList(genes);
        int c = 1;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            // We want to look at the top-right triangle values
            for (int i = c + 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA") || tokens[i].equals("TRUE")) {
                    //                        System.err.println(line);
                    continue; // Just ignore NA
                }
                Double value = new Double(tokens[i]);
                if (Math.abs(value) > cutoff) {
                    String gene1 = tokens[0];
                    String gene2 = geneList.get(i);
                    if (gene1.equals(gene2))
                        throw new IllegalStateException("Gene1 and Gene2 should not be the same: " + gene1);
                    String rel = InteractionUtilities.generateFIFromGene(gene1, gene2);
                    rels.add(rel);
                }
            }
            c ++;
        }
        fu.close();
        return rels;
    }
    
}