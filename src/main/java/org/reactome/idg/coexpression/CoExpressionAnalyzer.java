package org.reactome.idg.coexpression;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.r3.util.ProcessRunner;

/**
 * This class is used to do some analysis on the gene co-expression data generated from the TCGA 
 * and GTE expression data.
 * @author wug
 *
 */
public class CoExpressionAnalyzer {
    private static final Logger logger = Logger.getLogger(CoExpressionAnalyzer.class);
    private String rScriptCommand = "Rscript"; // Default value
    private final String rScript = "RSrc/GeneCoExpressionPlotter.R";
    private final int TOTAL_SAMPLE_SIZE = 20000;
    
    public static void main(String[] args) throws IOException {
        CoExpressionAnalyzer analyzer = new CoExpressionAnalyzer();
        if (args[0].equals("plotDistributions"))
            analyzer.plotDistributions(Arrays.copyOfRange(args, 1, args.length));
        else if (args[0].equals("performKSTests"))
            analyzer.performKSTests(Arrays.copyOfRange(args, 1, args.length));
    }
    
    /**
     * The args should an array of co-expression data directories.
     * @param args
     * @throws IOException
     */
    public void performKSTests(String[] args) throws IOException {
        List<File> files = Arrays.asList(args)
                                .stream()
                                .map(arg -> new File(arg))
                                .flatMap(dir -> Arrays.asList(dir.listFiles()).stream())
                                .filter(file -> file.getName().endsWith("Spearman_Adj.csv"))
                                .collect(Collectors.toList());
        // For performance reason, we need to sample values first and then cache them
        Map<File, double[]> fileToValues = new HashMap<>();
        for (File file : files) {
            logger.info("Sampling " + file.getAbsolutePath() + "...");
            double[] values = sampleCoExpression(file.getAbsolutePath(),
                                                 TOTAL_SAMPLE_SIZE);
            fileToValues.put(file, values);
        }
        FileUtility fu = new FileUtility();
        fu.setOutput("PairwiseKSTests.txt");
        fu.printLine("File1\tFile2\tKS-PValue");
        for (int i = 0; i < files.size(); i++) {
            File file1 = files.get(i);
            double[] values1 = fileToValues.get(file1);
            // Include itself for positive control
            for (int j = i; j < files.size(); j++) {
                File file2 = files.get(j);
                double[] values2 = fileToValues.get(file2);
                double ks = TestUtils.kolmogorovSmirnovTest(values1, values2);
                logger.info(file1.getName() + "\t" + 
                             file2.getName() + "\t" + 
                             ks);
                fu.printLine(file1.getName() + "\t" + 
                             file2.getName() + "\t" + 
                             ks);
            }
        }
        fu.close();
    }

    private void plotDistributions(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Provide parameters: plotDistributions(must have) data_source_dir, plot_output_dir, {debug}");
            System.exit(1);
        }
        File sourceDir = new File(args[0]);
        File resultDir = new File(args[1]);
        File[] files = sourceDir.listFiles();
        for (File file : files) {
            String fileName = file.getName();
            if (!fileName.endsWith("_Spearman_Adj.csv"))
                continue;
            logger.info("Processing file " + file.getAbsolutePath());
            File resultFile = new File(resultDir, fileName.split("\\.")[0] + ".pdf");
            plotDistributions(file.getAbsolutePath(),
                                       resultFile.getAbsolutePath());
            logger.info("Done.");
            if (args.length > 2 && args[2].equals("debug")) {
                logger.info("Processing one file only for debug");
                break;
            }
        }
    }
    
    public CoExpressionAnalyzer() {
    }
    
    /**
     * This method will call an R script to use ggplot2 to draw distribtion.
     * @param dataFileName
     * @param outFileName
     * @throws IOException
     */
    public void plotDistributions(String dataFileName,
                                  String outFileName) throws IOException {
        double[] sampleValues = sampleCoExpression(dataFileName, TOTAL_SAMPLE_SIZE);
        FileUtility fu = new FileUtility();
        // Temp file
        File file = new File(dataFileName);
        // Create a temp file for collected values
        file = new File(file.getName() + ".tmp");
        fu.setOutput(file.getAbsolutePath());
        fu.printLine("Coexpression");
        for (Double value : sampleValues)
            fu.printLine(value + "");
        fu.close();
        // Get the tile from the data file name
        File dataFile = new File(dataFileName);
        String title = dataFile.getName().split("_")[0];
        plotInR(file.getAbsolutePath(),
                outFileName,
                title);
        file.delete();
    }
    
    private void plotInR(String dataFileName,
                         String outFileName,
                         String title) throws IOException {
        ProcessRunner runner = new ProcessRunner();
        List<String> parameters = new ArrayList<String>();
        parameters.add(rScriptCommand);
        parameters.add(rScript);
        parameters.add(dataFileName);
        parameters.add(outFileName);
        parameters.add(title);
        String[] output = runner.runScript(parameters.toArray(new String[0]));
        if (output != null && output.length > 0) {
            for (int i = 0; i < output.length; i++) {
                if (output[i].trim().length() == 0)
                    continue;
                logger.info("Output " + i + ": " + output[i]);
            }
        }
    }
    
    public double performKSTest(String fileName) throws IOException {
        return performKSTest(fileName, fileName);
    }
    
    public double performKSTest(String fileName1, String fileName2) throws IOException {
        double[] values1 = sampleCoExpression(fileName1, TOTAL_SAMPLE_SIZE);
        double[] values2 = sampleCoExpression(fileName2, TOTAL_SAMPLE_SIZE);
        return TestUtils.kolmogorovSmirnovTest(values1, values2);
    }
    
    private double[] sampleCoExpression(String fileName,
                                        int numberOfValues) throws IOException {
        //        logger.info("Loading values in " + fileName + "...");
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split(",");
        long totalNumber = (tokens.length - 1) * (tokens.length - 2) / 2;
        double ratio = (double) numberOfValues / totalNumber;
        int lineIndex = 1;
        List<Double> values = new ArrayList<>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split(",");
            for (int i = lineIndex + 1; i < tokens.length; i++) {
                // For some reason, there is "TRUE" in the data set
                if (tokens[i].equals("NA") || tokens[i].equals("TRUE"))
                    continue;
                if (Math.random() < ratio) // Check if we need this value
                    values.add(new Double(tokens[i]));
            }
            lineIndex ++;
        }
        //        logger.info("Total loaded double values: " + values.size());
        fu.close();
        return values.stream()
                .mapToDouble(d -> d)
                .toArray();
    }
    
    @Test
    public void testPerformKSTest() throws IOException {
        String dirName = "/Users/wug/git/reactome-idg/idg-pairwise/examples";
        String fileName = dirName + File.separator + "TCGA-LIHC_Spearman_Adj.csv";
        // Let's run 10 times
        System.out.println("KS test:");
        for (int i = 0; i < 100; i++) {
            System.out.println(i + "\t" + performKSTest(fileName));
        }
    }
    
    @Test
    public void testPerformKSTest2() throws IOException {
        String dirName = "/Users/wug/git/reactome-idg/data-importer/results/coexpression";
        String fileName1 = dirName + File.separator + "Breast-MammaryTissue_Spearman_Adj.csv";
        String fileName2 = dirName + File.separator + "TCGA-BRCA_Spearman_Adj.csv";
        System.out.println("KS test:");
        for (int i = 0; i < 5; i++) {
            System.out.println(i + "\t" + performKSTest(fileName1, fileName2));
        }
    }
    
    @Test
    public void testPlotDistributions() throws IOException {
        String dirName = "/Users/wug/git/reactome-idg/idg-pairwise/examples";
        String fileName = dirName + File.separator + "TCGA-LIHC_Spearman_Adj.csv";
        String outFileName = dirName + File.separator + "Test.pdf";
        rScriptCommand = "/usr/local/bin/Rscript";
        File file = new File(rScript);
        System.out.println(file.getAbsolutePath());
        if (file.exists()) {
            System.out.println("R script exists!");
        }
        else
            System.out.println("R script doesn't exist!");
        plotDistributions(fileName, outFileName);
    }
    
    /**
     * This method is used to check the outliers flagged by customized R function, 
     * https://github.com/reactome-idg/gather-app/blob/master/python/scripts/R/functions.R,
     * and excluded in adjacency matrix calculation. The function is based on the first
     * dimension (PC1) of the loading matrix, which is generated by prcomp with samples as
     * columns. This approach may be arguable! 
     * @throws IOException
     */
    @Test
    public void checkExcludedOutliers() throws IOException {
        String dirName = "results/eda/";
        String[] sources = {"gtex", "tcga"};
        for (String source : sources) {
            System.out.println("Data source: " + source);
            File[] files = new File(dirName + source).listFiles();
            System.out.println("File\tTotal\tOutliers\tPercentage");
            for (File file : files) {
                if (!file.getName().endsWith(".csv"))
                    continue;
                int total = 0;
                int outliers = 0;
                List<String> lines = Files.readAllLines(Paths.get(file.getAbsolutePath()));
                for (int i = 1; i < lines.size(); i++) {
                    String line = lines.get(i); // First line is header
                    total ++;
                    if (line.split(",")[1].equals("\"1\""))
                        outliers ++;
                }
                String percentage = String.format("%.3f", (double) outliers / total);
                System.out.println(file.getName() + "\t" + total + "\t" + outliers + "\t" + percentage);
            }
            System.out.println();
        }
    }
    

}
