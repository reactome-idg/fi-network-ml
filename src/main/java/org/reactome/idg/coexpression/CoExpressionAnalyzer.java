package org.reactome.idg.coexpression;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

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
    
    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.err.println("Provide parameters: data_source_dir, plot_output_dir, {debug}");
            System.exit(1);
        }
        File sourceDir = new File(args[0]);
        File resultDir = new File(args[1]);
        File[] files = sourceDir.listFiles();
        CoExpressionAnalyzer analyzer = new CoExpressionAnalyzer();
        for (File file : files) {
            String fileName = file.getName();
            if (!fileName.endsWith("_Spearman_Adj.csv"))
                continue;
            logger.info("Processing file " + file.getAbsolutePath());
            File resultFile = new File(resultDir, fileName.split("\\.")[0] + ".pdf");
            analyzer.plotDistributions(file.getAbsolutePath(),
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
        FileUtility fu = new FileUtility();
        fu.setInput(dataFileName);
        // Temp file
        File file = new File(dataFileName);
        // Create a temp file for collected values
        file = new File(file.getName() + ".tmp");
        fu.setOutput(file.getAbsolutePath());
        fu.printLine("Coexpression");
        String line = fu.readLine();
        Random random = new Random();
        int c = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            int index = random.nextInt(tokens.length - 1);
            String token = tokens[index + 1];
            if (token.equals("NA"))
                continue; // Remove NA
            fu.printLine(tokens[index + 1]);
            c ++;
        }
        fu.close();
        logger.info("Collected data points from " + dataFileName + ": " + c);
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

}
