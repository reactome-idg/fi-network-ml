package org.reactome.idg.harmonizome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;

/**
 * This class is used to process a download similarity data and generate a set of selected
 * pairs based on some threshold values.
 * @author wug
 *
 */
public class DataProcessor {
    private final static Logger logger = LogManager.getLogger(DataProcessor.class);
    public static final int FIRST_INDEX = 3;
    public static final double VALUE_THRESHOLD = 0.50d;
    public static final int GENE_NUMBER_THRESHOLD = 12000; // Make sure a similarity file has more than 12,000 genes

    // Keep the target genes we want to investiage
    private Set<String> allGenes;
    
    public DataProcessor() {
    }
    
    @SuppressWarnings("unchecked")
    private Set<String> getAllGenes() throws Exception {
        if (allGenes != null)
            return allGenes;
        // Need to load all human genes from a Reactome database
        InputStream is = getClass().getClassLoader().getResourceAsStream("db.properties");
        Properties prop = new Properties();
        prop.load(is);
        is.close();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            prop.getProperty("reactome.db.name"), 
                                            prop.getProperty("mysql.user"), 
                                            prop.getProperty("mysql.password"));
        GKInstance human = dba.fetchInstance(48887L);
        Collection<GKInstance> refGenes = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                       ReactomeJavaConstants.species,
                                                                       "=",
                                                                       human);
        allGenes = new HashSet<>();
        for (GKInstance refGene : refGenes) {
            String gene = (String) refGene.getAttributeValue(ReactomeJavaConstants.geneName);
            if (gene != null)
                allGenes.add(gene);
        }
        logger.info("Total human genes from Reactome: " + allGenes.size());
        return allGenes;
    }

    /**
     * Process the download similarity file and generate two files: processed for
     * a smaller zipped file with decimal downed to .3f, filtered for a pair-wise
     * similarity file after thresholding. The src file will be checked. If the total
     * number of genes in src is less than the GENE_NUMBER_THRESHOLD, the file will be not
     * processed. Also it is possible that a similarity file doesn't provide any score
     * higher (absolute value) than the threshold, then no file will be generated.
     * Note: Don't use Scanner, which is much slower than using FileReader and BufferedReader.
     * Also don't zip on the fly, which is also slower than a pure write. We should zip
     * all files after the processing.
     * @param src
     * @param processed 
     * @param filterd
     * @throws FileNotFoundException
     */
    public boolean processCorrelations(File src,
                                       File processed,
                                       File filtered) throws Exception {
        FileReader fileReader = new FileReader(src);
        BufferedReader br = new BufferedReader(fileReader);
        String line = br.readLine();
        int currentRow = 0;
        String[] geneHeaders = line.split("\t");
        Set<String> allGenes = getAllGenes();
        int totalGenes = getTotalGenes(geneHeaders, allGenes);
        if (totalGenes < GENE_NUMBER_THRESHOLD) {
            br.close();
            logger.info("Not enough genes in the similariy file: " + src.getName());
            return false; // We don't want this download
        }
        // Need to escape two more line
        for (int i = 1; i < FIRST_INDEX; i++) {
            br.readLine();
            currentRow ++;
        }
        long time1 = System.currentTimeMillis();
        PrintWriter processedWriter = new PrintWriter(processed);
        PrintWriter filteredWriter = new PrintWriter(filtered);
        StringBuilder builder = new StringBuilder();
        // Generate the header
        builder.append("Gene");
        for (int i = FIRST_INDEX; i < geneHeaders.length; i++) {
            if (!allGenes.contains(geneHeaders[i]))
                continue;
            builder.append("\t").append(geneHeaders[i]);
        }
        processedWriter.println(builder.toString());
        builder.setLength(0);
        String[] tokens = null;
        boolean hasValue = false;
        while ((line = br.readLine()) != null) {
            currentRow ++;
//            if (currentRow == 10000)
//                break;
            tokens = line.split("\t");
            if (!allGenes.contains(tokens[0]))
                continue;
            builder.append(tokens[0]);
            // All empty cells
            for (int i = 0; i < currentRow; i++)
                builder.append("\t");
            // The matrix is square, therefore, we only need to handle the top-right half of the matrix
            for (int i = currentRow + 1; i < tokens.length; i++) {
                if (!allGenes.contains(tokens[0]))
                    continue;
                Double value = new Double(tokens[i]);
                if (Math.abs(value) >= VALUE_THRESHOLD) { // Make sure it is absolute
                    String fi = generateFI(tokens[0], geneHeaders[i]);
                    filteredWriter.println(fi + "\t" + (value > 0 ? "+" : "-"));
                    hasValue = true;
                }
                builder.append("\t").append(tokens[i]);
                // Don't bother to format. We should reduce the file size a lot.
                // Format itself costs a lot of time!
//                builder.append("\t").append(String.format("%.3f", value));
            }
            processedWriter.println(builder.toString());
            builder.setLength(0);
        }
        br.close();
        filteredWriter.close();
        processedWriter.close();
        long time2 = System.currentTimeMillis();
        logger.info("Time for processing correlations: " + (time2 - time1) / (60.0d * 1000) + " minutes.");
        if (!hasValue) {
            // Delete all files
            processed.delete();
            filtered.delete();
            logger.info("No value passed the threshold: " + src.getName());
            return false;
        }
        return true;
    }
    
    private String generateFI(String gene1, String gene2) {
        int compare = gene1.compareTo(gene2);
        if (compare < 0)
            return gene1 + "\t" + gene2;
        else
            return gene2 + "\t" + gene1;
    }

    private int getTotalGenes(String[] geneHeaders, Set<String> needed) {
        int count = 0;
        for (int i = FIRST_INDEX; i < geneHeaders.length; i++) {
            if (!needed.contains(geneHeaders[i]))
                continue;
            count ++;
        }
        return count;
    }

}
