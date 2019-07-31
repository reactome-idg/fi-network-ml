package org.reactome.idg.harmonizome;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import org.reactome.idg.model.Gene;
import org.reactome.idg.service.GeneCorrelationService;

/**
 * This class is used to load a download correlation matrix file into a local mysql database.
 * Beware that the data file may be huge for all human genes, over 4G. So it is strongly suggested
 * not to store all files locally.
 * @author wug
 *
 */
public class DataLoader {
    private static final int FIRST_INDEX = 3;
    
    public DataLoader() {
    }
    
    /**
     * This is the actual method to load the data file into a local database.
     * @param fileName
     * @param service
     * @throws Exception
     */
    public void loadCorrelation(String fileName, 
                                GeneCorrelationService service) throws FileNotFoundException {
        File file = new File(fileName);
        Scanner scanner = new Scanner(file);
        // This is the header. We need this.
        String line = scanner.nextLine();
        String[] geneHeaders = line.split("\t");
        Map<String, Gene> nameToGene = ensureGenesInPersistance(geneHeaders, service);
        scanner.close();
    }
    
    /**
     * Make sure all genes in the header are in the database. 
     * @param geneHeaders
     * @param service
     * @throws Exception
     */
    private Map<String, Gene> ensureGenesInPersistance(String[] geneHeaders, 
                                                       GeneCorrelationService service) {
        Map<String, Gene> nameToGene = new HashMap<>();
        for (int i = FIRST_INDEX; i < geneHeaders.length; i++) {
            Gene gene = service.fetchGene(geneHeaders[i]);
            nameToGene.put(geneHeaders[i], gene);
        }
        return nameToGene;
    }

}
