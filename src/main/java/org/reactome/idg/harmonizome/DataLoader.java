package org.reactome.idg.harmonizome;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;
import org.reactome.idg.model.Provenance;
import org.reactome.idg.service.GeneCorrelationService;

/**
 * This class is used to load a downloaded correlation matrix file into a local mysql database.
 * Beware that the data file may be huge for all human genes, over 4G. So it is strongly suggested
 * not to store all files locally.
 * @author wug
 *
 */
public class DataLoader {
    private static final int FIRST_INDEX = 3;
    private static final int BATCH_SIZE = 10000; // Dont' assign too large
    private static final double MIN_VALUE = 0.001d;
    
    public DataLoader() {
    }
    
    /**
     * This method is used to check if all Provenance objects used by the Harmonizome
     * are in the database. If not, they will be saved into the database.
     */
    public Map<String, Provenance> ensureProvenancesInPersistance(GeneCorrelationService service) throws Exception {
        List<Provenance> provenances = new DataDownloader().getProvenaces();
        Map<String, Provenance> datasetToProvenance = new HashMap<>();
        for (Provenance provenance : provenances) {
            Provenance persisted = service.fetchProvenance(provenance);
            datasetToProvenance.put(persisted.getName(), persisted);
        }
        return datasetToProvenance;
    }
    
    /**
     * In some files two genes may have different cases, but pointing to the same thing.
     * e.g. C1ORF21 and C1orf21 in file cheappi.
     * @param tokens
     * @return
     */
    private Set<String> getGeneNamesForEscaped(String[] tokens) {
        Set<String> escaped = new HashSet<>();
        for (int i = FIRST_INDEX; i < tokens.length; i++) {
            String gene = tokens[i];
            if (gene.equals(gene.toUpperCase()))
                continue;
            escaped.add(gene);
        }
        return escaped;
    }
    
    /**
     * This is the actual method to load the data file into a local database.
     * @param fileName
     * @param service
     * @throws Exception
     */
    public void loadCorrelation(String fileName, 
                                Provenance provenance,
                                GeneCorrelationService service) throws FileNotFoundException {
        File file = new File(fileName);
        Scanner scanner = new Scanner(file);
        // This is the header. We need this.
        String line = scanner.nextLine();
        int currentRow = 0;
        String[] geneHeaders = line.split("\t");
        Set<String> escaped = getGeneNamesForEscaped(geneHeaders);
        Map<String, Gene> nameToGene = ensureGenesInPersistance(geneHeaders,
                                                                escaped,
                                                                provenance,
                                                                service);
        // Need to escape two more line
        for (int i = 1; i < FIRST_INDEX; i++) {
            scanner.nextLine();
            currentRow ++;
        }
        String[] tokens = null;
        List<GenePairCorrelation> batch = new ArrayList<GenePairCorrelation>(BATCH_SIZE);
        long time1 = System.currentTimeMillis();
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();
            currentRow ++;
            tokens = line.split("\t");
            if (escaped.contains(tokens[0]))
                continue;
            // The matrix is square, therefore, we only need to handle the top-right half of the matrix
            for (int i = currentRow + 1; i < tokens.length; i++) {
                Double value = new Double(tokens[i]);
                if (Math.abs(value) < MIN_VALUE)
                    continue;
                if (escaped.contains(geneHeaders[i]))
                    continue;
                GenePairCorrelation correlation = createCorrelation(geneHeaders[i],
                                                                    tokens[0],
                                                                    new Double(tokens[i]),
                                                                    nameToGene,
                                                                    provenance);
//                service.saveCorrelation(correlation);
                batch.add(correlation);
                if (batch.size() == BATCH_SIZE) {
                    service.saveCorrelations(batch);
                    batch.clear();
                }
            }
//            if (currentRow == 10)
//                break;
        }
        if (batch.size() > 0)
            service.saveCorrelations(batch);
        scanner.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for loading correlations: " + (time2 - time1) / (60.0d * 1000));
    }
    
    private GenePairCorrelation createCorrelation(String name1,
                                                  String name2, 
                                                  double corr,
                                                  Map<String, Gene> nameToGene,
                                                  Provenance provenance) {
        GenePairCorrelation corrObj = new GenePairCorrelation();
        corrObj.setCorrelationValue(corr);
        corrObj.setProvenance(provenance);
        int compare = name1.compareTo(name2);
        if (compare < 0) {
            corrObj.setGene1(nameToGene.get(name1));
            corrObj.setGene2(nameToGene.get(name2));
        }
        else {
            corrObj.setGene1(nameToGene.get(name2));
            corrObj.setGene2(nameToGene.get(name1));
        }
        return corrObj;
    }
    //TODO: SOme files have duplicated genes with different cases. This caused duplicated key issues
    // Need to filter out all genes having lower cases!!!
    /**
     * Make sure all genes in the header are in the database. 
     * @param geneHeaders
     * @param service
     * @throws Exception
     */
    private Map<String, Gene> ensureGenesInPersistance(String[] geneHeaders, 
                                                       Set<String> escaped,
                                                       Provenance provenance,
                                                       GeneCorrelationService service) {
        Map<String, Gene> nameToGene = new HashMap<>();
        for (int i = FIRST_INDEX; i < geneHeaders.length; i++) {
            if (escaped.contains(geneHeaders[i]))
                continue;
            Gene gene = service.fetchGene(geneHeaders[i]);
            gene.addProvenance(provenance);
            service.updateGene(gene);
            nameToGene.put(geneHeaders[i], gene);
        }
        return nameToGene;
    }

}
