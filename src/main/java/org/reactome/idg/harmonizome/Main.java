package org.reactome.idg.harmonizome;

import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.reactome.idg.config.AppConfig;
import org.reactome.idg.model.Provenance;
import org.reactome.idg.service.GeneCorrelationService;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class Main {

    private static final Logger logger = Logger.getLogger(Main.class.getName());

    /**
     * Entrypoint to code that imports Harmonizome data. 
     * This build this application (using the HarmonizomeImporter profile):
     *  mvn clean package -DskipTests -PHarmonizomeImporter
     * to run: 
     *  java -jar HarmonizomeImporter-jar-with-dependencies.jar
     */
    public static void main(String[] args) {
        AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext(AppConfig.class);
        GeneCorrelationService service = context.getBean(GeneCorrelationService.class);
        
        DataLoader loader = new DataLoader();
        String dirName = "datasets/Harmonizome/download/gene_similarity_matrix_cosine";
        // Hard-coded for test
        String fileName = "cheappi.txt"; // 4.2 G
        String dataset = "CHEA Transcription Factor Targets";
//        String fileName = "achilles.txt"; // 219 MB
//        String dataset = "Achilles Cell Line Gene Essentiality Profiles";
//        String fileName = "clinvar.txt"; // 54.4MB
//        String dataset = "ClinVar Gene-Phenotype Associations";
        logger.info("Starting loading " + fileName + "...");
        long current = System.currentTimeMillis();
        try {
            Map<String, Provenance> datasetToProvenance = loader.ensureProvenancesInPersistance(service);
            Provenance provenance = datasetToProvenance.get(dataset);
            if (provenance == null) {
                context.close();
                throw new IllegalStateException("No Provenance availble for dataset, " + dataset + "!");
            }
            loader.loadCorrelation(dirName + "/" + fileName,
                                   provenance,
                                   service);
            context.close();
            long time = System.currentTimeMillis() - current;
            logger.info("Done: " + time + " millseconds.");
        }
        catch(Exception e) {
            if (context.isActive())
                context.close();
            logger.log(Level.SEVERE, "Exception " + e.getMessage(), e);
        }
    }

}
