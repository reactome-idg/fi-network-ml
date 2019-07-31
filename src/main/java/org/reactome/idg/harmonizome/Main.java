package org.reactome.idg.harmonizome;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.reactome.idg.config.AppConfig;
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
        String fileName = "achilles.txt"; 
        logger.info("Starting loading " + fileName + "...");
        long current = System.currentTimeMillis();
        try {
            loader.loadCorrelation(dirName + "/" + fileName,
                                   service);
            context.close();
            long time = System.currentTimeMillis() - current;
            logger.info("Done: " + time + " millseconds.");
        }
        catch(Exception e) {
            logger.log(Level.SEVERE, "Exception " + e.getMessage(), e);
        }
    }

}
