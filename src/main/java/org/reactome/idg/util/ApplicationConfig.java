package org.reactome.idg.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to load configurations from the resource folder. This is a singleton class.
 * @author wug
 *
 */
public class ApplicationConfig {
    private Logger logger = Logger.getLogger(ApplicationConfig.class);
    private static ApplicationConfig config;
    private final String APPLICATION_CONFIG = "application.properties";
    private Properties appConfig;
    private FIConfiguration fiConfig;
    
    private ApplicationConfig() {
        load();
    }
    
    public static ApplicationConfig getConfig() {
        if (config == null)
            config = new ApplicationConfig();
        return config;
    }
    
    public String getAppConfig(String propName) {
        return appConfig.getProperty(propName);
    }

    private void load() {
        try {
            appConfig = new Properties();
            InputStream is = getInputStream(APPLICATION_CONFIG);
            if (is == null) {
                logger.warn("Cannot find the configuration file: " + APPLICATION_CONFIG);
            }
            else {
                appConfig.load(is);
                is.close();
            }
        }
        catch(IOException e) {
            logger.error(e.getMessage(), e);
        }
    }
    
    public double getMaximumCutoff() {
        String value = getConfig().getAppConfig("maximum.cutoff");
        if (value == null || value.length() == 0)
            value = "0.999";
        return Double.parseDouble(value);
    }
    
    public Map<String, String> getGeneToUniProMap() throws IOException {
        String fileName = getConfig().getAppConfig("reactome.uniprot.to.gene");
        InputStream is = getInputStream(fileName);
        if (is == null) {
            logger.error("Cannot find file: " + fileName);
            return null;
        }
        try (Stream<String> lines = new BufferedReader(new InputStreamReader(is)).lines()) {
            return lines.skip(1)
                        .map(line -> line.split("\t"))
                        .collect(Collectors.toMap(tokens -> tokens[1], 
                                                  tokens -> tokens[0],
                                                  (protein1, protein2) -> protein1)); // In case two proteins are mapped to the same gene, choose the first protein.
        }
    }
    
    public Set<String> getAllGenes() throws IOException {
        return getGeneToUniProMap().keySet();
    }
    
    public Map<String, String> getUniProtToGeneMap() throws IOException {
        String fileName = getConfig().getAppConfig("reactome.uniprot.to.gene");
        InputStream is = getInputStream(fileName);
        if (is == null) {
            logger.error("Cannot find file: " + fileName);
            return null;
        }
        try (Stream<String> lines = new BufferedReader(new InputStreamReader(is)).lines()) {
            return lines.skip(1)
                        .map(line -> line.split("\t"))
                        .collect(Collectors.toMap(tokens -> tokens[0], 
                                                  tokens -> tokens[1],
                                                  (gene1, gene2) -> gene1)); // In case two genes are mapped to the same protein, choose the first one.
        }
    }
    
    public FIConfiguration getFIConfig() {
        if (fiConfig == null) {
            // Make sure the following code will be called once only.
            // Need to reset configuration
            String fileName = "configuration.prop";
            try {
                InputStream is = getInputStream(fileName);
                fiConfig = FIConfiguration.getConfiguration(is);
            }
            catch(IOException e) {
                logger.error(e.getMessage(), e);
            }
        }
        return fiConfig;
    }
    
    public Set<String> loadReactomeFIsInGenes() throws IOException {
        String fileName = getConfig().getAppConfig("reactome.fi.uniprot.file");
        InputStream is = getInputStream(fileName);
        if (is == null) {
            logger.error("Cannot find file: " + fileName);
            return null;
        }
        BufferedReader br = new BufferedReader(new InputStreamReader(is));
        List<String> fis = br.lines().collect(Collectors.toList());
        br.close();
        logger.info("Total loaded Reactome FIs in UniProt: " + fis.size());
        return mapFIsInUniProtToGenes(fis);
    }
    
    private Set<String> mapFIsInUniProtToGenes(List<String> fis) throws IOException {
        Map<String, String> uni2Gene = getConfig().getUniProtToGeneMap();
        Set<String> fisInGenes = fis.stream()
                .map(line -> line.split("\t"))
                .map(tokens -> {
                    String gene1 = uni2Gene.get(tokens[0]);
                    String gene2 = uni2Gene.get(tokens[1]);
                    if (gene1 == null || gene2 == null || gene1.equals(gene2))
                        return null;
                    return InteractionUtilities.generateFIFromGene(gene1, gene2);
                })
                .filter(fi -> fi != null)
                .distinct()
                .collect(Collectors.toSet());
        return fisInGenes;
    }
    
    public InputStream getInputStream(String fileName) throws IOException {
        File file = new File(fileName);
        if (file.exists()) {
            logger.info("Get " + fileName + " from " + file.getAbsolutePath());
            return new FileInputStream(file);
        }
        // Do another search
        file = new File("resources" + File.separator + fileName);
        if (file.exists()) {
            logger.info("Get " + fileName + " from " + file.getAbsolutePath());
            return new FileInputStream(file);
        }
        // Try to use the resource
        InputStream is = ApplicationConfig.class.getClassLoader().getResourceAsStream(fileName);
        return is; // We have done our best.
    }
    
}
