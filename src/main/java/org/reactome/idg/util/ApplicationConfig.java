package org.reactome.idg.util;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Properties;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;

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
        InputStream is = getClass().getClassLoader().getResourceAsStream(APPLICATION_CONFIG);
        appConfig = new Properties();
        try {
            appConfig.load(is);
        }
        catch(IOException e) {
            logger.error(e.getMessage(), e);
        }
    }
    
    public Map<String, String> getGeneToUniProMap() throws IOException {
        String fileName = getConfig().getAppConfig("reactome.uniprot.to.gene");
        URL url = getClass().getClassLoader().getResource(fileName);
        try (Stream<String> lines = Files.lines(Paths.get(url.getFile()))) {
            return lines.skip(1)
                        .map(line -> line.split("\t"))
                        .collect(Collectors.toMap(tokens -> tokens[1], 
                                                  tokens -> tokens[0],
                                                  (protein1, protein2) -> protein1)); // In case two proteins are mapped to the same gene, choose the first protein.
        }
    }
    
    public Map<String, String> getUniProtToGeneMap() throws IOException {
        String fileName = getConfig().getAppConfig("reactome.uniprot.to.gene");
        URL url = getClass().getClassLoader().getResource(fileName);
        try (Stream<String> lines = Files.lines(Paths.get(url.getFile()))) {
            return lines.skip(1)
                        .map(line -> line.split("\t"))
                        .collect(Collectors.toMap(tokens -> tokens[0], 
                                                  tokens -> tokens[1],
                                                  (gene1, gene2) -> gene1)); // In case two genes are mapped to the same protein, choose the first one.
        }
    }
    
}
