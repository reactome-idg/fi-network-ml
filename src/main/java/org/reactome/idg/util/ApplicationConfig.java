package org.reactome.idg.util;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

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
    
}
