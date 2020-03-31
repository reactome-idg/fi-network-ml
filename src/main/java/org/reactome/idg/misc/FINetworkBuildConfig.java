package org.reactome.idg.misc;

import java.net.URL;

import org.reactome.fi.util.FIConfiguration;

/**
 * This abstract class is used to load the configuration for the FINetworkConfiguration singleton.
 * It is used so that other classes in this package don't need to configure that singleton again
 * as long as they subclass to this class, which requires nothing.
 * @author wug
 *
 */
public abstract class FINetworkBuildConfig {

    static {
        // Need to reset configuration
        String fileName = "configuration.prop";
        URL url = FINetworkBuildConfig.class.getClassLoader().getResource(fileName);
        FIConfiguration.getConfiguration(url.getFile());
    }

    protected FINetworkBuildConfig() {
        
    }
    
}
