package org.reactome.idg.util;

/**
 * This abstract class is used to load the configuration for the FINetworkConfiguration singleton.
 * It is used so that other classes in this package don't need to configure that singleton again
 * as long as they subclass to this class, which requires nothing.
 * Note: We don't really need any configuration from the configuration file for FI network build.
 * However, since we use some of classes there, which require FINetworkConfiguration during construction
 * time. We may need to refactor the code there to avoid this dependency.
 * @author wug
 *
 */
public abstract class FINetworkBuildConfig {

    static {
        ApplicationConfig.getConfig().getFIConfig();
    }

    protected FINetworkBuildConfig() {
    }
    
}