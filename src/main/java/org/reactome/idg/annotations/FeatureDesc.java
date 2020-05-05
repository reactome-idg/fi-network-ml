package org.reactome.idg.annotations;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

import org.reactome.idg.model.FeatureSource;
import org.reactome.idg.model.FeatureSpecies;
import org.reactome.idg.model.FeatureType;

/**
 * Describe the feature (e.g. type, data source, etc)
 * @author wug
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface FeatureDesc {

    FeatureSource[] sources();
    FeatureType type();
    FeatureSpecies species() default FeatureSpecies.Homo_sapiens;
    
}
