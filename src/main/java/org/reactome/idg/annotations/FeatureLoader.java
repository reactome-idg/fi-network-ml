package org.reactome.idg.annotations;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Methods annotated with this Annotation are used to load features for ML.
 * @author wug
 *
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface FeatureLoader {

    /**
     * The full name of method, including package and class name. 
     * This is a messy usage: if the method name has two parts separated by
     * ",": the first part is the feature name and the second is the actual
     * method name.
     * @return
     */
    String[] methods();
    
}
