package org.reactome.idg.annotations;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

import org.reactome.idg.model.FeatureSource;

import net.bytebuddy.dynamic.loading.PackageDefinitionStrategy.Definition.Undefined;

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
    
    /**
     * Use feature source. Default is null.
     */
    FeatureSource source() default FeatureSource.UNDEFINED;
}
