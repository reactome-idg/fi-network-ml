package org.reactome.idg.dao;

import java.util.List;

import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;

/**
 * A minimum interface to handle persistance with GenePairCorrelation objects.
 * @author wug
 *
 */
public interface GenePairCorrelationDAO {
    
    public GenePairCorrelation get(Long id);
    
    public void save(GenePairCorrelation correlation);
    
    public List<GenePairCorrelation> fetch(Gene gene1, Gene gene2);

}
