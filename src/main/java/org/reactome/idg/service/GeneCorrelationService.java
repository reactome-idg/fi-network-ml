package org.reactome.idg.service;

import java.util.List;
import java.util.Map;

import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;
import org.reactome.idg.model.Provenance;

/**
 * The entry point to access to the database.
 * @author wug
 */
public interface GeneCorrelationService {

    public Map<Provenance, Double> getCorrelation(Gene gene1, Gene gene2);
	
	public Double getCorrelation(Gene gene1, Gene gene2, Provenance prov);
	
	public void loadGenePairsFromDataFile(String pathToFile);
	
	public void saveCorrelation(GenePairCorrelation correlation);
	
	public void saveCorrelations(List<GenePairCorrelation> correlations);
	
	public void updateGene(Gene gene);
	
	/**
	 * Fetch the gene for the provided symbol. If there is no gene in the database, then
	 * a new Gene object will be created and saved into the database.
	 * @param symbol
	 * @return
	 */
	public Gene fetchGene(String symbol);
	
	/**
	 * Fetch a matched Provenance object from the database. If there is no match,
	 * the template will be persisted into the database.
	 * @param template
	 * @return
	 */
	public Provenance fetchProvenance(Provenance template);
}
