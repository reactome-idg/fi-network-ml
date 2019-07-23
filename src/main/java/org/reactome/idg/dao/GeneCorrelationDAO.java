package org.reactome.idg.dao;

import java.util.Map;

import org.reactome.idg.model.Gene;
import org.reactome.idg.model.Provenance;

public interface GeneCorrelationDAO
{
	public void setCurrentProvenance(Provenance p);
	
	public Provenance getCurrentProvenance();
	
	public Long addGenePair(Gene gene1, Gene gene2, double correlationValue);
	
	public Map<Provenance, Double> getCorrelation(Gene gene1, Gene gene2);
	
	public Double getCorrelation(Gene gene1, Gene gene2, Provenance prov);
	
	public void setBatchSize(int batchSize);
	
	public int getBatchSize();
	
	public int getNumTxOps();
	
	public void loadGenePairsFromDataFile(String pathToFile);
}
