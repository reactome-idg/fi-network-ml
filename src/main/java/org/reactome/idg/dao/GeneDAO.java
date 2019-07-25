package org.reactome.idg.dao;

import java.util.List;
import java.util.Map;

import org.reactome.idg.model.Gene;

public interface GeneDAO
{
	public Long addGene(String symbol);
	public void addGenes(List<String> symbols);
	public List<Gene> getGene(String symbol);
	public List<Gene> getGene(Long id);
	Map<String, Long> getSymbolToIdMapping();
	List<Gene> getAllGenes();
}
