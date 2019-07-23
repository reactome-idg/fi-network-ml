package org.reactome.idg.dao;

import java.util.List;

import org.reactome.idg.model.Gene;

public interface GeneDAO
{
	public Long addGene(String symbol);
	public List<Gene> getGene(String symbol);
	public List<Gene> getGene(Long id);
}
