package org.reactome.idg.dao;

import java.util.List;
import java.util.Map;

import org.reactome.idg.model.Gene;

public interface GeneDAO
{
	/**
	 * Adds a gene symbol to the database.
	 * @param symbol
	 * @return
	 */
	public Integer addGene(String symbol);
	
	/**
	 * Adds a list of gene symbols to the database. Use this over {@link GeneDAO#addGene(String)} when you have a large number of genes to add,
	 * as this will use batch inserts, and should run much faster than adding genes one by one.
	 * @param symbols
	 */
	public void addGenes(List<String> symbols);
	
	/**
	 * Gets a Gene that match the given symbol.
	 * @param symbol - The symbol to look up.
	 * @return A list of Gene objects. There should only ever be one.
	 */
	public Gene getGene(String symbol);
	
	/**
	 * Gets a Gene that match the given symbol.
	 * @param id - The ID to look up.
	 * @return A list of Gene objects. There should only ever be one.
	 */
	public Gene getGene(Integer id);
	
	/**
	 * Gets a Symbol-to-ID mapping of Genes from the database. 
	 * @return A mapping of Gene symbols mapped to their IDs in the database.
	 */
	Map<String, Integer> getSymbolToIdMapping();
	
	/**
	 * Returns all genes.
	 * @return A List of all Genes.
	 */
	List<Gene> getAllGenes();
}
