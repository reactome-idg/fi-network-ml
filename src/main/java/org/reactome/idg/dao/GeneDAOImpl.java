package org.reactome.idg.dao;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.persistence.TypedQuery;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.idg.model.Gene;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;

@Repository
public class GeneDAOImpl implements GeneDAO {
	
	@Autowired
	private SessionFactory sessionFactory;

	public GeneDAOImpl() {
    }
	
	/**
	 * @see {@link GeneDAO#addGene(String)}
	 */
	@Override
	public Integer addGene(String symbol)
	{
		Session session = sessionFactory.getCurrentSession();
		Gene gene = new Gene();
		gene.setSymbol(symbol);
		Integer newGeneID = (Integer) session.save(gene);
		return newGeneID;
	}

	/**
	 * @see {@link GeneDAO#getGene(String)}
	 */
	@Override
	public Gene getGene(String symbol) {
	    Session session = sessionFactory.getCurrentSession();
	    TypedQuery<Gene> query = session.createQuery("FROM Gene where symbol = :symbol",
	                                                Gene.class);
	    List<Gene> results = query.setParameter("symbol", symbol).getResultList();
	    return results.isEmpty() ? null : results.get(0);
	}

	/**
	 * @see {@link GeneDAO#getGene(Long)}
	 */
	@Override
	public Gene getGene(Integer id)
	{
	    Session session = sessionFactory.getCurrentSession();
	    Gene gene = session.load(Gene.class, id);
	    return gene;
	}

	/**
	 * @see {@link GeneDAO#getSymbolToIdMapping()}
	 */
	@Override
	public Map<String, Integer> getSymbolToIdMapping()
	{
	    Session session = sessionFactory.getCurrentSession();
		Map<String, Integer> symbolsToIds = new HashMap<>();
		List<Gene> genes = session.createQuery("FROM " + Gene.class.getName(), Gene.class).getResultList();
		for (Gene g : genes)
		{
			symbolsToIds.put(g.getSymbol(), g.getId());
		}
		return symbolsToIds;
	}
	
	/**
	 * @see {@link GeneDAO#getAllGenes()}
	 */
	@Override
	public List<Gene> getAllGenes()
	{
	    Session session = sessionFactory.getCurrentSession();
		List<Gene> genes = session.createQuery("FROM Gene", Gene.class).getResultList();
		return genes;
	}
	
	/**
	 * @see {@link GeneDAO#addGenes(List)}
	 */
	@Override
	public void addGenes(List<String> symbols)
	{
	    Session session = sessionFactory.getCurrentSession();
		for (String symbol : symbols) {
			Gene gene = new Gene();
			gene.setSymbol(symbol);
			session.save(gene);
		}
	}

}
