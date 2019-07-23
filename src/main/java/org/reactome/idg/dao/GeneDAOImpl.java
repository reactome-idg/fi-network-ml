package org.reactome.idg.dao;


import java.util.List;

import org.hibernate.CacheMode;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.idg.model.Gene;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.transaction.annotation.Transactional;

public class GeneDAOImpl implements GeneDAO
{
//	private static final Logger logger = LogManager.getLogger();

	// The symbol-to-ID map is a cache intended to speed up the process of creating bulk-load files.
	// The bulk-load files will contain on each row: two gene IDs, and a correlation value.
	// Reading a correlation matrix will require lookups of gene symbols to gene IDs. This cache
	// should save on round-trips to the database.
//	private static Map<String, Long> symbolToIdCache = new HashMap<>();
	
	@Autowired
	private SessionFactory sessionFactory;

	private Session session;
	
	@Transactional(readOnly = false)
	@Override
	public Long addGene(String symbol)
	{
		// If the symbol is already in the local cache, return early!
//		if (symbolToIdCache.containsKey(symbol))
//		{
//			return symbolToIdCache.get(symbol);
//		}
//		else
		{
			// symbol has not already been loaded, so let's add it to the database.
			session = sessionFactory.getCurrentSession();
			if (!session.isOpen())
			{
				session = sessionFactory.openSession();
			}
			session.setCacheMode(CacheMode.NORMAL);
			Gene gene = new Gene();
			gene.setSymbol(symbol);
			
			Long newGeneID = (Long) session.save(gene);
			
//			symbolToIdCache.put(symbol, newGeneID);
			
			return newGeneID;
		}
	}

	@Transactional(readOnly = true)
	@Override
	public List<Gene> getGene(String symbol)
	{
		session = sessionFactory.getCurrentSession();
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
//		if (symbolToIdCache.containsKey(symbol))
//		{
//			Gene g = new Gene();
//			g.setId(symbolToIdCache.get(symbol));
//			g.setSymbol(symbol);
//			return Arrays.asList(g);
//		}
//		else
		{
			@SuppressWarnings("unchecked")
			List<Gene> results = session.createQuery("from Gene where symbol = :symbol").setParameter("symbol", symbol)
										.setCacheable(true)
										.setCacheMode(CacheMode.GET)
										.getResultList();
			
			return results;
		}
	}

	@Transactional(readOnly = true)
	@Override
	public List<Gene> getGene(Long id)
	{
		session = sessionFactory.getCurrentSession();
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		@SuppressWarnings("unchecked")
		List<Gene> results = session.createQuery("from Gene where id = :id").setParameter("id", id)
									.setCacheable(true)
									.setCacheMode(CacheMode.GET)
									.getResultList();
		
		return results;
	}

}
