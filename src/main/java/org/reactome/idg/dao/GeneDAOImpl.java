package org.reactome.idg.dao;


import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.CacheMode;
import org.hibernate.FlushMode;
import org.hibernate.JDBCException;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.exception.ConstraintViolationException;
import org.reactome.idg.model.Gene;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.transaction.annotation.Transactional;

public class GeneDAOImpl implements GeneDAO
{
	private static final Logger logger = LogManager.getLogger();

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
			getSession();
			session.setCacheMode(CacheMode.NORMAL);
			Gene gene = new Gene();
			gene.setSymbol(symbol);
			
			Long newGeneID = (Long) session.save(gene);
			
//			symbolToIdCache.put(symbol, newGeneID);
			
			return newGeneID;
		}
	}

	private void getSession()
	{
		session = sessionFactory.getCurrentSession();
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
	}

	@Transactional(readOnly = true)
	@Override
	public List<Gene> getGene(String symbol)
	{
		getSession();
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
		getSession();
		@SuppressWarnings("unchecked")
		List<Gene> results = session.createQuery("from Gene where id = :id").setParameter("id", id)
									.setCacheable(true)
									.setCacheMode(CacheMode.GET)
									.getResultList();
		
		return results;
	}

	@Transactional(readOnly = true)
	@Override
	public Map<String, Long> getSymbolToIdMapping()
	{
		getSession();
		Map<String, Long> symbolsToIds = new HashMap<>();
		
		@SuppressWarnings("unchecked")
		List<Gene> genes = session.createQuery("from Gene").getResultList();
		for (Gene g : genes)
		{
			symbolsToIds.put(g.getSymbol(), g.getId());
		}
		
		return symbolsToIds;
	}
	
	@Transactional(readOnly = true)
	@Override
	public List<Gene> getAllGenes()
	{
		getSession();
		@SuppressWarnings("unchecked")
		List<Gene> genes = session.createQuery("from Gene").getResultList();
		return genes;
	}
	
	@Override
	public void addGenes(List<String> symbols)
	{
		if (this.session == null || !this.session.isOpen())
		{
			this.session = sessionFactory.openSession();
		}
		this.session.setJdbcBatchSize(500);
//		this.session.setHibernateFlushMode(FlushMode.COMMIT);
		Transaction tran = this.session.beginTransaction();
		int i = 0;
		for (String symbol : symbols)
		{
			i++;
			Gene gene = new Gene();
			gene.setSymbol(symbol);
//			try
//			{
				session.save(gene);
//			}
//			catch (Error e)
//			{
//				// write out the message, but don't crash.
//				if (e.getMessage().contains("Duplicate entry"))
//				{
//					logger.warn(e.getMessage());
//				}
//				else
//				{
//					// any other constraint violation should be rethrown...
//					throw e;
//				}
//			}
			if (i % this.session.getJdbcBatchSize() == 0)
			{
				this.session.flush();
				this.session.clear();
			}
		}
		tran.commit();
		this.session.close();
	}

}
