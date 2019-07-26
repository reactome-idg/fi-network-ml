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
	
	@Autowired
	private SessionFactory sessionFactory;

	private Session session;
	
	/**
	 * @see {@link GeneDAO#addGene(String)}
	 */
	@Transactional(readOnly = false)
	@Override
	public Long addGene(String symbol)
	{
		getSession();
		session.setCacheMode(CacheMode.NORMAL);
		Gene gene = new Gene();
		gene.setSymbol(symbol);
		Long newGeneID = (Long) session.save(gene);
		return newGeneID;
	}

	private void getSession()
	{
		// TODO: Move this to some common class - all the other DAO Impls use this pattern somewhere.
		session = sessionFactory.getCurrentSession();
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
	}

	/**
	 * @see {@link GeneDAO#getGene(String)}
	 */
	@Transactional(readOnly = true)
	@Override
	public List<Gene> getGene(String symbol)
	{
		getSession();
		@SuppressWarnings("unchecked")
		List<Gene> results = session.createQuery("from Gene where symbol = :symbol").setParameter("symbol", symbol)
									.setCacheable(true)
									.setCacheMode(CacheMode.GET)
									.getResultList();
		return results;
	}

	/**
	 * @see {@link GeneDAO#getGene(Long)}
	 */
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

	/**
	 * @see {@link GeneDAO#getSymbolToIdMapping()}
	 */
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
	
	/**
	 * @see {@link GeneDAO#getAllGenes()}
	 */
	@Transactional(readOnly = true)
	@Override
	public List<Gene> getAllGenes()
	{
		getSession();
		@SuppressWarnings("unchecked")
		List<Gene> genes = session.createQuery("from Gene").getResultList();
		return genes;
	}
	
	/**
	 * @see {@link GeneDAO#addGenes(List)}
	 */
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
//			catch (ConstraintViolationException e)
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
