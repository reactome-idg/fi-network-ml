package org.reactome.idg.dao;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Transactional;

@Repository
public class ProvenanceDAOImpl implements ProvenanceDAO
{

	private static final Logger logger = LogManager.getLogger();
	
	@Autowired
	private SessionFactory sessionFactory;

	private Session session;

	
	@Override
	@Transactional(readOnly = true)
	public Provenance getProvenanceById(Long id)
	{
		session = sessionFactory.getCurrentSession();
		
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		
		Provenance provenance = (Provenance) session.createQuery("from Provenance where id = :id")
													.setParameter("id", id)
													.getSingleResult();
		return provenance;
	}

	@Override
	@Transactional(readOnly = true)
	public List<Provenance> getProvenanceByName(String name)
	{
		session = sessionFactory.getCurrentSession();
		
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		
		@SuppressWarnings("unchecked")
		List<Provenance> provenances = (List<Provenance>) session.createQuery("from Provenance where name = :name")
																.setParameter("name", name)
																.getResultList();
		return provenances;
	}

	/**
	 * @see {@link org.reactome.idg.dao.ProvenanceDAO#addProvenance(Provenance)}
	 */
	@Override
	@Transactional(readOnly = false)
	public Provenance addProvenance(Provenance p)
	{
		Provenance createdProvenance;
		
		session = sessionFactory.getCurrentSession();
		
		if (!session.isOpen())
		{
			session = sessionFactory.openSession();
		}
		
		// Before we try to persis this, let's make sure that it's not already there.
		@SuppressWarnings("unchecked")
		List<Provenance> results = session.createQuery("from Provenance where name = :name and url = :url and category = :cat and subcategory = :subcat")
											.setParameter("name", p.getName())
											.setParameter("url", p.getUrl())
											.setParameter("cat", p.getCategory())
											.setParameter("subcat", p.getSubcategory())
											.getResultList();
		

		if (null == results || results.size() == 0)
		{
			Long createdProvenanceId;
			createdProvenanceId = (Long) session.save(p);
			createdProvenance = (Provenance) session.createQuery("from Provenance where id = :id")
													.setParameter("id", createdProvenanceId)
													.getResultList().get(0);
		}
		else
		{
			logger.info("Provenance ({}) already exists, and will not be recreated.", results.get(0).toString());
			createdProvenance = results.get(0);
		}
		
		return createdProvenance;
	}

}
