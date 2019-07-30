package org.reactome.idg.service;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.query.NativeQuery;
import org.reactome.idg.model.Gene;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Isolation;
import org.springframework.transaction.annotation.Transactional;

@Service
public class GeneCorrelationServiceImpl implements GeneCorrelationService
{
    private static final Logger logger = LogManager.getLogger();

    private Provenance currentProvenance;
    // These should be configurable...
    private int numTxOps = 0;
    private int batchSize = 100;

    @Autowired
    private SessionFactory sessionFactory;

    public GeneCorrelationServiceImpl()
    {
    }

    @Override
    @Transactional(isolation = Isolation.READ_UNCOMMITTED)
    public void loadGenePairsFromDataFile(String pathToFile)
    {
        Session session = sessionFactory.getCurrentSession();
        // play with some tuning parameters... 
        //		this.session.createSQLQuery("set global innodb_buffer_pool_size=8053063680;").executeUpdate();
        //		this.session.createSQLQuery("set global innodb_io_capacity=5000;").executeUpdate();
        //		this.session.createSQLQuery("set global innodb_io_capacity_max=20000;").executeUpdate();
        //		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 2;").executeUpdate();
        session.createSQLQuery("SET global unique_checks=0;").executeUpdate();
        session.createSQLQuery("SET global autocommit=0;").executeUpdate();
        //		this.session.createSQLQuery("SET sql_log_bin ='OFF';").executeUpdate();
        LocalDateTime start = LocalDateTime.now();
        NativeQuery<?> nq = session.createSQLQuery("LOAD DATA LOCAL INFILE :file INTO TABLE gene_pair_correlation"
                + " LINES TERMINATED BY '\\n' "
                + " (gene_1_id, gene_2_id, correlation_value, provenance_id);").setParameter("file", pathToFile);
        long numRows = nq.executeUpdate();
        LocalDateTime end = LocalDateTime.now();
        //		this.session.createSQLQuery("SET sql_log_bin ='ON';").executeUpdate();
        session.createSQLQuery("SET global unique_checks=1;").executeUpdate();
        //		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 1;").executeUpdate();
        logger.info("Number of rows loaded: {}, time duration: {}", numRows, Duration.between(start, end).toString());
    }

    /**
     * Adds a gene-pair correlation value to the database. NOTE: The provenance for this gene-pair must be set separately, using <code>setCurrentProvenance()</code>.
     * If you try to insert the same gene-pair twice for the same provenance, you will get back the ID of pre-existing correlation.
     * @param gene1 - the first gene of the pair.
     * @param gene2 - the second gene of the pair.
     * @param correlationValue - the Correlation value.
     */
    @Override
    public Long addGenePair(Gene gene1, Gene gene2, double correlationValue)
    {
        Long corrId = null;
        //		// NOTE: This method is NOT @Transactional because I want to have control over when commits happen.
        //		// Executing a commit EVERY time a record is added is very slow. Executing a commit after some large number
        //		// of INSERTs speeds things up considerably.
        //			
        //		if (session == null || !this.session.isOpen())
        //		{
        //			session = sessionFactory.openSession();
        //		}
        //		
        //		this.session.setHibernateFlushMode(FlushMode.COMMIT);
        //
        //		Gene[] names = {gene1, gene2};
        //		Arrays.sort(names);
        //		
        //		GenePairCorrelation correlation = new GenePairCorrelation(names[0], names[1], correlationValue, this.currentProvenance);
        //		try
        //		{
        //			if (this.numTxOps == 0)
        //			{
        //				this.addGeneTx = session.beginTransaction();
        //			}
        //			corrId = (Long) session.save(correlation);
        //			this.numTxOps++;
        //
        //			// commit when we've executed enough operations.
        //			if (this.numTxOps > this.batchSize - 1)
        //			{
        //				this.numTxOps = 0;
        //				addGeneTx.commit();
        //			}
        //		}
        //		catch (ConstraintViolationException e)
        //		{
        //			if (this.addGeneTx.isActive())
        //			{
        //				addGeneTx.rollback();
        //			}
        //			if (e.getCause().getMessage().contains("Duplicate entry") && e.getCause().getMessage().contains("idx_gene_pair_provenance"))
        //			{
        //				logger.info("It looks like the gene-pair {},{} has already been loaded for that provenance (ID: {}).", gene1, gene2, this.currentProvenance.getId());
        //				correlation = (GenePairCorrelation) session.createQuery("from GenePairCorrelation where gene1 = :g1 AND gene2 = :g2 AND provenance = :p")
        //						.setParameter("g1", names[0])
        //						.setParameter("g2", names[1])
        //						.setParameter("p", this.getCurrentProvenance())
        //						.getSingleResult();
        //				corrId = correlation.getId();
        //			}
        //		}
        //		catch (Exception e)
        //		{
        //			e.printStackTrace();
        //			throw e;
        //		}
        return corrId;
    }

    /**
     * Get all correlations for a pair of gene names.
     */
    @Override
    @Transactional(readOnly = true)
    public Map<Provenance, Double> getCorrelation(Gene gene1, Gene gene2)
    {
        Map<Provenance, Double> provenanceCorrelationMap = new HashMap<>();

        //		getSession();
        //
        //		// We don't store both combinations of a gene pair in the database. A gene-pair is stored once, with genes sorted alphabetically.
        //		// Before we query, we must sort the gene names...
        //		Gene[] names = {gene1, gene2};
        //		Arrays.sort(names);
        //		
        //		@SuppressWarnings("unchecked")
        //		List<GenePairCorrelation> correlations = session.createQuery("from GenePairCorrelation where gene1 = :g1 AND gene2 = :g2 ")
        //														.setParameter("g1", names[0])
        //														.setParameter("g2", names[1])
        //														.getResultList();
        //		for (GenePairCorrelation correlation : correlations)
        //		{
        //			provenanceCorrelationMap.put(correlation.getProvenance(), correlation.getCorrelationValue());
        //		}
        return provenanceCorrelationMap;
    }

    /**
     * Get the correlation value for a pair of genes for a given Provenance.
     */
    @Override
    @Transactional(readOnly = true)
    public Double getCorrelation(Gene gene1, Gene gene2, Provenance prov)
    {
        //		getSession();
        //		
        //		// We don't store both combinations of a gene pair in the database. A gene-pair is stored once, with genes sorted alphabetically.
        //		// Before we query, we must sort the gene names...
        //		Gene[] names = {gene1, gene2};
        //		Arrays.sort(names);
        //		
        //		GenePairCorrelation correlation = (GenePairCorrelation) session.createQuery("from GenePairCorrelation where gene1 = :g1 AND gene2 = :g2 AND provenance = :p")
        //																		.setParameter("g1", names[0])
        //																		.setParameter("g2", names[1])
        //																		.setParameter("p", prov)
        //																		.getSingleResult();
        //		
        return null;
    }

    @Override
    public Provenance getCurrentProvenance()
    {
        return currentProvenance;
    }

    /**
     * Because many gene-pairs can be inserted for a single provenance, it makes sense that "current provenance"
     * is a part of the state of this DAO. Adding a gene pair will be done in the context of the current provenance.
     * Use this method to set the current Provenance. 
     * @param currentProvenance - the Provenance that all gene-pairs will be used when inserting gene-pairs.
     */
    @Override
    public void setCurrentProvenance(Provenance currentProvenance)
    {
        this.currentProvenance = currentProvenance;
    }

    /**
     * Batch size is the number of inserts before a COMMIT is executed.
     * @return the current batch size for this DAO.
     */
    public int getBatchSize()
    {
        return this.batchSize;
    }

    /**
     * Batch size is the number of inserts before a COMMIT is executed. This function sets the batch size.
     * @param batchSize - How many gene-pair correlation insertions should happen before issuing a COMMIT. 
     */
    public void setBatchSize(int batchSize)
    {
        this.batchSize = batchSize;
    }

    /**
     * How many opertaions have occurred since the last COMMIT.
     * @return The number of operations in the current gene-pair correlation insertion transaction. 
     */
    public int getNumTxOps()
    {
        return this.numTxOps;
    }
}
