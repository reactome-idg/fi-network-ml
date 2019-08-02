package org.reactome.idg.service;

import java.util.List;
import java.util.Map;

import org.reactome.idg.dao.GeneDAO;
import org.reactome.idg.dao.GenePairCorrelationDAO;
import org.reactome.idg.dao.ProvenanceDAO;
import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Isolation;
import org.springframework.transaction.annotation.Transactional;

@Service
public class GeneCorrelationServiceImpl implements GeneCorrelationService {
    
    @Autowired
    private GenePairCorrelationDAO correlationDAO;
    @Autowired
    private GeneDAO geneDAO;
    @Autowired
    private ProvenanceDAO provenanceDAO;

    public GeneCorrelationServiceImpl() {
    }
    
    @Override
    @Transactional
    public Provenance fetchProvenance(Provenance template) {
        return provenanceDAO.addProvenance(template);
    }
    
    @Override
    @Transactional
    public void updateGene(Gene gene) {
        geneDAO.updateGene(gene);
    }

    @Override
    @Transactional
    public Gene fetchGene(String symbol) {
        Gene gene = geneDAO.getGene(symbol);
        if (gene == null)
            gene = geneDAO.addGene(symbol);
        return gene;
    }

    @Override
    public Map<Provenance, Double> getCorrelation(Gene gene1, Gene gene2) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Double getCorrelation(Gene gene1, Gene gene2, Provenance prov) {
        // TODO Auto-generated method stub
        return null;
    }

    @Transactional
    @Override
    public void saveCorrelation(GenePairCorrelation correlation) {
        correlationDAO.save(correlation);
    }
    
    @Transactional
    @Override
    public void saveCorrelations(List<GenePairCorrelation> correlations) {
        correlationDAO.save(correlations);
    }

    @Override
    @Transactional(isolation = Isolation.READ_UNCOMMITTED)
    public void loadGenePairsFromDataFile(String pathToFile)
    {
//        Session session = sessionFactory.getCurrentSession();
//        // play with some tuning parameters... 
//        //		this.session.createSQLQuery("set global innodb_buffer_pool_size=8053063680;").executeUpdate();
//        //		this.session.createSQLQuery("set global innodb_io_capacity=5000;").executeUpdate();
//        //		this.session.createSQLQuery("set global innodb_io_capacity_max=20000;").executeUpdate();
//        //		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 2;").executeUpdate();
//        session.createSQLQuery("SET global unique_checks=0;").executeUpdate();
//        session.createSQLQuery("SET global autocommit=0;").executeUpdate();
//        //		this.session.createSQLQuery("SET sql_log_bin ='OFF';").executeUpdate();
//        LocalDateTime start = LocalDateTime.now();
//        NativeQuery<?> nq = session.createSQLQuery("LOAD DATA LOCAL INFILE :file INTO TABLE gene_pair_correlation"
//                + " LINES TERMINATED BY '\\n' "
//                + " (gene_1_id, gene_2_id, correlation_value, provenance_id);").setParameter("file", pathToFile);
//        long numRows = nq.executeUpdate();
//        LocalDateTime end = LocalDateTime.now();
//        //		this.session.createSQLQuery("SET sql_log_bin ='ON';").executeUpdate();
//        session.createSQLQuery("SET global unique_checks=1;").executeUpdate();
//        //		this.session.createSQLQuery("SET global innodb_autoinc_lock_mode = 1;").executeUpdate();
//        logger.info("Number of rows loaded: {}, time duration: {}", numRows, Duration.between(start, end).toString());
    }

}
