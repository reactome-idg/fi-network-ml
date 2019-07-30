package org.reactome.idg.dao;

import java.util.List;

import javax.persistence.TypedQuery;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;
import org.springframework.beans.factory.annotation.Autowired;


/**
 * An implementation of GenePairCorrelationDAO based on hibernate API.
 * @author wug
 *
 */
public class GenePairCorrelationDAOImpl implements GenePairCorrelationDAO {
    
    @Autowired
    private SessionFactory sessionFactory;
    
    public GenePairCorrelationDAOImpl() {
    }

    @Override
    public GenePairCorrelation get(Long id) {
        Session session = sessionFactory.getCurrentSession();
        GenePairCorrelation correlation = session.load(GenePairCorrelation.class, id);
        return correlation;
    }

    @Override
    public void save(GenePairCorrelation correlation) {
        Session session = sessionFactory.getCurrentSession();
        session.save(correlation);
    }

    @Override
    public List<GenePairCorrelation> fetch(Gene gene1, Gene gene2) {
        Session session = sessionFactory.getCurrentSession();
        TypedQuery<GenePairCorrelation> query = session.createQuery("FROM GenePairCorrelation as c WHERE c.gene1 = :gene1 and c.gene2 = :gene2",
                                                                    GenePairCorrelation.class);
        List<GenePairCorrelation> results = query.setParameter("gene1", gene1)
                                                 .setParameter("gene2", gene2)
                                                 .getResultList();
        return results;
    }
    

}
