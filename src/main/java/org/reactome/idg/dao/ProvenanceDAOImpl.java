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
    @Autowired
    private SessionFactory sessionFactory;

    @Override
    public Provenance getProvenanceById(Integer id) {
        Session session = sessionFactory.getCurrentSession();
        Provenance provenance = session.load(Provenance.class, id);
        return provenance;
    }


    @Override
    public List<Provenance> getProvenanceByName(String name)
    {
        Session session = sessionFactory.getCurrentSession();

        List<Provenance> provenances = session.createQuery("FROM Provenance where name = :name", Provenance.class)
                                               .setParameter("name", name)
                                               .getResultList();
        return provenances;
    }

    /**
     * @see {@link org.reactome.idg.dao.ProvenanceDAO#addProvenance(Provenance)}
     */
    @Override
    public Provenance addProvenance(Provenance p) {
        Session session = sessionFactory.getCurrentSession();
        List<Provenance> results = session.createQuery("from Provenance where name = :name and url = :url and category = :cat and subcategory = :subcat", Provenance.class)
                                          .setParameter("name", p.getName())
                                          .setParameter("url", p.getUrl())
                                          .setParameter("cat", p.getCategory())
                                          .setParameter("subcat", p.getSubcategory())
                                          .getResultList();

        if (null == results || results.size() == 0) {
            session.save(p);
            return p; // The transient is now persisted. There is no need to perform another call
        }
        else {
            return results.get(0);
        }
    }

}
