package org.reactome.idg.dao;

import java.io.File;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;

import javax.persistence.TypedQuery;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.jdbc.Work;
import org.reactome.idg.model.Gene;
import org.reactome.idg.model.GenePairCorrelation;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;


/**
 * An implementation of GenePairCorrelationDAO based on hibernate API.
 * @author wug
 *
 */
@Repository
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
        session.persist(correlation);
    }
    
    @Override
    public void save(List<GenePairCorrelation> correlations) {
        Session session = sessionFactory.getCurrentSession();
        Work work = new Work() {
            @Override
            public void execute(Connection connection) throws SQLException {
//                loadCorrelations(correlations, connection);
//                loadCorrelationsWithInsert(correlations, connection);
                loadCorrelationsWithFile(correlations, connection);
            }
        };
        session.doWork(work);
//        for (GenePairCorrelation corr : correlations)
//            session.persist(corr);
    }
    
    private void loadCorrelations(List<GenePairCorrelation> correlations, Connection connection) throws SQLException {
        PreparedStatement statement = connection.prepareStatement("insert into gene_pair_correlation (correlation_value, gene_1_id, gene_2_id, provenance_id) values (?, ?, ?, ?)");
        for (GenePairCorrelation corr : correlations) {
            statement.setDouble(1, corr.getCorrelationValue());
            statement.setInt(2, corr.getGene1().getId());
            statement.setInt(3, corr.getGene2().getId());
            statement.setInt(4, corr.getProvenance().getId());
            statement.addBatch();
        }
        statement.clearParameters();
        int[] results = statement.executeBatch();
        System.out.println("Size: " + results.length);
    }
    
    private void loadCorrelationsWithFile(List<GenePairCorrelation> correlations, Connection connection) throws SQLException {
        String fileName = "tmp.csv";
        try {
            File file = new File(fileName);
            PrintWriter writer = new PrintWriter(file);
            // Header
            StringBuilder builder = new StringBuilder();
            for (GenePairCorrelation corr : correlations) {
                builder.append(","); // Make sure the first is empty
                builder.append(corr.getCorrelationValue()).append(",");
                builder.append(corr.getGene1().getId()).append(",");
                builder.append(corr.getGene2().getId()).append(",");
                builder.append(corr.getProvenance().getId());
                writer.println(builder.toString());
                builder.setLength(0);
            }
            writer.close();
            Statement st = connection.createStatement();
            st.execute("SET foreign_key_checks=0");
            // Have to use "LOCAL". Otherwise, a security related error is thrown from mysql.
            st.execute("LOAD DATA LOCAL INFILE '" + fileName + "' INTO TABLE gene_pair_correlation FIELDS TERMINATED BY ','");
            st.execute("SET foreign_key_checks=1");
            st.close();
//            file.delete();
            System.out.println("Size: " + correlations.size());
        }
        catch(SQLException e) {
            throw e;
        }
        catch(Exception e) {
            throw new SQLException("Error in creating a temp file: " + e.getMessage()); // Wrap file related error as an SQLException for the method requirement
        }
    }
    
    private void loadCorrelationsWithInsert(List<GenePairCorrelation> correlations, Connection connection) throws SQLException {
        Statement st = connection.createStatement();
        // For parameters test, see: https://dev.mysql.com/doc/refman/5.7/en/optimizing-innodb-bulk-data-loading.html
//        st.execute("SET autocommit=0"); // Don't help improve performance
        st.execute("SET foreign_key_checks=0");
//        st.execute("SET innodb_autoinc_lock_mode=2"); // Read-only. Cannot set!
        StringBuilder builder = new StringBuilder();
        builder.append("INSERT INTO gene_pair_correlation (correlation_value, gene_1_id, gene_2_id, provenance_id) values ");
        for (GenePairCorrelation corr : correlations) {
            builder.append("(");
            builder.append(corr.getCorrelationValue()).append(",");
            builder.append(corr.getGene1().getId()).append(",");
            builder.append(corr.getGene2().getId()).append(",");
            builder.append(corr.getProvenance().getId()).append(")");
            builder.append(",");
        }
        builder.deleteCharAt(builder.length() - 1);
        st.execute(builder.toString());
//        st.execute("COMMIT");
        st.execute("SET foreign_key_checks=1");
//        st.execute("SET innodb_autoinc_lock_mode=1");
        st.close();
        System.out.println("Size: " + correlations.size());
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
