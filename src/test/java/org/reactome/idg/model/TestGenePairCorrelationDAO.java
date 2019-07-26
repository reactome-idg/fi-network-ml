package org.reactome.idg.model;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.sql.SQLException;
import java.util.Random;

import javax.sql.DataSource;

import org.junit.Test;
import org.reactome.harmonizome.config.AppConfig;
import org.reactome.idg.dao.GeneCorrelationDAO;
import org.reactome.idg.dao.GeneDAO;
import org.reactome.idg.dao.ProvenanceDAO;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

@SuppressWarnings("static-method")
public class TestGenePairCorrelationDAO
{
	@Test
	public void testAddGenePairCorrelationIT() throws SQLException
	{
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			GeneCorrelationDAO geneCorrelationDao = (GeneCorrelationDAO) context.getBean("dao");
			DataSource ds = (DataSource)context.getBean("dataSource");
			Random rand = new Random();
			String geneSymbol1 = "GENE_"+rand.nextInt(1000);
			String geneSymbol2 = "GENE_"+rand.nextInt(1000);
			double correlationValue = 0.1234;
			double longCorrelationValue = 0.123456789123456789; // deliberately longer than allowed precision.
			
			Provenance prov = new Provenance("TestProvenance_"+rand.nextInt(1000), "http://www.test.com", "CAT", "SUBCAT");
			
			ProvenanceDAO provDao = (ProvenanceDAO) context.getBean("provenanceDao");
			
			Provenance p = provDao.addProvenance(prov);
			
			assertNotNull(p);
			
			Gene gene1;
			Gene gene2;
			
			GeneDAO geneDao = (GeneDAO) context.getBean("geneDao");
			Long g1id = geneDao.addGene(geneSymbol1);
			Long g2id = geneDao.addGene(geneSymbol2);

			assertTrue(g1id > 0);
			assertTrue(g2id > 0);
			
			gene1 = geneDao.getGene(g1id).get(0);
			gene2 = geneDao.getGene(g2id).get(0);
			
			assertNotNull(gene1);
			assertNotNull(gene2);
			
//			GenePairCorrelation corr = new GenePairCorrelation(gene1, gene2, longCorrelationValue, p);
			geneCorrelationDao.setCurrentProvenance(p);
			geneCorrelationDao.setBatchSize(1);
			Long corrId = geneCorrelationDao.addGenePair(gene1, gene2, correlationValue);
			assertTrue(corrId > 0);
			
			Double corrFromDb = geneCorrelationDao.getCorrelation(gene1, gene2, p);
			assertTrue(corrFromDb.doubleValue() == correlationValue);
			System.out.println(corrFromDb);
			Double corrFromDb1 = geneCorrelationDao.getCorrelation(gene2, gene1, p);
			assertTrue(corrFromDb1.doubleValue() == correlationValue);
			System.out.println(corrFromDb1);
			// delete the correlation before recreating it with a different value
			ds.getConnection().setAutoCommit(true);
			ds.getConnection().prepareStatement("DELETE FROM gene_pair_correlation WHERE gene_1_id = "+gene1.getId()+" AND gene_2_id = "+gene2.getId()+" AND provenance_id = "+p.getId()).execute();
			ds.getConnection().setAutoCommit(false);
//			ds.getConnection().commit();
			Long corrId2 = geneCorrelationDao.addGenePair(gene1, gene2, longCorrelationValue);
			assertTrue(corrId2 > 0);

			Double corrFromDb2 = geneCorrelationDao.getCorrelation(gene1, gene2, p);
			
			assertNotNull(corrFromDb2);
			System.out.println(corrFromDb2);
			assertTrue(corrFromDb2.doubleValue() != longCorrelationValue);
			ds.getConnection().setAutoCommit(true);
			// Cleanup time!
			ds.getConnection().prepareStatement("DELETE FROM gene_pair_correlation WHERE gene_1_id = "+gene1.getId()+" AND gene_2_id = "+gene2.getId()+" AND provenance_id = "+p.getId()).execute();
			ds.getConnection().prepareStatement("DELETE FROM provenance WHERE id = "+p.getId()).execute();
			ds.getConnection().prepareStatement("DELETE FROM gene WHERE id = "+gene1.getId()).execute();
			ds.getConnection().prepareStatement("DELETE FROM gene WHERE id = "+gene2.getId()).execute();
			ds.getConnection().setAutoCommit(false);
//			ds.getConnection().commit();
		}
	}
}
