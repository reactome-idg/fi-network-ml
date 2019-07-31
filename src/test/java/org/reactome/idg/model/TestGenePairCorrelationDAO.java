package org.reactome.idg.model;

import java.sql.SQLException;

import org.junit.Test;

public class TestGenePairCorrelationDAO
{
	@Test
	public void testAddGenePairCorrelationIT() throws SQLException
	{
//		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
//		{
//			context.register(AppConfig.class);
//			context.refresh();
//			GeneCorrelationService geneCorrelationDao = (GeneCorrelationService) context.getBean("dao");
//			DataSource ds = (DataSource)context.getBean("dataSource");
//			Random rand = new Random();
//			String geneSymbol1 = "GENE_"+rand.nextInt(1000);
//			String geneSymbol2 = "GENE_"+rand.nextInt(1000);
//			double correlationValue = 0.1234;
//			double longCorrelationValue = 0.123456789123456789; // deliberately longer than allowed precision.
//			
//			Provenance prov = new Provenance("TestProvenance_"+rand.nextInt(1000), "http://www.test.com", "CAT", "SUBCAT");
//			
//			ProvenanceDAO provDao = (ProvenanceDAO) context.getBean("provenanceDao");
//			
//			Provenance p = provDao.addProvenance(prov);
//			
//			assertNotNull(p);
//			
//			Gene gene1;
//			Gene gene2;
//			
//			GeneDAO geneDao = (GeneDAO) context.getBean("geneDao");
//			Integer g1id = geneDao.addGene(geneSymbol1);
//			Integer g2id = geneDao.addGene(geneSymbol2);
//
//			assertTrue(g1id > 0);
//			assertTrue(g2id > 0);
//			
//			gene1 = geneDao.getGene(g1id);
//			gene2 = geneDao.getGene(g2id);
//			
//			assertNotNull(gene1);
//			assertNotNull(gene2);
//			
////			GenePairCorrelation corr = new GenePairCorrelation(gene1, gene2, longCorrelationValue, p);
//			geneCorrelationDao.setCurrentProvenance(p);
//			geneCorrelationDao.setBatchSize(1);
//			Long corrId = geneCorrelationDao.addGenePair(gene1, gene2, correlationValue);
//			assertTrue(corrId > 0);
//			
//			Double corrFromDb = geneCorrelationDao.getCorrelation(gene1, gene2, p);
//			assertTrue(corrFromDb.doubleValue() == correlationValue);
//			System.out.println(corrFromDb);
//			Double corrFromDb1 = geneCorrelationDao.getCorrelation(gene2, gene1, p);
//			assertTrue(corrFromDb1.doubleValue() == correlationValue);
//			System.out.println(corrFromDb1);
//			// delete the correlation before recreating it with a different value
//			ds.getConnection().setAutoCommit(true);
//			ds.getConnection().prepareStatement("DELETE FROM gene_pair_correlation WHERE gene_1_id = "+gene1.getId()+" AND gene_2_id = "+gene2.getId()+" AND provenance_id = "+p.getId()).execute();
//			ds.getConnection().setAutoCommit(false);
////			ds.getConnection().commit();
//			Long corrId2 = geneCorrelationDao.addGenePair(gene1, gene2, longCorrelationValue);
//			assertTrue(corrId2 > 0);
//
//			Double corrFromDb2 = geneCorrelationDao.getCorrelation(gene1, gene2, p);
//			
//			assertNotNull(corrFromDb2);
//			System.out.println(corrFromDb2);
//			assertTrue(corrFromDb2.doubleValue() != longCorrelationValue);
//			ds.getConnection().setAutoCommit(true);
//			// Cleanup time!
//			ds.getConnection().prepareStatement("DELETE FROM gene_pair_correlation WHERE gene_1_id = "+gene1.getId()+" AND gene_2_id = "+gene2.getId()+" AND provenance_id = "+p.getId()).execute();
//			ds.getConnection().prepareStatement("DELETE FROM provenance WHERE id = "+p.getId()).execute();
//			ds.getConnection().prepareStatement("DELETE FROM gene WHERE id = "+gene1.getId()).execute();
//			ds.getConnection().prepareStatement("DELETE FROM gene WHERE id = "+gene2.getId()).execute();
//			ds.getConnection().setAutoCommit(false);
//			ds.getConnection().commit();
//		}
	}
}
