package org.reactome.idg.model;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.sql.SQLException;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.dao.GeneDAO;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class TestGeneDAO
{

	
	@Test
	public void testAddGeneIT() throws SQLException
	{
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			
			context.register(AppConfig.class);
			context.refresh();
						
			GeneDAO dao = (GeneDAO) context.getBean("geneDao");
			
			Integer id = dao.addGene("TEST_GENE");
			assertTrue(id.longValue() > 0);
			
			Gene g = dao.getGene(id);
			assertNotNull(g);
			
			Gene g2 = dao.getGene("TEST_GENE");
			assertNotNull(g2);
		}
	}
}
