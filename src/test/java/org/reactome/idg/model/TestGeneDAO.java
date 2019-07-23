package org.reactome.idg.model;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.sql.SQLException;
import java.sql.Savepoint;
import java.util.List;

import javax.sql.DataSource;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.dao.GeneDAO;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

@SuppressWarnings("static-method")
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
			
			Long id = dao.addGene("TEST_GENE");
			assertTrue(id.longValue() > 0);
			
			List<Gene> g = dao.getGene(id);
			assertNotNull(g);
			assertTrue(g.size() > 0);
			assertTrue(g.get(0).getSymbol().equals("TEST_GENE"));
			
			List<Gene> g2 = dao.getGene("TEST_GENE");
			assertNotNull(g2);
			assertTrue(g2.size() > 0);
			assertTrue(g2.get(0).getSymbol().equals("TEST_GENE"));
			// restore the database to its previous state.
			DataSource ds = (DataSource)context.getBean("dataSource");
			ds.getConnection().prepareStatement("DELETE FROM gene WHERE ID = "+g2.get(0).getId()).execute();
		}
	}
}
