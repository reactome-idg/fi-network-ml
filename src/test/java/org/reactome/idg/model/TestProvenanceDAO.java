package org.reactome.idg.model;



import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import javax.print.attribute.standard.PrinterLocation;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.dao.ProvenanceDAO;
import org.reactome.idg.model.Provenance;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

@SuppressWarnings("static-method")
public class TestProvenanceDAO
{

	@Test
	public void testAddProvenanceIT()
	{
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			ProvenanceDAO provDao = (ProvenanceDAO) context.getBean("provenanceDao");
			
			Provenance p = new Provenance();
			p.setName("TEST" + (new Random()).nextInt(1000000)); //random number is to ensure that values inserted by previous test run will not affect this test run.
			p.setCategory("TEST Category");
			p.setSubcategory("Test SubCategory");
			p.setUrl("http://testUrl.com");
			Provenance p1 = provDao.addProvenance(p);
			assertNotNull(p1);
			System.out.println(p1.toString());
			Provenance pFromDB = provDao.getProvenanceById(p1.getId());
			assertEquals(pFromDB.getName(), p.getName());
			System.out.println(pFromDB.toString());
			List<Provenance> provList = provDao.getProvenanceByName(p.getName());
			assertTrue(provList.size() > 0);
			System.out.println(provList.size() + " record(s) found.");
			provList.stream().forEach( provenance -> assertEquals("Names don't match! From DB: " + provenance.getName() + " original value: "+p.getName(), provenance.getName(), p.getName()));
			// This should trigger a "Provenance already exists!" message.
			Provenance p2 = provDao.addProvenance(p);
			// trying to add the same provenance a second time should result in the first instance being returned, so the IDs should be the same.
			assertEquals(p1.getId(), p2.getId());
		}
	}
}
