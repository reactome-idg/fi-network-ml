package org.reactome.harmonizome;



import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;
import org.mockito.MockitoAnnotations;
import org.reactome.idg.loader.HarmonizomeLoader;
import org.reactome.idg.model.Provenance;

public class TestHarmonizomeDataRepository
{

	private HarmonizomeLoader loader = new HarmonizomeLoader("src/test/resources/test_gene_association_file.tsv");
	
	Provenance prov = new Provenance("Harmonizome test", null, "test", "test-test");
	
	private Map<Provenance, HarmonizomeLoader> dataLoaders = new HashMap<>(1);
	
	private HarmonizomeDataRepository repository;
	
	@Before
	public void setUp()
	{
		MockitoAnnotations.initMocks(this);
		dataLoaders.put(this.prov, this.loader);
		repository = HarmonizomeDataRepository.createDataRepository(this.dataLoaders);
	}
	
	@Test
	public void testGeneLookup()
	{
		repository.executeDataLoaders();
		
		Map<Provenance, List<Double>> correlations1 = repository.getGeneCorrelation("PRKAR1A", "AMICA1");
		// Should only be 1 thing.
		assertTrue(correlations1.get(prov).size() == 1);
		System.out.println(prov);
		System.out.println(correlations1.get(prov));
		Map<Provenance, List<Double>> correlations2 = repository.getGeneCorrelation("AMICA1", "PRKAR1A");
		// Should get the same values back.
		assertEquals(correlations1, correlations2);
	}
}
