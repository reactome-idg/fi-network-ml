package org.reactome.harmonizome;



import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;

import org.junit.Test;
import org.reactome.idg.loader.HarmonizomeLoader;

//Some JUnit test methods look like they could be made static, but that causes problems for JUnit.
@SuppressWarnings("static-method")
public class TestHarmonizomeLoader
{

	@Test
	public void testLoadData()
	{
		HarmonizomeLoader loader = new HarmonizomeLoader("src/test/resources/test_gene_association_file.tsv"/*, new Provenance("test", "test", "test", "test")*/);
		try
		{
			long numPairs = loader.loadGeneAssociationData();
			assertTrue(numPairs>0);
		}
		catch (IOException e)
		{
			e.printStackTrace();
			fail("An exception should NOT have been caught in this test.");
		}
	}
	
	@Test
	public void testLoadDataIT()
	{
		// test with a REAL file. Do NOT check that file in to git (it's kinda big). Just download it and place it in src/test/resources 
		HarmonizomeLoader loader = new HarmonizomeLoader("src/test/resources/gene_similarity_matrix_cosine_PathwaysInteractionDatabase.txt"/*, new Provenance("test", "test", "test", "test")*/);
		try
		{
			long numPairs = loader.loadGeneAssociationData();
			// Don't care about specific number, just want to know it is > 0.
			assertTrue(numPairs > 0);
		}
		catch (IOException e)
		{
			e.printStackTrace();
			fail("An exception should NOT have been caught in this test.");
		}
	}
}
