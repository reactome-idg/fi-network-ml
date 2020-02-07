package org.reactome.harmonizome;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Set;

import org.apache.http.client.ClientProtocolException;
import org.junit.Test;
import org.reactome.idg.harmonizome.DataProcessor;

public class TestHarmonizomeDownloader {
    
    @Test
    public void testLoadGeneList() throws Exception {
        DataProcessor processor = new DataProcessor();
        Set<String> genes = processor.getAllGenes();
        System.out.println("Total genes: " + genes.size());
    }

	@Test
	public void testDownloadHarmonizomeFileIT() throws URISyntaxException, ClientProtocolException, IOException
	{
//		URI uri = new URI("https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/kea/gene_similarity_matrix_cosine.txt.gz");
//		String downloadPath = "/media/sshorser/data/reactome/IDGFiles";
//		HarmonizomeDataDownloader downloader = new HarmonizomeDataDownloader(uri, "KEA_Substrates_of_Kinases", "", downloadPath);
//		downloader.downloadFile();
//		String[] parts = uri.getPath().split("/");
//		String fileName = parts[parts.length-1];
//		Path path = Paths.get(downloadPath+"/KEA_Substrates_of_Kinases"+fileName);
//		// make sure file exists.
//		assertTrue(Files.exists(path));
//		// make sure file has some content
//		assertTrue(Files.size(path) > 0);
	}
	
	@Test
	public void testListOfDownloaders()
	{
//		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
//		{
//			context.register(AppConfig.class);
//			context.refresh();
//			@SuppressWarnings("unchecked")
//			List<HarmonizomeDataDownloader> downloaders = (List<HarmonizomeDataDownloader>) context.getBean("harmonizomeDownloaders");
//			assertNotNull(downloaders);
//			assertTrue(downloaders.size() > 0);
//			for (HarmonizomeDataDownloader downloader : downloaders)
//			{
//				System.out.println(downloader.getDatasetName() + "\t" + downloader.getDataCategory() + "\t" + downloader.getUrlToFile().toString());
//			}
//		}
	}
	
	@Test
	public void testBatchIT()
	{
//		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
//		{
//			String pathToProps = "./src/test/resources/application.properties";
//
//			System.setProperty("pathToProperties", pathToProps);
//			context.register(AppConfig.class);
//			context.refresh();
//			
//			HarmonizomeBatch batch = (HarmonizomeBatch) context.getBean("harmonizomeBatch");
//			
//			assertNotNull(batch);
//			batch.downloadFiles();
//			
//			batch.loadFiles();
//		}
//		catch (Exception e)
//		{
//			e.printStackTrace();
//		}
	}
}
