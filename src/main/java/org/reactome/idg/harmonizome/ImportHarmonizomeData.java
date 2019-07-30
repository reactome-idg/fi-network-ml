package org.reactome.idg.harmonizome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.service.GeneCorrelationService;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class ImportHarmonizomeData {
	
    private static final Logger logger = LogManager.getLogger(ImportHarmonizomeData.class);
    
	/**
	 * Entrypoint to code that imports Harmonizome data. 
	 * This build this application (using the HarmonizomeImporter profile):
	 *  mvn clean package -DskipTests -PHarmonizomeImporter
	 * to run: 
	 *  java -jar target/HarmonizomeImporter-jar-with-dependencies.jar ./src/main/resources/application.properties
	 */
	public static void main(String[] args) {
	    AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext(AppConfig.class);
	    GeneCorrelationService service = context.getBean(GeneCorrelationService.class);
	    context.close();
	    logger.info("Finished test!");
//		// first argument should be path to properties file.
//		if (args.length > 0) {
//			pathToProps = args[0];
//		}
//		System.setProperty("pathToProperties", pathToProps);
//		Properties props = new Properties();
//		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();
//			FileInputStream fis = new FileInputStream(pathToProps);)
//		{
//			// Load the properties file.
//			props.load(fis);
//			List<HarmonizomeDataDownloader> downloaders = new ArrayList<>();
//			Path path = Paths.get(props.getProperty("pathToHarmonizomeFile"));
//			List<String> lines = Files.readAllLines(path);
//			String downloadPath = props.getProperty("harmonizomeDownloadPath");
//			
//			// Populate the list of Downloaders, based on the contents of the downloaders file, the location of which
//			// is specified in the properties file.
//			for (String line : lines)
//			{
//				String[] parts = line.split("\\t");
//				HarmonizomeDataDownloader downloader = new HarmonizomeDataDownloader(new URI(parts[2]), parts[0], parts[1], downloadPath);
//				downloaders.add(downloader);
//			}
//			// Create a new batch that will execute the downloaders.
//			HarmonizomeBatch batch = new HarmonizomeBatch();
//			batch.setHarmonizomeDownloaders(downloaders);
//			
//			// Get the matrix loader. We need to get this from the spring context because it depends
//			// on DAO components that are autowired by Spring.
//			context.register(AppConfig.class);
//			context.refresh();
//			CorrelationMatrixLoader correlationMatrixLoader = (CorrelationMatrixLoader) context.getBean("correlationMatrixLoader");
//			batch.setCorrelationMatrixLoader(correlationMatrixLoader );
//			// download the files.
//			batch.downloadFiles();
//			
//			try
//			{
//				// load to the database.
//				batch.loadFiles();
//			}
//			catch (Exception e)
//			{
//				e.printStackTrace();
//			}
//		}
//		catch (IOException e)
//		{
//			logger.error("Problem reading file!");
//			e.printStackTrace();
//		}
//		catch (URISyntaxException e)
//		{
//			logger.error("Bad URI for Harmonizome file!");
//			e.printStackTrace();
//		}
	}
}
