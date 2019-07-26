package org.reactome.harmonizome;

import org.reactome.idg.config.AppConfig;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class ImportHarmonizomeData
{

	/**
	 * Entrypoint to code that imports Harmonizome data. 
	 * This build this application (using the HarmonizomeImporter profile):
	 *  mvn clean package -DskipTests -PHarmonizomeImporter
	 * to run: 
	 *  java -jar target/HarmonizomeImporter-jar-with-dependencies.jar ./src/main/resources/application.properties
	 */
	public static void main(String[] args)
	{
		String pathToProps = "./application.properties";
		// first argument should be path to properties file.
		if (args.length > 0)
		{
			pathToProps = args[0];
		}
		System.setProperty("pathToProperties", pathToProps);
		
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();

			// Get the batch bean. It will be autowired with a list of HarmonizomeDownloaders, which it will execute.
			HarmonizomeBatch batch = (HarmonizomeBatch) context.getBean("harmonizomeBatch");
			// download the files.
			batch.downloadFiles();
			
			try
			{
				// load to the database.
				batch.loadFiles();
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
}
