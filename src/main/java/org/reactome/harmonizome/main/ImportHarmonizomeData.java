package org.reactome.harmonizome.main;

import java.io.FileInputStream;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.harmonizome.CorrelationMatrixLoader;
import org.reactome.harmonizome.HarmonizomeBatch;
import org.reactome.harmonizome.HarmonizomeDataDownloader;
import org.reactome.harmonizome.config.AppConfig;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class ImportHarmonizomeData
{
	private static final Logger logger = LogManager.getLogger();
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
		Properties props = new Properties();
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
//			try
			{
				props.load(new FileInputStream(pathToProps));
	
				List<HarmonizomeDataDownloader> downloaders = new ArrayList<>();
				Path path = Paths.get(props.getProperty("pathToHarmonizomeFile"));
				List<String> lines = Files.readAllLines(path);
				String downloadPath = props.getProperty("harmonizomeDownloadPath");
				for (String line : lines)
				{
					String[] parts = line.split("\\t");
					HarmonizomeDataDownloader downloader = new HarmonizomeDataDownloader(new URI(parts[2]), parts[0], parts[1], downloadPath);
					downloaders.add(downloader);
				}
				
	
				// Get the batch bean. It will be autowired with a list of HarmonizomeDownloaders, which it will execute.
	//			HarmonizomeBatch batch = (HarmonizomeBatch) context.getBean("harmonizomeBatch");
				HarmonizomeBatch batch = new HarmonizomeBatch();
				batch.setHarmonizomeDownloaders(downloaders);
				context.register(AppConfig.class);
				context.refresh();

				CorrelationMatrixLoader correlationMatrixLoader = (CorrelationMatrixLoader) context.getBean("correlationMatrixLoader");
				batch.setCorrelationMatrixLoader(correlationMatrixLoader );
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
		catch (IOException e)
		{
			logger.error("Problem reading file!");
			e.printStackTrace();
		}
		catch (URISyntaxException e)
		{
			logger.error("Bad URI for Harmonizome file!");
			e.printStackTrace();
		}
	}
}
