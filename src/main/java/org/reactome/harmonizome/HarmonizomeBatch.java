package org.reactome.harmonizome;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.zip.GZIPInputStream;
import java.util.zip.InflaterInputStream;

import org.apache.http.client.ClientProtocolException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * Processes Harmonizome data in batches. This class will execute set of HarmonizomeDataDownloaders, and then 
 * load the data from each download into the database.
 * 
 * @author sshorser
 *
 */
@Component("harmonizomeBatch")
public class HarmonizomeBatch
{
	private static final Logger logger = LogManager.getLogger();
	
	@Autowired
	private List<HarmonizomeDataDownloader> harmonizomeDownloaders;

	@Autowired
	private CorrelationMatrixLoader correlationMatrixLoader;
	
	private List<String> harmonizomeFilePaths = new ArrayList<>();
	
	public void loadFiles() throws Exception
	{
		for (String filePath : harmonizomeFilePaths)
		{
			correlationMatrixLoader.setFilePath(filePath);
			correlationMatrixLoader.setDelimeter("\\t");
			correlationMatrixLoader.setHeaderSize(3);
			Provenance provenance = new Provenance(filePath, "", "", "");
			correlationMatrixLoader.setProvenance(provenance);
			try
			{
				correlationMatrixLoader.loadData();
			}
			catch (FileNotFoundException e)
			{
				e.printStackTrace();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
	}
	
	public void downloadFiles()
	{
		for (HarmonizomeDataDownloader downloader : this.harmonizomeDownloaders)
		{
			try
			{
				String filePath;
				filePath = downloader.downloadFile();
				String extension = filePath.substring(filePath.lastIndexOf(".")).toLowerCase();
				String outFileName = filePath.replace(extension, "");
				// This code only considers GZip files because that is what is provided by Harmonizome's download page.
				try (InflaterInputStream inflaterStream =  new GZIPInputStream(new FileInputStream(filePath));)
				{
					// This Consumer will write a zip stream to a file.
					// I copied this from a part of AddLinks that could handle zip and gzip files.
					// In this context, we are only expecting gzip so having a data-writer object like this
					// looks like overkill. TODO: Rewrite to be simpler, if possible.
					BiConsumer<InflaterInputStream, String> dataWriter = (inStream, outputFileName) -> {
						byte[] buffer = new byte[1024];
						try
						{
							int count = 0;
							try (FileOutputStream outStream = new FileOutputStream(new File(outputFileName));)
							{
								while ((count = inStream.read(buffer)) > 0)
								{
									outStream.write(buffer, 0, count);
								}
							}
						}
						catch (Exception e)
						{
							logger.error("Error writing zip file: {}", e.getMessage());
						}
					};
					
					dataWriter.accept(inflaterStream, outFileName);
					this.harmonizomeFilePaths.add(outFileName);
				}
				catch (IOException e)
				{
					logger.error("Exception was caught while trying to unzip a file: {}", e.getMessage());
					e.printStackTrace();
				}
			}
			catch (ClientProtocolException e1)
			{
				e1.printStackTrace();
			}
			catch (IOException e1)
			{
				e1.printStackTrace();
			}
		}
	}

	public List<HarmonizomeDataDownloader> getHarmonizomeDownloaders()
	{
		return harmonizomeDownloaders;
	}

	public void setHarmonizomeDownloaders(List<HarmonizomeDataDownloader> harmonizomeDownloaders)
	{
		this.harmonizomeDownloaders = harmonizomeDownloaders;
	}

	public CorrelationMatrixLoader getCorrelationMatrixLoader()
	{
		return correlationMatrixLoader;
	}

	public void setCorrelationMatrixLoader(CorrelationMatrixLoader correlationMatrixLoader)
	{
		this.correlationMatrixLoader = correlationMatrixLoader;
	}
}
