package org.reactome.idg.util;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.http.HttpEntity;
import org.apache.http.HttpStatus;
import org.apache.http.NoHttpResponseException;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.client.utils.URIBuilder;
import org.apache.http.entity.ContentType;
import org.apache.http.entity.mime.MultipartEntityBuilder;
import org.apache.http.entity.mime.content.StringBody;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class UniprotFileRetriever
{
	// To be used to wait 500 ms to retry if a URL from UniProt returns nothing.
	private static final int RETRY_DELAY_MS = 3000;
	private static final int MAX_NUM_ATTEMPTS = 5;
	private String mapFromDb = "";
	private String mapToDb = "";
	private BufferedInputStream inStream;
	// A list of paths that were actually downloaded to.
	private static List<String> actualFetchDestinations = Collections.synchronizedList(new ArrayList<String>());

	private final int maxAttemptCount = 5;
	private URI uri;
	private String destination;
	private Logger logger = LogManager.getLogger();

	/**
	 * This enum provides a mapping between Reactome names for reference databases and the Uniprot ID that is used by their
	 * mapping service.
	 *
	 * @author sshorser
	 *
	 */
	public enum UniprotDB
	{
		// For a list of database IDs that Uniprot can map from, see: https://www.uniprot.org/help/api_idmapping
		OMIM("MIM_ID"), PDB("PDB_ID"), RefSeqPeptide("P_REFSEQ_AC"), RefSeqRNA("REFSEQ_NT_ID"), ENSEMBL("ENSEMBL_ID"), ENSEMBLProtein("ENSEMBL_PRO_ID"), ENSEMBLGenomes("ENSEMBLGENOME_ID"),
		ENSEMBLTranscript("ENSEMBL_TRS_ID"), Ensembl("ENSEMBL_ID"), Wormbase("WORMBASE_ID"), Entrez_Gene("P_ENTREZGENEID"), UniprotGeneName("GENENAME"), KEGG("KEGG_ID"), STRINGDB("STRING_ID"),
		UniProtAcc("ACC"), UniProtID("ID"), UniProt("ACC+ID"), UCSC("UCSC_ID"),
		UniProtGeneName("GENENAME");

		private String uniprotName;
		private static Map<String, UniprotDB> mapToEnum;

		private UniprotDB(String s)
		{
			this.uniprotName = s;
			updateMap(s);
		}

		private void updateMap(String s)
		{
			if (mapToEnum == null)
			{
				mapToEnum = new HashMap<>(10);
			}
			mapToEnum.put(s, this);
		}

		public String getUniprotName()
		{
			return this.uniprotName;
		}

		public static UniprotDB uniprotDBFromUniprotName(String uniprotName)
		{
			return mapToEnum.get(uniprotName);
		}
	}

	public UniprotFileRetriever()
	{
		super();
	}

//	public UniprotFileRetriever(String retrieverName)
//	{
//		super(retrieverName);
//	}

	private URIBuilder uriBuilderFromDataLocation(String location) throws URISyntaxException
	{
		URIBuilder builder = new URIBuilder();
		if (location == null || location.length() == 0)
		{
			logger.error("Location was NULL/0-length!");
			return null;
		}
		String[] parts = location.split("\\?");
		String[] schemeAndHost = parts[0].split("://");
		builder.setScheme(schemeAndHost[0]);
		if (schemeAndHost.length > 1)
		{
			builder.setHost(schemeAndHost[1]);
		}
		else
		{
			logger.warn("schemeAndHost had length < 2: {} ; will use try to host from original URI: {}", schemeAndHost, this.uri.getHost());
			builder.setScheme(this.uri.getScheme());
			builder.setHost(this.uri.getHost() + "/" + schemeAndHost[0]);
		}

		if (parts.length > 1)
		{
			// If the Location header string contains query information, we need to properly reformat that before requesting it.
			String[] params = parts[1].split("&");
			for (String s : params)
			{
				String[] nameAndValue = s.split("=");
				builder.addParameter(nameAndValue[0], nameAndValue[1]);
			}
		}
		else
		{
			// Add .tab to get table.
			if (!builder.getHost().endsWith(".tab"))
				builder.setHost(builder.getHost() + ".tab");
		}
		return builder;
	}

	/**
	 * Attempt to GET data from UniProt.
	 *
	 * @param get - the HttpGet object.
	 * @return A byte array of the result's content.
	 * @throws IOException
	 * @throws URISyntaxException
	 * @throws InterruptedException
	 */
	private byte[] attemptGetFromUniprot(HttpGet get) throws IOException, URISyntaxException, InterruptedException
	{
		byte[] result = null;
		logger.trace("getting from: {}", get.getURI());
		try (CloseableHttpClient getClient = HttpClients.createDefault(); CloseableHttpResponse getResponse = getClient.execute(get);)
		{
			switch (getResponse.getStatusLine().getStatusCode())
			{
			case HttpStatus.SC_SERVICE_UNAVAILABLE:
			case HttpStatus.SC_INTERNAL_SERVER_ERROR:
			case HttpStatus.SC_BAD_GATEWAY:
			case HttpStatus.SC_GATEWAY_TIMEOUT:
				logger.error("Error {} detected! Message: {}", getResponse.getStatusLine().getStatusCode(), getResponse.getStatusLine().getReasonPhrase());
				break;

			case HttpStatus.SC_OK:
			case HttpStatus.SC_MOVED_PERMANENTLY:
			case HttpStatus.SC_MOVED_TEMPORARILY:
				logger.trace("HTTP Status: {}", getResponse.getStatusLine().toString());
				result = EntityUtils.toByteArray(getResponse.getEntity());
				if (result == null)
				{
					logger.warn("Response did not contain data.");
				}
				break;

			default:
				logger.warn("Nothing was downloaded due to an unexpected status code and message: {} / {} ", getResponse.getStatusLine().getStatusCode(), getResponse.getStatusLine());
				break;
			}
		}
		return result;
	}

	/**
	 * Attempt to POST to UniProt. If successful, the URL to the *actual* data will be returned.
	 *
	 * @param post - the POST object.
	 * @return The URL to find the mapped data at, as a string.
	 * @throws ClientProtocolException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private String attemptPostToUniprot(HttpPost post) throws ClientProtocolException, IOException, InterruptedException
	{
		boolean done = false;
		int attemptCount = 0;
		String mappingLocationURI = null;
		while (!done)
		{

			{
				try (CloseableHttpClient postClient = HttpClients.createDefault(); CloseableHttpResponse postResponse = postClient.execute(post);)
				{
					switch (postResponse.getStatusLine().getStatusCode())
					{
					case HttpStatus.SC_SERVICE_UNAVAILABLE:
					case HttpStatus.SC_INTERNAL_SERVER_ERROR:
					case HttpStatus.SC_BAD_GATEWAY:
					case HttpStatus.SC_GATEWAY_TIMEOUT:
						logger.error("Error {} detected! Message: {}", postResponse.getStatusLine().getStatusCode(), postResponse.getStatusLine().getReasonPhrase());
						break;

					case HttpStatus.SC_OK:
					case HttpStatus.SC_MOVED_PERMANENTLY:
					case HttpStatus.SC_MOVED_TEMPORARILY:
						if (postResponse.containsHeader("Location"))
						{
							mappingLocationURI = postResponse.getHeaders("Location")[0].getValue();
							logger.trace("Location of data: {}", mappingLocationURI);
							if (mappingLocationURI != null && !mappingLocationURI.equals("http://www.uniprot.org/502.htm"))
							{
								done = true;
							}
							else
							{
								logger.warn("Response did not contain data.");
							}
						}
						else
						{
							logger.warn("Status was {}, \"Location\" header was not present. Other headers are: {}", postResponse.getStatusLine().toString(),
									Arrays.stream(postResponse.getAllHeaders()).map((h -> h.toString())).collect(Collectors.joining(" ; ")));
						}
						break;
					default:
						logger.warn("Nothing was downloaded due to an unexpected status code and message: {} / {} ", postResponse.getStatusLine().getStatusCode(), postResponse.getStatusLine());
						break;
					}
					attemptCount++;

				}
				if (attemptCount > this.maxAttemptCount)
				{
					logger.error("Reached max attempt count! No more attempts.");
					done = true;
				}
				else
				{
					if (attemptCount < this.maxAttemptCount && !done)
					{
						logger.info("Re-trying... {} attempts made, {} allowed", attemptCount, this.maxAttemptCount);
						Thread.sleep(Duration.ofSeconds(5).toMillis());
					}
				}
			}
		}
		return mappingLocationURI;
	}

	public List<String> downloadAndReturnDataLines()
	{
		try
		{
			this.destination = Files.createTempFile("uniprotData", ".tab").toString();
			this.downloadData();
			// After the download, the destination variable actually points to the "notMapped" file,
			// but this function is only interested in the mapped values, so snip the ".notMapped" out of the filename.
			this.destination = this.destination.replace(".notMapped", "");
			return Files.readAllLines(Paths.get(this.destination));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Getting data from UniProt is a 3-stage process: 1) POST a list of identifiers to UniProt. The response received
	 * contains a URL to the mapped data. 2) GET the data from the URL in the response from 1). 3) GET the "not" mapped data
	 * from the URL in the response from 1).
	 */
	public void downloadData()
	{
		// Check inputs:
		if (this.inStream == null)
		{
			throw new RuntimeException("inStream is null! You must provide an data input stream!");
		}
		else if (this.mapFromDb.trim().length() == 0)
		{
			throw new RuntimeException("You must provide a database name to map from!");
		}
		else if (this.mapToDb.trim().length() == 0)
		{
			throw new RuntimeException("You must provide a database name to map to!");
		}

		try
		{
			String location = getDataLocation();
			if (location != null)
			{
				// Get values that Uniprot was able to map.
				URI mappedUri = createURI(location, true);
				getUniprotValues(mappedUri);
				// Now get the unmapped values.
				URI unmappedUri = createURI(location, false);
				getUniprotValues(unmappedUri);
			}
			else
			{
				logger.error("We could not determine the location of the data, file was not downloaded.");
			}
		}
		catch (URISyntaxException | ClientProtocolException | InterruptedException e)
		{
			logger.error("A problem occured while trying to get the data location, or the data: {}", e.getMessage());
			e.printStackTrace();
		}
		catch (IOException e)
		{
			// some of the exceptions in the first catch-block are also covered by IOException. They are in the
			// first catch block so we can treat them differently, if necessary.
			logger.error("IOException was caught: {}", e.getMessage());
			e.printStackTrace();
		}
		catch (Exception e)
		{
			logger.error("Exception occurred: {}", e.getMessage());
			e.printStackTrace();
		}
	}

	/**
	 * This function gets the location of the Uniprot-mapped data.
	 *
	 * @return The location of the data, as a string.
	 * @throws IOException
	 * @throws MalformedURLException
	 * @throws ClientProtocolException
	 * @throws InterruptedException
	 */
	private String getDataLocation() throws IOException, MalformedURLException, ClientProtocolException, InterruptedException
	{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		byte[] buffer = new byte[2048];
		int len;
		while ((len = inStream.read(buffer)) > -1)
		{
			baos.write(buffer, 0, len);
		}
		String location = null;
		int attemptCount = 0;
		while (location == null && attemptCount < maxAttemptCount)
		{
			InputStream fileData = new ByteArrayInputStream(baos.toByteArray());
			HttpPost post = new HttpPost(this.uri);
			logger.trace("URI: {}", post.getURI().toURL());
			HttpEntity attachment = MultipartEntityBuilder.create().addBinaryBody("file", fileData, ContentType.TEXT_PLAIN, "uniprot_ids.txt")
					.addPart("format", new StringBody("tab", ContentType.MULTIPART_FORM_DATA)).addPart("from", new StringBody(this.mapFromDb, ContentType.MULTIPART_FORM_DATA))
					.addPart("to", new StringBody(this.mapToDb, ContentType.MULTIPART_FORM_DATA)).build();
			post.setEntity(attachment);

			try
			{
				location = this.attemptPostToUniprot(post);
			}
			catch (NoHttpResponseException e)
			{
				// If we don't catch this here, but let it go to "catch (IOException e)" in the outer try-block,
				// then we won't be able to retry. Catching it here lets us continue processing: increment the attempt counter, and loop
				// through again.
				logger.error("No HTTP Response! Message: {}", e.getMessage());
				e.printStackTrace();
			}
			attemptCount++;
			if (location == null)
			{
				if (attemptCount < maxAttemptCount)
				{
					Random r = new Random(System.nanoTime());
					long delay = (long) (3000 + (attemptCount * r.nextFloat()));
					logger.warn("Attempt {} out of {}, next attempt in {} ms", attemptCount, maxAttemptCount, delay);
					Thread.sleep(delay);
				}
				else if (attemptCount >= maxAttemptCount)
				{
					logger.error("Could not get the Location of the data in {} attempts.", attemptCount);
				}
			}

		}
		return location;
	}

	/**
	 * Get values from Uniprot.
	 *
	 * @param location - the URL that the data will be at. This is returned from the inital query to UniProt.
	 * @param mapped   - Set to true to get the list of successfully mapped values. Set to false to get the list of
	 *                 identifiers which UniProt could not map.
	 * @throws URISyntaxException
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws Exception
	 */
	private void getUniprotValues(URI uri) throws URISyntaxException, IOException, InterruptedException, Exception
	{
		int numAttempts = 0;
		boolean done = false;

		while (!done)
		{
			HttpGet get = new HttpGet(uri);
			byte[] result = this.attemptGetFromUniprot(get);
			numAttempts++;
			Path path = Paths.get(new URI("file://" + this.destination));
			// The loop is "done" if the number of attempts exceeds the max allowed, OR if the result is not null/empty.
			done = (numAttempts >= MAX_NUM_ATTEMPTS) || (result != null && result.length > 0);
			if (result != null && result.length > 0)
			{
				logger.debug(".tab result size: {}", result.length);
				Files.createDirectories(path.getParent());
				// Files.write(path, result);
				BufferedWriter writer = Files.newBufferedWriter(path);
				writer.write(new String(result));
				writer.flush();
				writer.close();
				if (!Files.isReadable(path))
				{
					throw new Exception("The new file " + path + " is not readable!");
				}
				UniprotFileRetriever.actualFetchDestinations.add(path.toString());
			}
			else
			{
				handleNullResult(numAttempts, done, path);
			}
		}
	}

	/**
	 * Create the URI to get data from. <code>mapped</code> is used to determine if ".not" will be in the URI, which is used
	 * to get identifiers which were not mapped.
	 *
	 * @param location - The URL to the data, for the mapped values. UniProt will return this by default.
	 * @param mapped   - Set to true if you want the URI for mapped values. Set to false if you want values that UniProt
	 *                 couldn't map. <b>NOTE:</b> Setting <code>mapped</code> to false will also modify
	 *                 <code>this.destination</code> to include "notMapped" in the filename.
	 * @return A URI that can be used to get data from UniProt.
	 * @throws URISyntaxException
	 */
	private URI createURI(String location, boolean mapped) throws URISyntaxException
	{
		URI uri;
		if (mapped)
		{
			uri = this.uriBuilderFromDataLocation(location).build();
		}
		else
		{
			URIBuilder builder = uriBuilderFromDataLocation(location);
			uri = builder.setHost(builder.getHost().replace(".tab", ".not")).build();
			String[] filenameParts = this.destination.split("\\.");
			this.destination = this.destination.replace(filenameParts[filenameParts.length - 1], "notMapped." + filenameParts[filenameParts.length - 1]);
		}
		return uri;
	}

	/**
	 * Handle a NULL result.
	 *
	 * @param numAttempts - The number of attempts that have been made (so far), only needed for logging.
	 * @param done        - Are we done? If this is true, an exception will be thrown since the maximum number of attempts
	 *                    have been made and the result is still null. If false, a log message will be printed and this
	 *                    thread will sleep for RETRY_DELAY_MS milliseconds.
	 * @param path        - the path to the file that we were *trying* to download. Used for logging purposes.
	 * @throws InterruptedException
	 * @throws Exception
	 */
	private void handleNullResult(int numAttempts, boolean done, Path path) throws InterruptedException, Exception
	{
		if (!done)
		{
			logger.info("Result was NULL. Will sleep for a bit, and try again. {} attempts remaining, {} total are allowed.", MAX_NUM_ATTEMPTS - numAttempts, MAX_NUM_ATTEMPTS);
			// Sleep a half a second...
			Thread.sleep(RETRY_DELAY_MS);
		}
		else
		{
			throw new Exception("Result for .tab file (" + path.toString() + ") was null/empty!");
		}
	}

	public String getMapFromDb()
	{
		return this.mapFromDb;
	}

	public String getMapToDb()
	{
		return this.mapToDb;
	}

	public void setMapFromDbEnum(UniprotDB mapFromDb)
	{
		this.mapFromDb = mapFromDb.getUniprotName();
	}

	public void setMapToDbEnum(UniprotDB mapToDb)
	{
		this.mapToDb = mapToDb.getUniprotName();
	}

	public void setMapFromDb(String mapFromDb)
	{
		this.mapFromDb = mapFromDb;
	}

	public void setMapToDb(String mapToDb)
	{
		this.mapToDb = mapToDb;
	}

	public void setDataInputStream(BufferedInputStream inStream)
	{
		this.inStream = inStream;
	}

	public String getFetchDestination()
	{
		return this.destination;
	}

	public List<String> getActualFetchDestinations()
	{
		return UniprotFileRetriever.actualFetchDestinations;
	}

	public URI getUri()
	{
		return uri;
	}

	public void setUri(URI uri)
	{
		this.uri = uri;
	}
}