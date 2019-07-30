package org.reactome.harmonizome;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.dao.GeneCorrelationDAO;
import org.reactome.idg.dao.GeneDAO;
import org.reactome.idg.dao.ProvenanceDAO;
import org.reactome.idg.model.Provenance;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * Loads a gene-pair correlation matrix into a MySQL database.
 * @author sshorser
 *
 */
@Component("correlationMatrixLoader")
public class CorrelationMatrixLoader
{
	private static AtomicInteger fileNum = new AtomicInteger(0);
	
	private int headerSize = 1;
	
	private String delimeter = ",";

	private static final Logger logger = LogManager.getLogger();
	
	private Map<String, Long> symbolToId = new HashMap<>();
	
	private static final String PATH_TO_TMP_FILE = "/tmp/data_for_idg";

	@Autowired
	private int chunkSize;
	
	private String filePath;
	
	@Autowired
	private GeneCorrelationDAO dao;
	
	@Autowired
	private GeneDAO geneDao;
	
	@Autowired
	private ProvenanceDAO provenanceDao;
	
	private Provenance provenance;
	
	public CorrelationMatrixLoader()
	{
		
	}
	
	public CorrelationMatrixLoader(String filePath)
	{
		this.filePath = filePath;
	}
	
	public void loadData() throws Exception
	{
		try(FileInputStream fis = new FileInputStream(this.filePath);
			Scanner scanner = new Scanner(fis);)
		{
			// Using a HashSet ensures that in each chunk, there are duplicate pairs.
			Set<String> lineBuffer = new HashSet<>(chunkSize);
			
			if (scanner.hasNext())
			{
				String geneSymbolHeader = scanner.nextLine();
				String[] parts = geneSymbolHeader.split(delimeter);
				String[] allGeneSymbols = new String[parts.length-1];
				Map<String, Long> allGenesMap = geneDao.getSymbolToIdMapping();
				List<String> symbolsToAdd = new ArrayList<>();
				// start at headerSize, because headers go in BOTH directions. If there is more than 1 header row, there will
				// also be more than 1 column, so you must start extracting gene symbols headerSize elements from the beginning of the array.
				for (int i = headerSize; i < parts.length; i++)
				{
					String gene = parts[i].replaceAll("\"", "");
					allGeneSymbols[i - headerSize] = gene;
					if (!allGenesMap.containsKey(gene))
					{
						symbolsToAdd.add(gene);
					}
				}
				
				logger.info("{} gene symbols in the header, {} are new and will be added to the database.", allGeneSymbols.length, symbolsToAdd.size());
				geneDao.addGenes(symbolsToAdd);
				this.symbolToId = geneDao.getSymbolToIdMapping();
				int maxPairs = (allGeneSymbols.length * (allGeneSymbols.length + 1))/2;
//				this.chunkSize = maxPairs; // to force a single-file dump & load. This could be something that is configurable flag: "single-file-data-load" or something...
				logger.info("{} possible gene-pair correlations.",  maxPairs);
				// global key buffer ensures that we don't add duplicate keys to the database.
				Set<String> globalKeyBuffer = new TreeSet<>();
				final LocalDateTime startTime = LocalDateTime.now();
				AtomicInteger linesInFile = new AtomicInteger(0);
				
				String tempFileName = PATH_TO_TMP_FILE;
				
				Provenance provenanceToUse = provenanceDao.addProvenance(this.provenance);
				AtomicInteger recordCount = new AtomicInteger(0);
				// lineStartOffset will be incremented for each line, moving the "start" pointer to the right. This way, we only capture the "upper" half of the matrix.
				// We don't need to capture the lower half because this is a symmetric matrix.
				int lineStartOffset = this.headerSize;
				
				// Read the rest of the header, but don't consume it. The gene symbols are always the first row, other rows may not 
				// be populated (such as containing "NA" for all Uniprot or NCBI identifiers)
				for (int i = 1; i < headerSize; i++)
				{
					scanner.nextLine();
				}
				// now read the actual lines.
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					parts = line.split(delimeter);
					// if there are 3 header rows/columns, then the offset is 2 (because the arrays are 0-indexed), to skip the header and go straight to the data.
					int headerOffset = headerSize - 1;
					
					String currentGeneSymbol = parts[0];

					int lineWidth = parts.length - lineStartOffset;
					final int startIndex = headerOffset ;
					final int endIndex = lineWidth;
					int i = -1; // initialised to SOMEthing because it gets read outside the scope of the for-loop; in the catch-block.
					try
					{
						// Iterate through the values in the defined range.
						for (i = startIndex; i < parts.length - 1; i++)
						{
							int lookupIndex = i - headerOffset;
							String correlationValue = "";
							correlationValue = parts[i + 1];
							String otherGeneSymbol = allGeneSymbols[lookupIndex];
							String key = generateKey(currentGeneSymbol, otherGeneSymbol);
							String keyParts[] = key.split("\\|");
							boolean newKey = globalKeyBuffer.add(key);
							// if the key is new to this Loader, then add it to the file...
							if (newKey)
							{
								if (this.symbolToId.get(keyParts[0]) == null)
								{
									logger.warn("Map has no ID for symbol {}", keyParts[0]);
								}
								if (this.symbolToId.get(keyParts[1]) == null)
								{
									logger.warn("Map has no ID for symbol {}", keyParts[1]);
								}
								String g1 = this.symbolToId.get(keyParts[0]).toString();
								String g2 = this.symbolToId.get(keyParts[1]).toString();
								
								int count = recordCount.getAndIncrement();
								lineBuffer.add(""+g1+"\t"+g2+"\t"+correlationValue+"\t"+provenanceToUse.getId()+"\n");
								int lineCount = linesInFile.incrementAndGet();
								if (count % 1000000 == 0)
								{
									logger.info("{} million gene-pairs...", count/1000000);
								}
								// When the file is as big as chunkSize, load it to the database. In my experience, 1,000,000 seems to be a good number of
								// records for a single bulk-load. Too big of a chunk sice and the number of inserts/second seems to drop a little.
								if (lineCount % chunkSize == 0)
								{
									loadFileToDatabase(lineBuffer, tempFileName);
								}
							}
						}
					}
					catch (ArrayIndexOutOfBoundsException e)
					{
						logger.info("ArrayOutOfBounds: startIndex: {} endIndex: {} i: {}", startIndex, endIndex, i);
						e.printStackTrace();
						throw e;
					}
					catch (Exception e)
					{
						e.printStackTrace();
						throw e;
					}
				}
				lineStartOffset++;
				// Now process the remainders.
				logger.info("{} gene-pairs will be written to tmp file.", lineBuffer.size());
				if (Files.notExists(Paths.get(tempFileName)))
				{
					writeCorrelationsToFile(lineBuffer, tempFileName);
				}
				// load the last file, probably is smaller than CHUNK_SIZE.
				dao.loadGenePairsFromDataFile(tempFileName);
				Files.move(Paths.get(tempFileName), Paths.get(tempFileName + "_" + fileNum.getAndIncrement()));
				LocalDateTime endTime = LocalDateTime.now();
				logger.info("{} time spent loading the data.", Duration.between(startTime, endTime).toString());
			}
		}
	}

	private void loadFileToDatabase(Set<String> lineBuffer, String tempFileName) throws IOException
	{
		writeCorrelationsToFile(lineBuffer, tempFileName);
		// Now it's time to load the file to the database.
		dao.loadGenePairsFromDataFile(tempFileName);
		// Shuffle the files - current file gets renamed via a move operation, instead of simply being over-written.
		// This way, if something goes wrong, you have ALL of the files on disk.
		Files.move(Paths.get(tempFileName), Paths.get(tempFileName + "_" + fileNum.getAndIncrement()));
		lineBuffer.clear();
	}

	/**
	 * Writes the correlations to a file.
	 * @param lineBuffer - a buffer containing the lines.
	 * @param tempFileName - the name of the file to write.
	 * @throws IOException
	 */
	private static void writeCorrelationsToFile(Set<String> lineBuffer, String tempFileName) throws IOException
	{
		// Sort the buffer. MySQL bulk insert runs faster if the rows are pre-sorted w.r.t. key/index fields.
		// Ideally, the *entire* data set would be sorted before we start chunking, but that might not be possible
		// if the set is large and memory is limited. Which is why we are chunking the data in the first place. :/
		for (String outLine : lineBuffer.parallelStream().sorted().collect(Collectors.toList()))
		{
			Files.write(Paths.get(tempFileName), outLine.getBytes() , StandardOpenOption.WRITE, StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		}
	}


	private static String generateKey(String symbol1, String symbol2)
	{
		return symbol1.compareTo(symbol2) <= 0
				? symbol1 + "|" + symbol2
				: symbol2 + "|" + symbol1;
	}
	
	public Provenance getProvenance()
	{
		return provenance;
	}

	public void setProvenance(Provenance provenance)
	{
		this.provenance = provenance;
	}

	public String getDelimeter()
	{
		return delimeter;
	}

	public void setDelimeter(String delimeter)
	{
		this.delimeter = delimeter;
	}

	public int getHeaderSize()
	{
		return headerSize;
	}

	public void setHeaderSize(int headerSize)
	{
		this.headerSize = headerSize;
	}

	public String getFilePath()
	{
		return filePath;
	}

	public void setFilePath(String filePath)
	{
		this.filePath = filePath;
	}
}
