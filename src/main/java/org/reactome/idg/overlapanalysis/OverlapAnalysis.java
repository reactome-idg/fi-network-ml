package org.reactome.idg.overlapanalysis;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.reactome.idg.util.UniprotFileRetriever;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

public class OverlapAnalysis
{

	private static final String STRING_DB_PPIS_FILE = "output/overlaps/string-db-ppis.txt";
	private static final String BIOGRID_PPIS_FILE = "output/overlaps/biogrid-ppis.txt";
	private static final String REACTOME_PPIS_FILE = "output/overlaps/reactome-ppis.txt";
	private static final String MATRIXDB_PPIS_FILE = "output/overlaps/matrixdb-ppis.txt";
	private static final String BIOPLEX_PPIS_FILE = "output/overlaps/bioplex-ppis.txt";
	private static final String EMPTY_TOKEN = "<EMPTY>";
//	private static final String matrixDBFile = "matrixdb_FULL.tab";

	private static String putProteinsInOrder(String p1, String p2)
	{
		// Use String's compareTo to put the proteins in order, and put a tab between them.
		return p1.compareTo(p2) < 0
				? p1 + "\t" + p2
				: p2 + "\t" + p1;
	}

	public static void main(String[] args) throws IOException, URISyntaxException
	{
		System.out.println("Extracting interactions from various sources...");
//		extractFromStringDBFile();
//		MatrixDBDataExtractor.extractFromMatrixDBFile();
//		BioPlexDataExtractor.extractFromBioPlexFile();
		BioGridDataExtractor bioGridExtractor = new BioGridDataExtractor();
		bioGridExtractor.extractFromBioGridFile();
		StringDBDataExtractor stringDBExtractor = new StringDBDataExtractor();
		stringDBExtractor.extractFromStringDBFile();
//		extractFromReactomeFile();
		// I wish Java had an nicer way to create a Map *and* populate it at the same time.
		// I think this is the best there is, unless we migrate beyond Java 8
		Map<String, String> files = new HashMap<String, String>() {{
			put("StringDB", STRING_DB_PPIS_FILE);
//			put("MatrixDB", MATRIXDB_PPIS_FILE);
//			put("Bioplex", BIOPLEX_PPIS_FILE);
			put("Biogrid", BIOGRID_PPIS_FILE);
//			put("Reactome", REACTOME_PPIS_FILE);
		}};
		analyzeOverlap(files);

//		Number of interactions in "StringDB":	387052
//		Number of interactions in "StringDB, MatrixDB":	9330

//		Number of interactions in "Bioplex":	117326
//		Number of interactions in "StringDB":	380830
//		Number of interactions in "StringDB, MatrixDB, Bioplex":	1392
//		Number of interactions in "MatrixDB, Bioplex":	6222
//		Number of interactions in "StringDB, MatrixDB":	7938
	}

	public static void analyzeOverlap(Map<String, String> interactorFiles) throws IOException
	{
		System.out.println("Generating overlap file");

		// In this map, the KEY is the interaction itself, and the value is a string containing the names of databases where the interaction is found.
		Map<String, Set<String>> interactionsToDBLists = new HashMap<>();

		for (String dbName : interactorFiles.keySet())
		{
			try(BufferedReader br = new BufferedReader(new FileReader(interactorFiles.get(dbName))))
			{
				boolean done = false;
				while (!done)
				{
					String line = br.readLine();
					if (line != null)
					{
						String interaction = line.replace("\n", "");
						if (!interactionsToDBLists.containsKey(interaction))
						{
							interactionsToDBLists.put(interaction, new HashSet<>(Arrays.asList(dbName + "\t")));
						}
						else if (!interactionsToDBLists.get(interaction).contains(dbName))
						{
							Set<String> dbSet = interactionsToDBLists.get(interaction);
							dbSet.add(dbName + "\t");
							interactionsToDBLists.put(interaction, dbSet);
						}
					}
					else
					{
						done = true;
					}
				}
			}
		}

		TreeMap<String, Set<String>> sorted = new TreeMap<>();

		//Now, sort the array in order of most overlap to least. Sort by length of the VALUE that the keys point to.
//		interactionsToDBLists.keySet().stream().forEach(key -> sorted.computeIfAbsent(setToString.apply(interactionsToDBLists.get(key)), k -> new ArrayList<>()).add(key) );

		// For each interaction, add to "sorted tree...
		interactionsToDBLists.keySet().stream().forEach( interaction -> {
//			String db = setToString.apply(interactionsToDBLists.get(interaction));
			String db = (interactionsToDBLists.get(interaction)).stream().sorted().reduce("", String::concat);
			Set<String> tmpSet;
			if (sorted.containsKey(db))
			{
				tmpSet = sorted.get(db);
			}
			else
			{
				tmpSet = new HashSet<>();
			}
			tmpSet.add(interaction);
			sorted.put(db, tmpSet);
		});

		try(FileWriter fw = new FileWriter("output/overlaps/overlaps.txt"))
		{
			sorted.keySet().stream().forEach(k -> sorted.get(k).stream().forEach(interaction -> {
				try
				{
					fw.write(k + interaction + "\n");
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}));
		}

		Map<String, Integer> summaryCounts = new HashMap<>();
		// now for a summary report:
		sorted.keySet().stream().sorted().forEach(k -> summaryCounts.put(k, sorted.get(k).size()));

		summaryCounts.keySet().stream()
			.sorted(new Comparator<String>()
			{
				@Override
				public int compare(String o1, String o2)
				{
					return Integer.compare(o1.length(), o2.length());
				}
			})
			.forEach(k -> System.out.println("Number of interactions in \""+k.trim().replace("\t", ", ")+"\":\t"+summaryCounts.get(k)));

	}

	private static void extractFromReactomeFile() throws FileNotFoundException, IOException
	{
		final String reactomeFile = "reactome.homo_sapiens.interactions.tab-delimited.txt";
		Set<String> interactions = new HashSet<>();
		try(FileReader reader = new FileReader(reactomeFile);
				FileWriter writer = new FileWriter(REACTOME_PPIS_FILE))
		{
			CSVFormat format = CSVFormat.DEFAULT
					.withDelimiter('\t')
					.withFirstRecordAsHeader();
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				List<CSVRecord> records = parser.getRecords();
//				int totalNumRecords = records.size();
//				int unMappedCount = 0;
				for (CSVRecord record : records)
				{
					String interactor1 = record.get("# Interactor 1 uniprot id");
					String interactor2 = record.get("Interactor 2 uniprot id");
					if (interactor1.contains("uniprotkb:") && interactor2.contains("uniprotkb:"))
					{
						interactions.add(putProteinsInOrder(interactor1.replace("uniprotkb:", "") , interactor2.replace("uniprotkb:", "")) + "\n");
					}
				}
			}
			for (String interaction : interactions)
			{
				writer.write(interaction);
			}
		}
	}

	private static int getMappingsFromUniProt(Map<String, String> mappings, UniprotDB targetDB) throws URISyntaxException
	{
		int unMappedCount;
		System.out.println("Retrieving UniProt mappings, please wait...");
		UniprotFileRetriever retriever = new UniprotFileRetriever();
		retriever.setMapToDbEnum(UniprotDB.UniProtAcc);
		retriever.setMapFromDbEnum(targetDB);
		retriever.setUri(new URI("https://www.uniprot.org/uploadlists/"));

		byte[] buf = mappings.keySet().stream()
								.filter(key -> mappings.get(key).equals(EMPTY_TOKEN))
								.reduce("", (t, u) -> t + " " + u).getBytes();
		ByteArrayInputStream inStream = new ByteArrayInputStream(buf);
		BufferedInputStream stringStream = new BufferedInputStream(inStream);
		retriever.setDataInputStream(stringStream );

		List<String> dataLines = retriever.downloadAndReturnDataLines();
		System.out.println(dataLines.size() + " " + targetDB.toString() + "-to-Uniprot mappings retrieved.");
		// Just in case...
		if (dataLines!=null)
		{
			for (String line : dataLines)
			{
				if (!line.equals("From	To") && !line.equals("not mapped"))
				{
					String[] parts = line.split("\\t");
					String ensemblIDentifier = parts[0];
					String uniProtIdentifier = parts[1];

					mappings.put(ensemblIDentifier, uniProtIdentifier);
				}
			}
		}
		unMappedCount = 0;
		return unMappedCount;
	}
}
