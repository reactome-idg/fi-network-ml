package org.reactome.idg.overlapanalysis;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.management.relation.RelationServiceNotRegisteredException;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.reactome.idg.util.UniprotFileRetriever;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

public class OverlapAnalysis
{

	private static final String STRING_DB_PPIS_FILE = "string-db-ppis.txt";
	private static final String BIOGRID_PPIS_FILE = "biogrid-ppis.txt";
	private static final String REACTOME_PPIS_FILE = "reactome-ppis.txt";
	private static final String MATRIXDB_PPIS_FILE = "matrixdb-ppis.txt";
	private static final String BIOPLEX_PPIS_FILE = "bioplex-ppis.txt";
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

//		try(BufferedReader br = new BufferedReader(new FileReader(STRING_DB_PPIS_FILE)))
//		{
//			boolean done = false;
//			while (!done)
//			{
//				String line = br.readLine();
//				if (line != null)
//				{
//					String key = line.replace("\n", "");
//					interactionsToDBLists.put(key, "StringDB\t");
//				}
//				else
//				{
//					done = true;
//				}
//			}
//		}
//		try(BufferedReader br = new BufferedReader(new FileReader(MATRIXDB_PPIS_FILE)))
//		{
//			boolean done = false;
//			while (!done)
//			{
//				String line = br.readLine();
//				if (line != null)
//				{
//					String key = line.replace("\n", "");
//					if (!interactionsToDBLists.containsKey(key))
//					{
//						interactionsToDBLists.put(key, "MatrixDB\t");
//					}
//					else if (!interactionsToDBLists.get(line).contains("MatrixDB"))
//					{
//						String tmpVal = interactionsToDBLists.get(line) + "MatrixDB\t";
//						interactionsToDBLists.put(key, tmpVal);
//					}
//				}
//				else
//				{
//					done = true;
//				}
//			}
//		}
//		try(BufferedReader br = new BufferedReader(new FileReader(BIOPLEX_PPIS_FILE)))
//		{
//			boolean done = false;
//			while (!done)
//			{
//				String line = br.readLine();
//				if (line != null)
//				{
//					String key = line.replace("\n", "");
//					if (!interactionsToDBLists.containsKey(key))
//					{
//						interactionsToDBLists.put(key, "Bioplex\t");
//					}
//					else if (!interactionsToDBLists.get(line).contains("Bioplex"))
//					{
//						String tmpVal = interactionsToDBLists.get(line) + "Bioplex\t";
//						interactionsToDBLists.put(key, tmpVal);
//					}
//				}
//				else
//				{
//					done = true;
//				}
//			}
//		}
//		try(BufferedReader br = new BufferedReader(new FileReader(BIOGRID_PPIS_FILE)))
//		{
//			boolean done = false;
//			while (!done)
//			{
//				String line = br.readLine();
//				if (line != null)
//				{
//					String key = line.replace("\n", "");
//					if (!interactionsToDBLists.containsKey(key))
//					{
//						interactionsToDBLists.put(key, "BioGrid\t");
//					}
//					else if (!interactionsToDBLists.get(line).contains("BioGrid"))
//					{
//						String tmpVal = interactionsToDBLists.get(line) + "BioGrid\t";
//						interactionsToDBLists.put(key, tmpVal);
//					}
//				}
//				else
//				{
//					done = true;
//				}
//			}
//		}
		TreeMap<String, Set<String>> sorted = new TreeMap<>(/*new Comparator<String>()
				{
					@Override
					public int compare(String o1, String o2)
					{
						return -1 * Integer.compare(o1.length(), o2.length());
					}
				}*/);
//		Function<Set<String>, String> setToString = set -> {
//			return set.stream().sorted().reduce("", String::concat);
//		};

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

		try(FileWriter fw = new FileWriter("overlaps.txt"))
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


//	private static void extractFromBioPlexFile() throws FileNotFoundException, IOException
//	{
//		String bioPlexV3File1 = "BioPlex_2.3_interactionList.tsv";
//		String bioPlexV3File2 = "BioPlex_unpublishedInteractions_May_2019.tsv";
//		Set<String> interactors = new HashSet<>();
//		try(FileReader reader1 = new FileReader(bioPlexV3File1);
//			FileReader reader2 = new FileReader(bioPlexV3File2);
//				FileWriter writer = new FileWriter(BIOPLEX_PPIS_FILE))
//		{
//			CSVFormat format = CSVFormat.DEFAULT
//								.withDelimiter('\t')
//								.withFirstRecordAsHeader();
//
//			Function<FileReader, Integer> processBioPlexFile = (reader) -> {
//				int numInFile = 0;
//
//				try(CSVParser parser = new CSVParser(reader, format);)
//				{
//					for (CSVRecord record : parser.getRecords())
//					{
//						if (Double.parseDouble(record.get("p(Interaction)")) > 0.75d
//								&& (!record.get("UniprotA").equals("UNKNOWN") && !record.get("UniprotB").equals("UNKNOWN")))
//						{
//							String protein1 = record.get("UniprotA");
//							String protein2 = record.get("UniprotB");
//							numInFile++;
//							interactors.add( putProteinsInOrder(protein1, protein2) );
//						}
//					}
//					System.out.println("Number of Bioplex interactions: "+numInFile + "; out of a total of "+parser.getRecordNumber());
//				}
//				catch (IOException e)
//				{
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//				return numInFile;
//			};
//			processBioPlexFile.apply(reader1);
//			processBioPlexFile.apply(reader2);
//			for (String interactor : interactors)
//			{
//				writer.write(interactor + "\n");
//			}
//		}
//	}

//	private static void extractFromMatrixDBFile() throws FileNotFoundException, IOException
//	{
//		final String uniprotPrefix = "uniprotkb:";
//		Set<String> interactors = new HashSet<>();
//		try(FileReader reader = new FileReader(matrixDBFile);
//				FileWriter writer = new FileWriter(MATRIXDB_PPIS_FILE))
//		{
//			CSVFormat format = CSVFormat.DEFAULT
//					.withDelimiter('\t')
//					.withQuote('`')
//					.withFirstRecordAsHeader();
//			try(CSVParser parser = new CSVParser(reader, format);)
//			{
//				for (CSVRecord record : parser.getRecords())
//				{
////					if (!record.get("Confidence value(s)").equals("-"))
//					{
//						String protein1 = record.get("#ID(s) interactor A");
//						String protein2 = record.get("ID(s) interactor B");
//						// Only include MatrixDB interactions where both identifiers are UniProt.
//						if (protein1.contains(uniprotPrefix) && protein2.contains(uniprotPrefix))
//						{
////							writer.write( putProteinsInOrder(protein1, protein2)+"\n" );
//							interactors.add(putProteinsInOrder(protein1.replace(uniprotPrefix, ""), protein2.replace(uniprotPrefix, "")));
//						}
//					}
//				}
//				for (String interactor : interactors)
//				{
//					writer.write(interactor + "\n");
//				}
//			}
//		}
//
//	}

//	private static void convertStringDBFileToUniProtIdentifiers() throws FileNotFoundException, IOException
//	{
//		String stringDBPPIFile = "string-db-ppis.txt";
//		try(BufferedReader bis = new BufferedReader(new FileReader(stringDBPPIFile)))
//		{
//			Set<String> ensemblIdentifiers = new HashSet<>();
//			boolean done = false;
//			while (!done)
//			{
//				String line = bis.readLine();
//				// the line has two identifiers
//				String[] parts = line.split("\\t");
//				ensemblIdentifiers.add(parts[0]);
//				ensemblIdentifiers.add(parts[1]);
//			}
//		}
//	}


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

//	private static void extractFromBioGridFile() throws FileNotFoundException, IOException, URISyntaxException
//	{
//		String bioGridFile = "BIOGRID-ALL-3.5.176.tab2.txt";
//		List<String> interactors = new ArrayList<>();
//		Map<String, String> entrezGene2UniProt = new HashMap<>();
//		// TAXONOMY ID is 9606
//		final String humanOranismTaxonID = "9606";
//		try(FileReader reader = new FileReader(bioGridFile);
//				FileWriter writer = new FileWriter(BIOGRID_PPIS_FILE))
//		{
//			CSVFormat format = CSVFormat.DEFAULT
//								.withDelimiter('\t')
//								.withFirstRecordAsHeader();
//			try(CSVParser parser = new CSVParser(reader, format);)
//			{
//				List<CSVRecord> records = parser.getRecords();
//				int totalNumRecords = records.size();
//				int unMappedCount = 0;
//				for (CSVRecord record : records)
//				{
//					String interactorA = record.get("Organism Interactor A");
//					String interactorB = record.get("Organism Interactor B");
//
//					if (interactorA != null && interactorA.trim().equals(humanOranismTaxonID)
//						&& interactorB != null && interactorB.trim().equals(humanOranismTaxonID))
//					{
//						String score = record.get("Score");
//						if (score != null && !score.equals("-") && Double.parseDouble(score) > 0.75d)
//						{
//							String entrezGeneA = record.get("Entrez Gene Interactor A");
//							String entrezGeneB = record.get("Entrez Gene Interactor B");
//							if (!entrezGene2UniProt.containsKey(entrezGeneA))
//							{
//								entrezGene2UniProt.put(entrezGeneA, EMPTY_TOKEN);
//								unMappedCount++;
//							}
//							if (!entrezGene2UniProt.containsKey(entrezGeneB))
//							{
//								entrezGene2UniProt.put(entrezGeneB, EMPTY_TOKEN);
//								unMappedCount++;
//							}
//							if (unMappedCount == 1000 || record.getRecordNumber() >= totalNumRecords)
//							{
//								unMappedCount = getMappingsFromUniProt(entrezGene2UniProt, UniprotDB.Entrez_Gene);
//							}
//							interactors.add(entrezGeneA + "|" + entrezGeneB);
//						}
//
//					}
//
//				}
//			}
//			int unmappedCount = 0;
//			Set<String> interactorSet = new HashSet<>();
//			int mappedCount = 0;
//			// Now, we go through the interactors, map them to UniProt and write them to a file.
//			for (String interactorPair : interactors)
//			{
//				String[] parts = interactorPair.split("\\|");
//				String interactorA = parts[0];
//				String interactorB = parts[1];
//
//				String mappedInteractorAId = entrezGene2UniProt.get(interactorA);
//				String mappedInteractorBId = entrezGene2UniProt.get(interactorB);
//
//
//				if (!mappedInteractorAId.equals(EMPTY_TOKEN) && !mappedInteractorBId.equals(EMPTY_TOKEN))
//				{
//					mappedCount++;
//					interactorSet.add(putProteinsInOrder(mappedInteractorAId, mappedInteractorBId));
////					writer.write( putProteinsInOrder(mappedInteractorAId, mappedInteractorBId)+"\n" );
//				}
//				else
//				{
//					unmappedCount ++;
//					if (mappedInteractorAId.equals(EMPTY_TOKEN))
//					{
////						System.out.println("Could not map interactor " + interactorA + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
//					}
//					if (mappedInteractorBId.equals(EMPTY_TOKEN))
//					{
////						System.out.println("Could not map interactor " + interactorB + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
//					}
//				}
//			}
//			for (String interactor : interactorSet)
//			{
//				writer.write( interactor + "\n" );
//			}
//			System.out.println("StrignDB Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedCount);
//		}
//
//	}

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


//	private static void extractFromStringDBFile() throws IOException, FileNotFoundException, URISyntaxException
//	{
//		Map<String, String> ensembl2UniprotMappings = new HashMap<>();
//		List<String> interactors = new ArrayList<>();
//		int mappedCount = 0;
//		int unMappedCount = 0;
//		String stringDBFile = "9606.protein.actions.v11.0.txt";
//
//		try(FileReader reader = new FileReader(stringDBFile);
//				FileWriter writer = new FileWriter(STRING_DB_PPIS_FILE))
//		{
//			CSVFormat format = CSVFormat.DEFAULT
//								.withDelimiter('\t')
//								.withFirstRecordAsHeader();
//			try(CSVParser parser = new CSVParser(reader, format);)
//			{
//				int numBinding = 0;
//				List<CSVRecord> records = parser.getRecords();
//				int totalNumRecords = records.size();
//				for (CSVRecord record : records)
//				{
//					if (record.get("mode").equals("binding"))
//					{
//						String protein1 = record.get("item_id_a").replace("9606.", "");
//						String protein2 = record.get("item_id_b").replace("9606.", "");
//
//						// Ok, now we need to start getting mappings from UniProt ACC. If we find an ENSEMBL identifier that's not in the map, add it, but with
//						// the value "<EMPTY>"
//						if (!ensembl2UniprotMappings.containsKey(protein1))
//						{
//							ensembl2UniprotMappings.put(protein1, EMPTY_TOKEN);
//							unMappedCount++;
//						}
//						if (!ensembl2UniprotMappings.containsKey(protein2))
//						{
//							ensembl2UniprotMappings.put(protein2, EMPTY_TOKEN);
//							unMappedCount++;
//						}
//
//						// every 1000 proteins, we'll get a mapping from UniProt;
//						if (unMappedCount == 1000 || record.getRecordNumber() >= totalNumRecords)
//						{
////							System.out.println("Retrieving UniProt mappings, please wait...");
////							UniprotFileRetriever retriever = new UniprotFileRetriever();
////							retriever.setMapToDbEnum(UniprotDB.UniProtAcc);
////							retriever.setMapFromDbEnum(UniprotDB.ENSEMBLProtein);
////							retriever.setUri(new URI("https://www.uniprot.org/uploadlists/"));
//////							BinaryOperator<String> accumulator = ;
////							byte[] buf = ensembl2UniprotMappings.keySet().stream()
////																	.filter(key -> ensembl2UniprotMappings.get(key).equals(EMPTY_TOKEN))
////																	.reduce("", (t, u) -> t + " " + u).getBytes();
////							ByteArrayInputStream inStream = new ByteArrayInputStream(buf);
////							BufferedInputStream stringStream = new BufferedInputStream(inStream);
////							retriever.setDataInputStream(stringStream );
////
////							List<String> dataLines = retriever.downloadAndReturnDataLines();
////							System.out.println(dataLines.size() + "Ensenbl-to-Uniprot mappings retrieved.");
////							// Just in case...
////							if (dataLines!=null)
////							{
////								for (String line : dataLines)
////								{
////									if (!line.equals("From	To") && !line.equals("not mapped"))
////									{
////										String[] parts = line.split("\\t");
////										String ensemblIDentifier = parts[0];
////										String uniProtIdentifier = parts[1];
////
////										ensembl2UniprotMappings.put(ensemblIDentifier, uniProtIdentifier);
////									}
////								}
////							}
////							unMappedCount = 0;
//							unMappedCount = getMappingsFromUniProt(ensembl2UniprotMappings, UniprotDB.ENSEMBLProtein);
//						}
//						interactors.add(protein1 + "|" + protein2);
//
//						numBinding++;
////						writer.write( putProteinsInOrder(protein1.replace("9606.", ""), protein2.replace("9606.", ""))+"\n" );
//					}
//				}
//				System.out.println("Number of \"binding\" StringDB interactions: "+numBinding + "; out of a total of "+parser.getRecordNumber());
//			}
//			int unmappedCount = 0;
//			Set<String> interactorSet = new HashSet<>();
//			// Now, we go through the interactors, map them to UniProt and write them to a file.
//			for (String interactorPair : interactors)
//			{
//				String[] parts = interactorPair.split("\\|");
//				String interactorA = parts[0];
//				String interactorB = parts[1];
//
//				String mappedInteractorAId = ensembl2UniprotMappings.get(interactorA);
//				String mappedInteractorBId = ensembl2UniprotMappings.get(interactorB);
//
//
//				if (!mappedInteractorAId.equals(EMPTY_TOKEN) && !mappedInteractorBId.equals(EMPTY_TOKEN))
//				{
//					mappedCount++;
//					interactorSet.add(putProteinsInOrder(mappedInteractorAId, mappedInteractorBId));
////					writer.write( putProteinsInOrder(mappedInteractorAId, mappedInteractorBId)+"\n" );
//				}
//				else
//				{
//					unmappedCount ++;
////					if (mappedInteractorAId.equals(EMPTY_TOKEN))
////					{
////						System.out.println("Could not map interactor " + interactorA + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
////					}
////					if (mappedInteractorBId.equals(EMPTY_TOKEN))
////					{
////						System.out.println("Could not map interactor " + interactorB + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
////					}
//				}
//			}
//			for (String interactor : interactorSet)
//			{
//				writer.write( interactor + "\n" );
//			}
//			System.out.println("StrignDB Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedCount);
//		}
//	}

}
