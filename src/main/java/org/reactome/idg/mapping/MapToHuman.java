package org.reactome.idg.mapping;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.util.UniprotFileRetriever;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

public class MapToHuman
{
	private static final Logger logger = LogManager.getLogger(MapToHuman.class);
	private static final String HUMAN = "HUMAN";
	private static final String PATH_TO_DATA_FILES = "src/main/resources/data/";
	//	private static final String EMPTY_TOKEN = "<EMTPY>";
	private static final String PPIS_MAPPED_TO_HUMAN_FILE = "PPIS_mapped_to_human.txt";
	public static void main(String[] args) throws IOException, URISyntaxException
	{
		mapOtherSpeciesToHuman("YEAST", "4932", true);
//		mapOtherSpeciesToHuman("SCHPO", "4896", true);
	}


	private static void mapOtherSpeciesToHuman(String speciesName, String stringDBSpeciesCode, boolean allowBidirectionalMappings) throws FileNotFoundException, IOException, URISyntaxException
	{
//		boolean allowBidirectionalMappings = false;
		// Mapping of UniProt identifiers from the Other species to Human, according to PANTHER.
		Map<String, Set<String>> otherSpeciesMappedToHuman = new HashMap<>();
		// String species name
//		String speciesName = "SCHPO";//"YEAST"; //
		// String species code - used in StringDB file
//		String stringDBSpeciesCode = "4896";//"4932"; //
		// Path to input StringDB file.
		String stringDBFile;
		// Path to input Panther file. Obtain the file here: ftp://ftp.pantherdb.org/ortholog/current_release/Orthologs_HCOP.tar.gz
		// Panther file contains UniProt ACCESSIONS - the StringDB file for YEAST contains UniProt GENE NAMES!!
		String pantherFile = PATH_TO_DATA_FILES + "Orthologs_HCOP";
		// Output file: non-human PPI's mapped to human (via Panther mappings).
		String mappedToHumanFile = speciesName + "_mapped_to_human.tsv";
		// Read the Panther file - it's big so maybe filter out the lines we want and rewrite those to a temp file.
		// Then, read the StringDB file and extract the PPIs for the species we're interested in.
		// Use the Panther mappings to map the PPIs.
		int lineCount = 0;
		Map<String,Set<String>> humanGeneFamilies = new HashMap<>();
		Map<String,Set<String>> otherSpeciesGeneFamilies = new HashMap<>();
		try(BufferedReader br = new BufferedReader(new FileReader(pantherFile)))
		{
			String line = br.readLine();
//			lineCount++;
//			if (lineCount % 100000 == 0)
//			{
//				logger.info(lineCount + " lines read from PANTHER Orthologs.");
//			}
			while (line != null)
			{
				String[] parts = line.split("\\t");
				// Only the first two parts are useful, and only if they contain HUMAN or YEAST (or whatever speciesName is).
				String species1 = (parts[0].split("\\|"))[0];
				String species2 = (parts[1].split("\\|"))[0];
				String geneFamily = parts[4];
				boolean lineIsUseful = !species1.equals(species2) && ((allowBidirectionalMappings && species1.equals(HUMAN) && species2.equals(speciesName))
																	|| (species1.equals(speciesName) && species2.equals(HUMAN)));
//				boolean lineIsUseful = ((parts[0].startsWith("HUMAN") && parts[1].startsWith(speciesName)) || (parts[0].startsWith(speciesName) && parts[1].startsWith("HUMAN")));
				if (lineIsUseful)
				{
					String[] humanParts;
					String[] otherSpeciesParts;

					if (species1.equals(HUMAN))
					{
						humanParts = parts[0].split("\\|");
						otherSpeciesParts = parts[1].split("\\|");
					}
					else
					{
						humanParts = parts[1].split("\\|");
						otherSpeciesParts = parts[0].split("\\|");
					}
					// Now... get the identifier values for Human and the other species.
					String s = "";
					boolean done = false;
					int i = 0;
					String humanIdentifier = null;
					String otherSpeciesIdentifier = null;
					while(!done && i < humanParts.length)
					{
						s = humanParts[i];
						if (s.startsWith("UniProt"))
						{
							humanIdentifier = s.substring(s.indexOf('=') + 1);
							done = true;
						}
						i++;
					}
					done = false;
					i = 0;
					while(!done && i < otherSpeciesParts.length)
					{
						s = otherSpeciesParts[i];
						if (s.startsWith("UniProt"))
						{
							otherSpeciesIdentifier = s.substring(s.indexOf('=') + 1);
							done = true;
						}
						i++;
					}
					// Now, set up the gene families map. This will be used to create permutations of gene mappings (from $SPECIES -> HUMAN).
					if (line.startsWith(HUMAN))
					{
						Set<String> genes ;
						if (humanGeneFamilies.containsKey(geneFamily))
						{
							genes = humanGeneFamilies.get(geneFamily);
						}
						else
						{
							genes = new HashSet<>();
						}
						genes.add(humanIdentifier);
						humanGeneFamilies.put(geneFamily, genes);
					}
					else
					{
						Set<String> genes ;
						if (otherSpeciesGeneFamilies.containsKey(geneFamily))
						{
							genes = otherSpeciesGeneFamilies.get(geneFamily);
						}
						else
						{
							genes = new HashSet<>();
						}
						genes.add(otherSpeciesIdentifier);
						otherSpeciesGeneFamilies.put(geneFamily, genes);
					}

					if (otherSpeciesIdentifier == null || humanIdentifier == null)
					{
						logger.error("Error, one or more identifiers could not be mapped to UniProt!");
					}
					else
					{
						Set<String> identifiers ;
						if (otherSpeciesMappedToHuman.containsKey(otherSpeciesIdentifier))
						{
							identifiers = otherSpeciesMappedToHuman.get(otherSpeciesIdentifier);
						}
						else
						{
							identifiers = new HashSet<>();
						}
						identifiers.add(humanIdentifier);
						otherSpeciesMappedToHuman.put(otherSpeciesIdentifier, identifiers);
						otherSpeciesMappedToHuman.computeIfAbsent(otherSpeciesIdentifier, identifier -> { Set<String> set = new HashSet<>(); set.add(identifier); return set; }  );
					}
				}
				line = br.readLine();
			}
		}
		logger.info("Number of keys in mapping: " + otherSpeciesMappedToHuman.size());
		logger.info("Number of values in mapping: " + otherSpeciesMappedToHuman.values().stream().map( set -> set.size()).reduce( 0, (t , u) -> {return t + u;} ));

		logger.info("Number of gene families (human): " + humanGeneFamilies.size());
		logger.info("Number of gene families (other species): " + otherSpeciesGeneFamilies.size());

		Set<String> commonFamilies = otherSpeciesGeneFamilies.keySet().parallelStream().filter(fam -> humanGeneFamilies.containsKey(fam)).collect(Collectors.toSet());

		int extraMapped = 0;
		logger.info("Number of common gene families: " + commonFamilies.size());
		// TODO: create permutations across gene families...
		for (String commonFamily : commonFamilies)
		{
			Set<String> humanGenes = humanGeneFamilies.get(commonFamily);
			Set<String> nonHumanGenes = otherSpeciesGeneFamilies.get(commonFamily);
			// Now create a mapping for everything in nonHumanGenes to humanGenes
			for (String nonHumanGene : nonHumanGenes)
			{
				for (String humanGene : humanGenes)
				{
					Set<String> mapped = otherSpeciesMappedToHuman.get(nonHumanGene);
					if(mapped.add(humanGene))
					{
						extraMapped++;
						otherSpeciesMappedToHuman.put(nonHumanGene, mapped);
					}
				}
			}
		}
		logger.info(extraMapped + " extra mappings were created, based on common gene families.");
		logger.info("Number of keys in mapping: " + otherSpeciesMappedToHuman.size());
		logger.info("Number of values in mapping: " + otherSpeciesMappedToHuman.values().stream().map( set -> set.size()).reduce( 0, (t , u) -> {return t + u;} ));

		// Now we need to process StringDB files.
		int experimentsAndDBScore = 0;
		String stringDBProteinActionsFile = PATH_TO_DATA_FILES + stringDBSpeciesCode + ".protein.actions.v11.0.txt";
		String stringDBProteinLinksFile = PATH_TO_DATA_FILES + stringDBSpeciesCode + ".protein.links.full.v11.0.txt";
		String ppisWithExperimentsScoreFile = stringDBSpeciesCode + "_PPIs_with_experiments.tsv";
		Set<String> interactionsWithExperiments = new HashSet<>();
		// First we have to filter for interactions in protein.links.full where experiment > 0
		try (FileReader reader = new FileReader(stringDBProteinLinksFile);
				CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
										.withDelimiter(' ')
										.withFirstRecordAsHeader());
				FileWriter writer = new FileWriter(ppisWithExperimentsScoreFile);)
		{
			List<CSVRecord> records = parser.getRecords();
			logger.info(records.size() + " records will be parsed.");
			for (CSVRecord  record : records)
			{
				int experiments = Integer.parseInt(record.get("experiments"));
				int dbScore = Integer.parseInt(record.get("database"));
				if (experiments > 0 /*|| dbScore > 0*/)
				{
					String protein1 = record.get("protein1").replace(stringDBSpeciesCode+".","");
					String protein2 = record.get("protein2").replace(stringDBSpeciesCode+".","");

					// For StringDB S. Pombe, some additional cleanup is needed...
					if (stringDBSpeciesCode.equals("4896"))
					{
						if (protein1.endsWith(".1"))
						{
							protein1 = protein1.substring(0,protein1.length()-2);
						}

						if (protein2.endsWith(".1"))
						{
							protein2 = protein2.substring(0,protein2.length()-2);
						}
					}

					interactionsWithExperiments.add(putProteinsInOrder(protein1, protein2));
					writer.write(protein1 + "\t" + protein2+"\n");
				}
			}
		}
		logger.info(interactionsWithExperiments.size() + " interactions had experiments. " /*+ experimentsAndDBScore + " had DB Score > 0 AND Experiments Score > 0 AND mapped to HUMAN."*/);

		Map<String, Set<String>> uniProtGeneNameToAccessionMapping = mapToUniProtAccessions(stringDBSpeciesCode, interactionsWithExperiments, UniprotDB.UniprotGeneName);
		// Before mapping from StringDB to UniProt Accession, prepend the species code to the identifier.
		Map<String, Set<String>> uniProtStringDBsToAccessionMapping = mapToUniProtAccessions(stringDBSpeciesCode, interactionsWithExperiments.parallelStream().map(identifier -> { return stringDBSpeciesCode + "."  + identifier;}).collect(Collectors.toSet()), UniprotDB.STRINGDB);


		//re-assign, using the mapped values.
//		interactionsWithExperiments = mappedUniProtIdentifiers.keySet().stream().map(k -> putProteinsInOrder(k, mappedUniProtIdentifiers.get(k))).collect(Collectors.toSet());
		logger.info("number of gene names (from interactions with experiments > 0) that mapped to accessions: " + uniProtGeneNameToAccessionMapping.size());
		logger.info("number of StringDB identifiers (from interactions with experiments > 0) that mapped to accessions: " + uniProtStringDBsToAccessionMapping.size());
//		int unMappedCount = 0, mappedCount = 0;
//		List<String> interactors = new ArrayList<>();
//		Map<String, String> ensembl2UniprotMappings = new HashMap<>();
//		Set<String> identifiersToMapToUniprot = new HashSet<>();
		try(FileReader reader = new FileReader(stringDBProteinActionsFile);
				FileWriter writer = new FileWriter(speciesName + "_" + PPIS_MAPPED_TO_HUMAN_FILE);
				FileWriter unmappedWriter = new FileWriter("unmapped_from_stringdb.txt"))
		{
			Set<String> unmappedIdentifiers = new HashSet<>();
			Set<String> noHumanMappedIdentifiers = new HashSet<>();
			Set<String> npeIdentifiers = new HashSet<>();

			CSVFormat format = CSVFormat.DEFAULT
								.withDelimiter('\t')
								.withFirstRecordAsHeader();
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				int numBinding = 0;
				int numBindingAndExperimentsGt0 = 0;
				int numPPIsInMapping = 0;
				int numMappingsInStringDBOnly = 0;
				List<CSVRecord> records = parser.getRecords();
				int totalNumRecords = records.size();
				Set<String> outLines = new HashSet<>();
				int selfInteractions = 0;
				for (CSVRecord record : records)
				{
					boolean useStringDBMapping = false;
					if (record.get("mode").equals("binding"))
					{
						numBinding++;
						String protein1 = record.get("item_id_a").replace(stringDBSpeciesCode + ".", "");
						String protein2 = record.get("item_id_b").replace(stringDBSpeciesCode + ".", "");
						if (stringDBSpeciesCode.equals("4896"))
						{
							if (protein1.endsWith(".1"))
							{
								protein1 = protein1.substring(0,protein1.length()-2);
							}

							if (protein2.endsWith(".1"))
							{
								protein2 = protein2.substring(0,protein2.length()-2);
							}
						}
						// Check if this PPI has experiment score > 0 (it will be in interactionsWithExperiments)
						if (interactionsWithExperiments.contains(putProteinsInOrder(protein1, protein2)))
						{
							boolean protein1MapsGeneNameToAccession = uniProtGeneNameToAccessionMapping.containsKey(protein1)
																		&& uniProtGeneNameToAccessionMapping.get(protein1).stream().anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein));

							boolean protein2MapsGeneNameToAccession = uniProtGeneNameToAccessionMapping.containsKey(protein2)
																		&& uniProtGeneNameToAccessionMapping.get(protein2).stream().anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein));


							boolean protein1MapsStringDBToAccession = uniProtStringDBsToAccessionMapping.containsKey(stringDBSpeciesCode + "." + protein1)
																		&& uniProtStringDBsToAccessionMapping.get(stringDBSpeciesCode + "." + protein1).stream()
																											.anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein.replace(stringDBSpeciesCode + ".", "")));


							boolean protein2MapsStringDBToAccession = uniProtStringDBsToAccessionMapping.containsKey(stringDBSpeciesCode + "." + protein2)
																		&& uniProtStringDBsToAccessionMapping.get(stringDBSpeciesCode + "." + protein2).stream()
																											.anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein.replace(stringDBSpeciesCode + ".", "")));


							boolean proteinsAreInMapping = (protein1MapsGeneNameToAccession || protein1MapsStringDBToAccession)
															&& (protein2MapsGeneNameToAccession || protein2MapsStringDBToAccession);

							if ((protein1MapsStringDBToAccession && !protein1MapsGeneNameToAccession)
								|| (protein2MapsStringDBToAccession && !protein2MapsGeneNameToAccession))
							{
								useStringDBMapping = true;
								numMappingsInStringDBOnly++;
							}

							numBindingAndExperimentsGt0++;
							if (proteinsAreInMapping)
							{
								numPPIsInMapping++;
								// Ok we need to get the accession that is also in the PANTHER species mapping...
								String protein1Accession = protein1MapsGeneNameToAccession
															? uniProtGeneNameToAccessionMapping.get(protein1).stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p)).findFirst().get()
															: uniProtStringDBsToAccessionMapping.get(stringDBSpeciesCode + "." + protein1)
																								.stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p.replace(stringDBSpeciesCode + ".", ""))).findFirst().get();

								String protein2Accession = protein2MapsGeneNameToAccession
															? uniProtGeneNameToAccessionMapping.get(protein2).stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p)).findFirst().get()
															: uniProtStringDBsToAccessionMapping.get(stringDBSpeciesCode + "." + protein2)
																								.stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p.replace(stringDBSpeciesCode + ".", ""))).findFirst().get();

								// Self-interactions occur when two of the non-human proteins both map to the exact same human protein (and they do not map to any other proteins).
								// These are ignored in the output. When two non-human proteins map to multiple human proteins where some of them are the same,
								// it's OK since there are permutations of those human proteins that can exclude the self-interaction.
								Set<String> humanProteins1 = otherSpeciesMappedToHuman.get(protein1Accession);
								Set<String> humanProteins2 = otherSpeciesMappedToHuman.get(protein2Accession);
								boolean selfInteraction = false;
								if (humanProteins1.size() == 1 && humanProteins2.size() == 1)
								{
									if (humanProteins1.stream().findFirst().equals(humanProteins2.stream().findFirst()))
									{
										selfInteraction = true;
										selfInteractions++;
									}
								}
								if (!selfInteraction)
								{
									outLines.add( putProteinsInOrder(protein1Accession + " (Human: " + String.join(",", humanProteins1.stream().sorted().collect(Collectors.toList())) + ")",
										 						protein2Accession + " (Human: " + String.join(",", humanProteins2.stream().sorted().collect(Collectors.toList())) + ")" ) + "\n" );
								}
							}
							// DEBUG: let's get some REASONS!
							else
							{
								logger.info("PPI ("+putProteinsInOrder(protein1, protein2)+") could not be fully mapped because: ");
								if (!uniProtGeneNameToAccessionMapping.containsKey(protein1))
								{
									logger.info("  No gene-to-accession mapping for " + protein1);
									unmappedIdentifiers.add(protein1);
								}
								if (!uniProtGeneNameToAccessionMapping.containsKey(protein2))
								{
									logger.info("  No gene-to-accession mapping for " + protein2);
									unmappedIdentifiers.add(protein2);
								}
								try
								{
									if (!uniProtGeneNameToAccessionMapping.get(protein1).stream().anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein)))
									{
										logger.info("  " + protein1 + " could not map to Human.");
										noHumanMappedIdentifiers.add(protein1);
									}
								}
								catch (NullPointerException e)
								{
									npeIdentifiers.add(protein1);
									logger.info("  NPE caught while trying to check for a mapping from " + protein1 + " to Human");
								}
								try
								{
									if (!uniProtGeneNameToAccessionMapping.get(protein2).stream().anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein)))
									{
										logger.info("  " + protein2 + " could not map to Human.");
										noHumanMappedIdentifiers.add(protein2);
									}
								}
								catch (NullPointerException e)
								{
									npeIdentifiers.add(protein2);
									logger.info("  NPE caught while trying to check for a mapping from " + protein2 + " to Human");
								}
							}
						}
					}
				}
				logger.info(numBinding + " PPIs were binding.");
				logger.info(numBindingAndExperimentsGt0 + " PPIs had experiments > 0 AND were binding.");
				logger.info(numPPIsInMapping + " PPIs were in the mapping");
				logger.info(selfInteractions + " self-interactions were omitted from output.");
				for (String line : outLines.stream().sorted().collect(Collectors.toList()))
				{
					writer.write(line);
				}
				logger.info("Number of \"binding\" StringDB interactions (PPIs) that mapped between species: "+outLines.size() + "; out of a total of "+parser.getRecordNumber() + " records.");
				logger.info("Number of proteins that were only mapped as StringDB-to-UniProt Accession: " + numMappingsInStringDBOnly);
				logger.info("Number of proteins that could not map to UniProt accessions: " + unmappedIdentifiers.size());
				for (String unmapped : unmappedIdentifiers)
				{
					unmappedWriter.write(unmapped + "\n");
				}
				logger.info("Number of proteins that could not map to human: " + noHumanMappedIdentifiers.size());
				logger.info("Number of proteins identifies that caused NPEs: " + npeIdentifiers.size());
			}


//			logger.info("StrignDB Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedInteractorsCount);
		}
	}

	private static Map<String, Set<String>> mapToUniProtAccessions(String stringDBSpeciesCode, Set<String> interactionsWithExperiments, UniprotDB mappingSource) throws URISyntaxException, IOException
	{
		// Now map the UniProt Gene Names (from StringDB) to UniProt Accessions
		Set<String> identifiersToMapToUniprot = new HashSet<>();
		for (String interactorPair : interactionsWithExperiments)
		{
			String[] parts = interactorPair.split("\\t");
			String interactor1 = parts[0];
			String interactor2 = parts[1];
			identifiersToMapToUniprot.add(interactor1);
			identifiersToMapToUniprot.add(interactor2);
		}
		Map<String, Set<String>> uniProtGeneNameToAccessionMapping = getMappingsFromUniProt(identifiersToMapToUniprot, mappingSource);

		try(FileWriter writer = new FileWriter(stringDBSpeciesCode + "_mappedToUniProt.tsv"))
		{
			for (Entry<String, Set<String>> entry : uniProtGeneNameToAccessionMapping.entrySet())
			{
				String mapping = String.join(",", entry.getValue().stream().sorted().collect(Collectors.toList()));
				writer.write(entry.getKey() + "\t" + mapping+"\n");
			}
		}
		return uniProtGeneNameToAccessionMapping;
	}

	protected static String putProteinsInOrder(String p1, String p2)
	{
		// Use String's compareTo to put the proteins in order, and put a tab between them.
		return p1.compareTo(p2) < 0
				? p1 + "\t" + p2
				: p2 + "\t" + p1;
	}

	protected static Map<String, Set<String>> getMappingsFromUniProt(Set<String> identifiersToMapToUniprot, UniprotDB targetDB) throws URISyntaxException
	{
		Map<String, Set<String>> identifierToUniprotMap = new HashMap<>();
		int i = 0, j = 0;
		int k = 2000;
		String inputString = "";
		for (String identifierToMap : identifiersToMapToUniprot)
		{
			i++;
			j++;
			inputString += identifierToMap + " ";

			if (i % k == 0 || identifiersToMapToUniprot.size() == j)
			{
				UniprotFileRetriever retriever = new UniprotFileRetriever();
				retriever.setMapToDbEnum(UniprotDB.UniProtAcc);
				retriever.setMapFromDbEnum(targetDB);
				retriever.setUri(new URI("https://www.uniprot.org/uploadlists/"));
				retriever.setDataInputStream(new BufferedInputStream(new ByteArrayInputStream(inputString.getBytes())));

				List<String> dataLines = retriever.downloadAndReturnDataLines();
				logger.info(dataLines.size() + " " + targetDB.toString() + "-to-Uniprot mappings retrieved.");
				// Just in case...
				if (dataLines!=null)
				{
					for (String line : dataLines)
					{
						if (!line.equals("From	To") && !line.equals("not mapped"))
						{
							String[] parts = line.split("\\t");
							String sourceIDentifier = parts[0];
							String uniProtIdentifier = parts[1];
							Set<String> set;
							if (!identifierToUniprotMap.containsKey(sourceIDentifier))
							{
								set = new HashSet<>();
							}
							else
							{
								set = identifierToUniprotMap.get(sourceIDentifier);
							}
							set.add(uniProtIdentifier);
							identifierToUniprotMap.put(sourceIDentifier, set);
						}
					}
				}
				i = 0;
				inputString = "";
			}
		}

		return identifierToUniprotMap;
	}
}
