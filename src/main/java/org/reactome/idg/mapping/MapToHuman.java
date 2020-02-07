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
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
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
	private static final String PPIS_MAPPED_TO_HUMAN_FILE = "PPIS_mapped_to_human.txt";

	private String outputPath = "";
	private static URI uniprotMappingServiceURI = null;

	private String stringDBSpeciesCode;

	public MapToHuman() throws URISyntaxException
	{
		try
		{
			uniprotMappingServiceURI = new URI("https://www.uniprot.org/uploadlists/");
		}
		catch (URISyntaxException e)
		{
			e.printStackTrace();
			throw e;
		}
	}

	public static void main(String[] args) throws IOException, URISyntaxException
	{
		MapToHuman mapper = new MapToHuman();
		mapper.stringDBSpeciesCode = "4932";
		mapper.mapOtherSpeciesToHuman("YEAST", true);
		mapper.stringDBSpeciesCode = "4896";
//		mapper.mapOtherSpeciesToHuman("SCHPO", true);

	}

	private void mapOtherSpeciesToHuman(String speciesName, boolean allowBidirectionalMappings) throws FileNotFoundException, IOException, URISyntaxException
	{
		this.outputPath = "output/"+speciesName+"_results/";
		Files.createDirectories(Paths.get(this.outputPath));
		// Mapping of UniProt identifiers from the Other species to Human, according to PANTHER.
		Map<String, Set<String>> otherSpeciesMappedToHuman = new HashMap<>();
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
					String humanIdentifier = extractUniprotIdentifier(humanParts);
					String otherSpeciesIdentifier = extractUniprotIdentifier(otherSpeciesParts);
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
//						otherSpeciesMappedToHuman.computeIfAbsent(otherSpeciesIdentifier, identifier -> { Set<String> set = new HashSet<>(); set.add(identifier); return set; }  );
					}
				}
				line = br.readLine();
			}
		}
		logger.info("Number of keys in mapping: {}", otherSpeciesMappedToHuman.size());
		logger.info("Number of values in mapping: {}", otherSpeciesMappedToHuman.values().stream().map(Set::size).reduce( 0, (t , u) -> t + u ));
		logger.info("Number of gene families (human): {}", humanGeneFamilies.size());
		logger.info("Number of gene families (other species): {}", otherSpeciesGeneFamilies.size());

		Map<String, Set<String>> extraMappings = generateMappingsFromCommonGeneFamilies(otherSpeciesMappedToHuman, humanGeneFamilies, otherSpeciesGeneFamilies);
		for (Entry<String, Set<String>> extraMapping : extraMappings.entrySet())
		{
			if (otherSpeciesMappedToHuman.containsKey(extraMapping.getKey()))
			{
				otherSpeciesMappedToHuman.get(extraMapping.getKey()).addAll(extraMapping.getValue());
			}
			else
			{
				otherSpeciesMappedToHuman.put(extraMapping.getKey(), extraMapping.getValue());
			}
		}

		logger.info("Number of keys in mapping: {}", otherSpeciesMappedToHuman.size());
		logger.info("Number of values in mapping: {}", otherSpeciesMappedToHuman.values().stream().map(Set::size).reduce( 0, (t , u) -> t + u ));

		// Now we need to process StringDB files.
		int experimentsAndDBScore = 0;
		String stringDBProteinActionsFile = PATH_TO_DATA_FILES + this.stringDBSpeciesCode + ".protein.actions.v11.0.txt";
		String stringDBProteinLinksFile = PATH_TO_DATA_FILES + this.stringDBSpeciesCode + ".protein.links.full.v11.0.txt";
		String ppisWithExperimentsScoreFile = outputPath + this.stringDBSpeciesCode + "_PPIs_with_experiments.tsv";
		Set<String> interactionsWithExperiments = new HashSet<>();
		// First we have to filter for interactions in protein.links.full where experiment > 0
		try (FileReader reader = new FileReader(stringDBProteinLinksFile);
				CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
										.withDelimiter(' ')
										.withFirstRecordAsHeader());
				FileWriter writer = new FileWriter(ppisWithExperimentsScoreFile);)
		{
			List<CSVRecord> records = parser.getRecords();
			logger.info("{} records will be parsed.", records.size());
			for (CSVRecord  record : records)
			{
				int experiments = Integer.parseInt(record.get("experiments"));
				int dbScore = Integer.parseInt(record.get("database"));
				if (experiments > 0 /*|| dbScore > 0*/)
				{
					String protein1 = record.get("protein1").replace(this.stringDBSpeciesCode+".","");
					String protein2 = record.get("protein2").replace(this.stringDBSpeciesCode+".","");
					// For StringDB S. Pombe, some additional cleanup is needed...
					if (this.stringDBSpeciesCode.equals("4896"))
					{
						String[] fixedProteins = MapToHuman.fixSPombeProteins(protein1, protein2);
						protein1 = fixedProteins[0];
						protein2 = fixedProteins[1];
					}
					interactionsWithExperiments.add(putProteinsInOrder(protein1, protein2));
					writer.write(protein1 + "\t" + protein2+"\n");
				}
			}
		}
		logger.info("{} interactions had experiments. " /*+ experimentsAndDBScore + " had DB Score > 0 AND Experiments Score > 0 AND mapped to HUMAN."*/, interactionsWithExperiments.size() );

		// Before we make mapping calls to UniProt, let's see what we can get out of the SGD file.
		MAPPING_TYPE from = MAPPING_TYPE.UNIPROT_GENE_NAME;
		MAPPING_TYPE to = MAPPING_TYPE.UNIPROT_ACCESSION;
		Map<String, Set<String>> mappingsFromSGDFile = getSGDMappings(from, to);

		logger.info("Mappings from SGD ({} -> {}) have {} keys", from, to, mappingsFromSGDFile.size());

		Map<String, Set<String>> uniProtGeneNameToAccessionMapping = mapToUniProtAccessions(this.stringDBSpeciesCode, mappingsFromSGDFile, interactionsWithExperiments, UniprotDB.UniprotGeneName);
		// Before mapping from StringDB to UniProt Accession, prepend the species code to the identifier.
		Map<String, Set<String>> uniProtStringDBsToAccessionMapping = mapToUniProtAccessions(this.stringDBSpeciesCode, mappingsFromSGDFile, interactionsWithExperiments.parallelStream().map(identifier -> this.stringDBSpeciesCode + "."  + identifier).collect(Collectors.toSet()), UniprotDB.STRINGDB);



		//re-assign, using the mapped values.
//		interactionsWithExperiments = mappedUniProtIdentifiers.keySet().stream().map(k -> putProteinsInOrder(k, mappedUniProtIdentifiers.get(k))).collect(Collectors.toSet());
		logger.info("number of gene names (from interactions with experiments > 0) that mapped to accessions: {}", uniProtGeneNameToAccessionMapping.size());
		logger.info("number of StringDB identifiers (from interactions with experiments > 0) that mapped to accessions: {}", uniProtStringDBsToAccessionMapping.size());

		// Merge all the maps.
		Map<String, Set<String>> globalMap = new HashMap<>();
		uniProtGeneNameToAccessionMapping.keySet().stream().forEach(k ->
		{
			globalMap.merge(k, uniProtGeneNameToAccessionMapping.get(k), (s,t) -> {t.addAll(s); return t;});
		});

		uniProtStringDBsToAccessionMapping.keySet().stream().forEach(k ->
		{
			globalMap.merge(k, uniProtStringDBsToAccessionMapping.get(k), (s,t) -> {t.addAll(s); return t;});
		});

		mappingsFromSGDFile.keySet().stream().forEach(k ->
		{
			globalMap.merge(k, mappingsFromSGDFile.get(k), (s,t) -> {t.addAll(s); return t;});
		});

		logger.info("number of keys in global mapping: {}", globalMap.size());

		try(FileReader reader = new FileReader(stringDBProteinActionsFile);
				FileWriter writer = new FileWriter(outputPath + speciesName + "_" + PPIS_MAPPED_TO_HUMAN_FILE);
				FileWriter mappingFailureReasons = new FileWriter(outputPath + speciesName + "_mapping_errors_details.log");
				FileWriter unmappedWriter = new FileWriter(outputPath + speciesName + "_unmapped_from_stringdb.txt");
				FileWriter noMappingToHumanWriter = new FileWriter(outputPath + speciesName + "_no_mapping_to_human.txt"))
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
				Set<String> mappedPPILines = new HashSet<>();
				int selfInteractions = 0;
				for (CSVRecord record : records)
				{
					boolean useStringDBMapping = false;
					if (record.get("mode").equals("binding"))
					{
						numBinding++;
						String protein1 = record.get("item_id_a").replace(this.stringDBSpeciesCode + ".", "");
						String protein2 = record.get("item_id_b").replace(this.stringDBSpeciesCode + ".", "");
						if (this.stringDBSpeciesCode.equals("4896"))
						{
							String[] fixedProteins = MapToHuman.fixSPombeProteins(protein1, protein2);
							protein1 = fixedProteins[0];
							protein2 = fixedProteins[1];
						}

						// Check if this PPI has experiment score > 0 (it will be in interactionsWithExperiments)
						if (interactionsWithExperiments.contains(putProteinsInOrder(protein1, protein2)))
						{
							boolean protein1InGlobalMap = globalMap.containsKey(protein1) && globalMap.get(protein1).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);
							boolean protein2InGlobalMap = globalMap.containsKey(protein2) && globalMap.get(protein2).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);;

							// TODO: cleanup this mess...
//							boolean protein1MapsGeneNameToAccession = uniProtGeneNameToAccessionMapping.containsKey(protein1)
//																		&& uniProtGeneNameToAccessionMapping.get(protein1).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);
//
//							boolean protein2MapsGeneNameToAccession = uniProtGeneNameToAccessionMapping.containsKey(protein2)
//																		&& uniProtGeneNameToAccessionMapping.get(protein2).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);


//							boolean protein1MapsStringDBToAccession = uniProtStringDBsToAccessionMapping.containsKey(this.stringDBSpeciesCode + "." + protein1)
//																		&& uniProtStringDBsToAccessionMapping.get(this.stringDBSpeciesCode + "." + protein1).stream()
//																											.anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein.replace(stringDBSpeciesCode + ".", "")));
//
//
//							boolean protein2MapsStringDBToAccession = uniProtStringDBsToAccessionMapping.containsKey(this.stringDBSpeciesCode + "." + protein2)
//																		&& uniProtStringDBsToAccessionMapping.get(this.stringDBSpeciesCode + "." + protein2).stream()
//																											.anyMatch(protein -> otherSpeciesMappedToHuman.containsKey(protein.replace(stringDBSpeciesCode + ".", "")));

//							boolean protein1MapsViaSGDMApping = mappingsFromSGDFile.containsKey(protein1)
//																&& mappingsFromSGDFile.get(protein1).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);
//
//							boolean protein2MapsViaSGDMApping = mappingsFromSGDFile.containsKey(protein2)
//																&& mappingsFromSGDFile.get(protein2).stream().anyMatch(otherSpeciesMappedToHuman::containsKey);

//							boolean proteinsAreInMapping = (protein1MapsGeneNameToAccession || protein1MapsStringDBToAccession || protein1MapsViaSGDMApping)
//															&& (protein2MapsGeneNameToAccession || protein2MapsStringDBToAccession || protein2MapsViaSGDMApping);

							boolean proteinsAreInMapping = protein1InGlobalMap && protein2InGlobalMap;
//							if ((protein1MapsStringDBToAccession /*&& !protein1MapsGeneNameToAccession*/)
//								|| (protein2MapsStringDBToAccession /*&& !protein2MapsGeneNameToAccession*/))
//							{
//								useStringDBMapping = true;
//								numMappingsInStringDBOnly++;
//							}

							numBindingAndExperimentsGt0++;
							if (proteinsAreInMapping)
							{
								numPPIsInMapping++;
								// Ok we need to get the accession that is also in the PANTHER species mapping...
								String protein1Accession = getUniProtAccessionForProtein(otherSpeciesMappedToHuman, globalMap, uniProtStringDBsToAccessionMapping, protein1);
								String protein2Accession = getUniProtAccessionForProtein(otherSpeciesMappedToHuman, globalMap, uniProtStringDBsToAccessionMapping, protein2);

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
									mappedPPILines.add( putProteinsInOrder(protein1Accession + " (Human: " + String.join(",", humanProteins1.stream().sorted().collect(Collectors.toList())) + ")",
										 						protein2Accession + " (Human: " + String.join(",", humanProteins2.stream().sorted().collect(Collectors.toList())) + ")" ) + "\n" );
								}
							}
							// DEBUG: let's get some REASONS!
							else
							{
								mappingFailureReasons.write("PPI ("+putProteinsInOrder(protein1, protein2)+") could not be fully mapped because:" );
								logErrorsWithProtein(uniProtGeneNameToAccessionMapping, mappingFailureReasons, unmappedIdentifiers, noHumanMappedIdentifiers, npeIdentifiers, otherSpeciesMappedToHuman, protein1);
								logErrorsWithProtein(uniProtGeneNameToAccessionMapping, mappingFailureReasons, unmappedIdentifiers, noHumanMappedIdentifiers, npeIdentifiers, otherSpeciesMappedToHuman, protein2);
								mappingFailureReasons.write('\n');
							}
						}
					}
				}
				logger.info("{} PPIs were binding.", numBinding );
				logger.info("{} PPIs had experiments > 0 AND were binding.", numBindingAndExperimentsGt0);
				logger.info("{} PPIs were in the mapping", numPPIsInMapping);
				logger.info("{} self-interactions were omitted from output.", selfInteractions);
				for (String line : mappedPPILines.stream().sorted().collect(Collectors.toList()))
				{
					writer.write(line);
				}
				logger.info("Number of \"binding\" StringDB interactions (PPIs) that mapped between species: {}; out of a total of {} records.", mappedPPILines.size(), parser.getRecordNumber());
				logger.info("Number of proteins that were only mapped as StringDB-to-UniProt Accession: {}", numMappingsInStringDBOnly);
				logger.info("Number of proteins that could not map to UniProt accessions: {}", unmappedIdentifiers.size());
				for (String unmapped : unmappedIdentifiers)
				{
					unmappedWriter.write(unmapped + "\n");
				}
				logger.info("Number of proteins that could not map to human: {}", noHumanMappedIdentifiers.size());
				for (String unmapped : noHumanMappedIdentifiers)
				{
					noMappingToHumanWriter.write(unmapped + "\n");
				}
				logger.info("Number of proteins identifies that caused NPEs: {}", npeIdentifiers.size());
			}


//			logger.info("StrignDB Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedInteractorsCount);
		}
	}

	private void logErrorsWithProtein(Map<String, Set<String>> uniProtGeneNameToAccessionMapping, FileWriter mappingFailureReasons, Set<String> unmappedIdentifiers, Set<String> noHumanMappedIdentifiers, Set<String> npeIdentifiers, Map<String, Set<String>> otherSpeciesMappedToHuman, String protein1) throws IOException
	{
		if (!uniProtGeneNameToAccessionMapping.containsKey(protein1))
		{
			mappingFailureReasons.write("\t No gene-to-accession mapping for " + protein1);
			unmappedIdentifiers.add(protein1);
		}
		try
		{
			if (!uniProtGeneNameToAccessionMapping.get(protein1).stream().anyMatch(otherSpeciesMappedToHuman::containsKey))
			{
				mappingFailureReasons.write("\t " + protein1 + " could not map to Human.");
				noHumanMappedIdentifiers.add(protein1);
			}
		}
		catch (NullPointerException e)
		{
			npeIdentifiers.add(protein1);
			mappingFailureReasons.write("\t NPE caught while trying to check for a mapping from " + protein1 + " to Human");
		}
	}

	private String getUniProtAccessionForProtein(Map<String, Set<String>> otherSpeciesMappedToHuman, Map<String, Set<String>> globalMapping, Map<String, Set<String>> uniProtStringDBsToAccessionMapping, String protein)
	{
		// well, this code is rather ugly... but not as bad as before.
		String proteinAccession = null;

		// Let's try the StringDB mappings (also from UniProt) because that requires special identifier manipulation
		if (uniProtStringDBsToAccessionMapping.containsKey(protein))
		{
			proteinAccession = uniProtStringDBsToAccessionMapping.get(this.stringDBSpeciesCode + "." + protein).stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p.replace(this.stringDBSpeciesCode + ".", ""))).findFirst().get();
		}
		else
		{
			try
			{
				proteinAccession = globalMapping.get(protein).stream().filter(p -> otherSpeciesMappedToHuman.containsKey(p)).findFirst().get();
			}
			catch (NoSuchElementException e)
			{
				logger.error("NoSuchElementException caught when trying to get mapping for protein: {}", protein);
			}
		}
		return proteinAccession;
	}

	private Map<String, Set<String>> generateMappingsFromCommonGeneFamilies(Map<String, Set<String>> otherSpeciesMappedToHuman, Map<String, Set<String>> humanGeneFamilies, Map<String, Set<String>> otherSpeciesGeneFamilies)
	{
		Set<String> commonFamilies = otherSpeciesGeneFamilies.keySet().parallelStream().filter(humanGeneFamilies::containsKey).collect(Collectors.toSet());
		Map<String, Set<String>> extraMappings = new HashMap<>();
		int extraMapped = 0;
		logger.info("Number of common gene families: {}", commonFamilies.size());
		// create permutations across gene families...
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
						extraMappings.put(nonHumanGene, mapped);
					}
				}
			}
		}
		logger.info("{} extra mappings were created, based on common gene families.", extraMapped);
		return extraMappings;
	}

	private String extractUniprotIdentifier(String[] humanParts)
	{
		String identifier = null;
		String s;
		boolean done = false;
		int i = 0;
		while(!done && i < humanParts.length)
		{
			s = humanParts[i];
			if (s.startsWith("UniProt"))
			{
				identifier = s.substring(s.indexOf('=') + 1);
				done = true;
			}
			i++;
		}
		return identifier;
	}

	private Map<String, Set<String>> mapToUniProtAccessions(String stringDBSpeciesCode, Map<String, Set<String>> preexistingMappings, Set<String> interactionsWithExperiments, UniprotDB mappingSource) throws IOException
	{
		int interactorsAlreadyMapped = 0;
		// Now map the UniProt Gene Names (from StringDB) to UniProt Accessions
		Set<String> identifiersToMapToUniprot = new HashSet<>();
		for (String interactorPair : interactionsWithExperiments)
		{

			String[] parts = interactorPair.split("\\t");
			String interactor1 = parts[0];
			String interactor2 = parts[1];
			if (preexistingMappings.containsKey(interactor1))
			{
				interactorsAlreadyMapped++;
			}
			else
			{
				identifiersToMapToUniprot.add(interactor1);
			}

			if (preexistingMappings.containsKey(interactor2))
			{
				interactorsAlreadyMapped++;
			}
			else
			{
				identifiersToMapToUniprot.add(interactor2);
			}
		}
		logger.info("{} interactors were already found in a pre-existing mapping, so will not be sent to uniprot. {} will be mapped.", interactorsAlreadyMapped, identifiersToMapToUniprot.size());
		Map<String, Set<String>> uniProtGeneNameToAccessionMapping = getMappingsFromUniProt(identifiersToMapToUniprot, mappingSource);

		try(FileWriter writer = new FileWriter(this.outputPath + stringDBSpeciesCode + "_"+mappingSource+"_mappedToUniProt.tsv"))
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

	protected static Map<String, Set<String>> getMappingsFromUniProt(Set<String> identifiersToMapToUniprot, UniprotDB targetDB)
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
				retriever.setUri(uniprotMappingServiceURI);
				retriever.setDataInputStream(new BufferedInputStream(new ByteArrayInputStream(inputString.getBytes())));

				List<String> dataLines = retriever.downloadAndReturnDataLines();
				logger.info("{} {}-to-Uniprot mappings retrieved.", dataLines.size(), targetDB.toString());
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

	private static String[] fixSPombeProteins(String protein1, String protein2)
	{
		String[] proteins = new String[2];
		if (protein1.endsWith(".1"))
		{
			proteins[0] = protein1.substring(0,protein1.length()-2);
		}

		if (protein2.endsWith(".1"))
		{
			proteins[1] = protein2.substring(0,protein2.length()-2);
		}
		return proteins;
	}

	private enum MAPPING_TYPE { UNIPROT_GENE_NAME, SGD, UNIPROT_ACCESSION; }

	private Map<String, Set<String>> getSGDMappings(MAPPING_TYPE from, MAPPING_TYPE to)
	{
		final String uniProtName = "UniProtKB";
		Map<String, Set<String>> sgdMappings = new HashMap<>();
		// See https://downloads.yeastgenome.org/curation/chromosomal_feature/ to get dbxref.tab
		// SGD file will have identifier such as Q0140 in the "UniProtKB/Swiss-Prot ID" field.
		// This can be used to help map identifiers in the StringDB file to UniProt accession numbers.
		// Q0140 -> P02381 in the SGD file, and this reduces the number of lookups you'll need to do with the UniProt mapping service.
		// You could also get the SGD ID and do lookups with that.
		try(BufferedReader reader = Files.newBufferedReader(Paths.get("src", "main", "resources", "data", "dbxref.tab"));)
		{
			String line = reader.readLine();
			while (line != null)
			{
				String[] lineParts = line.split("\\t");
				if (lineParts[1].equals(uniProtName))
				{
					String uniprotAccession = lineParts[0];
					String uniprotSwissprotID = lineParts[3];
					String sgdID = lineParts[4];
					String fromValue = null, toValue = null;
					switch (from)
					{
						// This seems to be a good place to use Java's case-expressions - when they are no longer considered a "preview" feature.
						case UNIPROT_GENE_NAME: fromValue = uniprotSwissprotID; break;
						case UNIPROT_ACCESSION: fromValue = uniprotAccession; break;
						case SGD: fromValue = sgdID; break;
					}
					switch (to)
					{
						case UNIPROT_GENE_NAME: toValue = uniprotSwissprotID; break;
						case UNIPROT_ACCESSION: toValue = uniprotAccession; break;
						case SGD: toValue = sgdID; break;
					}
					Set<String> existingMappings;
					if (sgdMappings.containsKey(fromValue))
					{
						existingMappings = sgdMappings.get(fromValue);
					}
					else
					{
						existingMappings = new HashSet<>();
					}
					existingMappings.add(toValue);
					sgdMappings.put(fromValue, existingMappings);
				}
				line = reader.readLine();
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return sgdMappings;
	}
}
