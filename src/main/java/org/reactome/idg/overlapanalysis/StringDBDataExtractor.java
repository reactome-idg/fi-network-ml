package org.reactome.idg.overlapanalysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.util.Species;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

public class StringDBDataExtractor extends DataExtractor
{
	private static final String STRING_DB_PPIS_FILE = "output/overlaps/string-db-ppis.txt";
	private static final Logger logger = LogManager.getLogger(StringDBDataExtractor.class);

	public void extractFromStringDBFile() throws IOException, FileNotFoundException, URISyntaxException
	{
		System.out.println("Extracting StringDB Data...");
		Map<String, String> ensembl2UniprotMappings = new HashMap<>();
		List<String> interactors = new ArrayList<>();
		int mappedCount = 0;
		int unMappedCount = 0;
		int experimentsAndDBScore = 0;
		String stringDBProteinActionsFile = "src/main/resources/data/9606.protein.actions.v11.0.txt";
		String stringDBProteinLinksFile = "src/main/resources/data/9606.protein.links.full.v11.0.txt";
		Set<String> interactionsWithExperiments = new HashSet<>();
		// First we have to filter for interactions in protein.links.full where experiment > 0
		try (FileReader reader = new FileReader(stringDBProteinLinksFile);
				CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
										.withDelimiter(' ')
										.withFirstRecordAsHeader());)
		{
			List<CSVRecord> records = parser.getRecords();
			System.out.println(records.size() + " records will be parsed.");
			for (CSVRecord  record : records)
			{
				int experiments = Integer.parseInt(record.get("experiments"));
				int dbScore = Integer.parseInt(record.get("database"));
				if (experiments > 0/* && dbScore > 0*/)
				{
					String protein1 = record.get("protein1");
					String protein2 = record.get("protein2");

					interactionsWithExperiments.add(putProteinsInOrder(protein1, protein2));

					// Guanming was curious about this:
					if (Integer.parseInt(record.get("database")) > 0)
					{
						experimentsAndDBScore ++;
					}
				}
			}
		}
		// TODO: the MapToHuman project tries to get the same data from StringDB files, so: extract the common code to a new "StringDBUtilies" class.
		System.out.println(interactionsWithExperiments.size() + " interactions had experiments. " + experimentsAndDBScore + " had DB Score > 0 AND Experiments Score > 0.");
		Set<String> identifiersToMapToUniprot = new HashSet<>();
		try(FileReader reader = new FileReader(stringDBProteinActionsFile);
				FileWriter writer = new FileWriter(STRING_DB_PPIS_FILE);
				FileWriter unmappedWriter = new FileWriter("unmapped_from_stringdb.txt"))
		{
			Set<String> unmappedIdentifiers = new HashSet<>();
			CSVFormat format = CSVFormat.DEFAULT
								.withDelimiter('\t')
								.withFirstRecordAsHeader();
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				int numBinding = 0;
				List<CSVRecord> records = parser.getRecords();
				int totalNumRecords = records.size();
				for (CSVRecord record : records)
				{
					if (record.get("mode").equals("binding"))
					{
						String protein1 = record.get("item_id_a");
						String protein2 = record.get("item_id_b");

						if (interactionsWithExperiments.contains(putProteinsInOrder(protein1, protein2)))
						{
							// Ok, now we need to start getting mappings from UniProt ACC. If we find an ENSEMBL identifier that's not in the map, add it, but with
							// the value "<EMPTY>"
							if (!ensembl2UniprotMappings.containsKey(protein1))
							{
								ensembl2UniprotMappings.put(protein1, EMPTY_TOKEN);
								unMappedCount++;
							}
							if (!ensembl2UniprotMappings.containsKey(protein2))
							{
								ensembl2UniprotMappings.put(protein2, EMPTY_TOKEN);
								unMappedCount++;
							}

							identifiersToMapToUniprot.add(protein1);
							identifiersToMapToUniprot.add(protein2);
							interactors.add(protein1 + "|" + protein2);

							numBinding++;
						}
					}
				}
				System.out.println("Number of \"binding\" StringDB interactions: "+numBinding + "; out of a total of "+parser.getRecordNumber());
			}
			ensembl2UniprotMappings.clear();
			// TODO: StringDB provides its own mapping files, use those instead of the UniProt service (StringDB's file will provide better mappings)
			ensembl2UniprotMappings = getMappingsFromUniProt(identifiersToMapToUniprot , UniprotDB.ENSEMBLProtein);
			int unmappedInteractorsCount = 0;
			Set<String> interactorSet = new HashSet<>();
			// Now, we go through the interactors, map them to UniProt and write them to a file.
			for (String interactorPair : interactors)
			{
				String[] parts = interactorPair.split("\\|");
				String interactorA = parts[0];
				String interactorB = parts[1];

				String mappedInteractorAId = ensembl2UniprotMappings.get(interactorA);
				String mappedInteractorBId = ensembl2UniprotMappings.get(interactorB);

				if (mappedInteractorAId != null && mappedInteractorBId!=null &&
					!mappedInteractorAId.equals(EMPTY_TOKEN) && !mappedInteractorBId.equals(EMPTY_TOKEN))
				{
					mappedCount++;
					interactorSet.add(putProteinsInOrder(mappedInteractorAId, mappedInteractorBId));
				}
				else
				{
					if (mappedInteractorAId == null)
					{
						unmappedIdentifiers.add(interactorA);
						unmappedInteractorsCount++;
					}
					if (mappedInteractorBId == null)
					{
						unmappedIdentifiers.add(interactorB);
						unmappedInteractorsCount++;
					}
				}
			}
			for (String interactor : interactorSet)
			{
				writer.write( interactor + "\n" );
			}
			for (String unmappedIdentifier : unmappedIdentifiers)
			{
				unmappedWriter.write(unmappedIdentifier + "\n");
			}
			System.out.println("StrignDB Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedInteractorsCount);
		}
	}

	private Map<String, Set<String>> getStringDBToUniProtMappings(String pathToMappings, Species species)
	{
		// TODO: Move this method to a StringDBUtilies class
		Map<String, Set<String>> mappings = new HashMap<>();
		// StringDB has a mappings file that shows StringDB -> UniProt
		try(BufferedReader br = new BufferedReader(new FileReader(pathToMappings)))
		{
			String line = br.readLine();
			// skip the header line.
			line = br.readLine();
			while (line != null)
			{
				String[] parts = line.split("\\t");
				// first field is species, we should filter on that.
				if (Integer.parseInt(parts[0]) == species.getSpeciesCode())
				{
					// 2nd and 3rd columns are the ones we're interested in.
					String stringDBIdentifier = parts[2];
					String uniProtIdentifier = parts[1].split("\\|")[0];

					Set<String> s = mappings.containsKey(stringDBIdentifier) ? mappings.get(stringDBIdentifier) : new HashSet<>();
					s.add(uniProtIdentifier);
					mappings.put(stringDBIdentifier, s);
				}
				line = br.readLine();
			}
		}
		catch (IOException e)
		{
			logger.error("File Error!", e);
		}
		return mappings;
	}


	@Override
	protected Map<String, String> getMappingsFromUniProt(Set<String> identifiersToMapToUniprot, UniprotDB targetDB) throws URISyntaxException
	{
		Map<String, Set<String>> mappings = this.getStringDBToUniProtMappings("src/main/resources/data/all_organisms.uniprot_2_string.2018.tsv", Species.HUMAN);
		Map<String, String> simpleMappings = new HashMap<>(mappings.size());
		for (Entry<String, Set<String>> mapping : mappings.entrySet())
		{
			simpleMappings.put(mapping.getKey(), mapping.getValue().stream().findFirst().get());
		}
		return simpleMappings;
	}

}
