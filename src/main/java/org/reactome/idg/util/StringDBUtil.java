package org.reactome.idg.util;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

public class StringDBUtil
{
	private StringDBUtil() {}

	public static Set<String> getBindingPPIs(String stringDBProteinActionsFile, Set<String> interactionsWithExperiments)
	{
		Set<String> interactors = new HashSet<>();
		try(FileReader reader = new FileReader(stringDBProteinActionsFile);
//				FileWriter writer = new FileWriter(STRING_DB_PPIS_FILE);
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

						if (interactionsWithExperiments.contains(putProteinsInOrder(protein1, protein2)) && !protein1.equals(protein2))
						{
							// Ok, now we need to start getting mappings from UniProt ACC. If we find an ENSEMBL identifier that's not in the map, add it, but with
							// the value "<EMPTY>"
//							if (!ensembl2UniprotMappings.containsKey(protein1))
//							{
//								ensembl2UniprotMappings.put(protein1, EMPTY_TOKEN);
//								unMappedCount++;
//							}
//							if (!ensembl2UniprotMappings.containsKey(protein2))
//							{
//								ensembl2UniprotMappings.put(protein2, EMPTY_TOKEN);
//								unMappedCount++;
//							}
//
//							identifiersToMapToUniprot.add(protein1);
//							identifiersToMapToUniprot.add(protein2);
							interactors.add(StringDBUtil.putProteinsInOrder(protein1 , protein2));

							numBinding++;
						}
					}
				}
//				System.out.println("Number of \"binding\" StringDB interactions: "+numBinding + "; out of a total of "+parser.getRecordNumber());
			}
		}
		catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return interactors;
	}

	public static Set<String> getPPIsWithExperiments(int threshold, String stringDBProteinLinksFile)
	{
		Set<String> interactionsWithExperiments = new HashSet<>();
		try (FileReader reader = new FileReader(stringDBProteinLinksFile);
				CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
										.withDelimiter(' ')
										.withFirstRecordAsHeader());)
		{
			parser.forEach( record -> {
				int experiments = Integer.parseInt(record.get("experiments"));
//				int dbScore = Integer.parseInt(record.get("database"));
				if (experiments > threshold/* && dbScore > 0*/)
				{
					String protein1 = record.get("protein1");
					String protein2 = record.get("protein2");

					interactionsWithExperiments.add(putProteinsInOrder(protein1, protein2));

				}
			});

		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return interactionsWithExperiments;
	}

	public static String putProteinsInOrder(String p1, String p2)
	{
		// Use String's compareTo to put the proteins in order, and put a tab between them.
		return p1.compareTo(p2) < 0
				? p1 + "\t" + p2
				: p2 + "\t" + p1;
	}

}
