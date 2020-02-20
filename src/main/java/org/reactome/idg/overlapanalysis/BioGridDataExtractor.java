package org.reactome.idg.overlapanalysis;

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
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

public class BioGridDataExtractor extends DataExtractor
{
	private static final String BIOGRID_PPIS_FILE = "biogrid-ppis.txt";

	public void extractFromBioGridFile() throws FileNotFoundException, IOException, URISyntaxException
	{
		System.out.println("Extracting BioGrid Data...");
		String bioGridFile = "src/main/resources/data/BIOGRID-ORGANISM/BIOGRID-ORGANISM-Homo_sapiens-3.5.178.tab2.txt";
//		String bioGridFile = "BIOGRID-ALL-3.5.178.tab2.txt";
		List<String> interactors = new ArrayList<>();
		Map<String, String> entrezGene2UniProt = new HashMap<>();
		Set<String> identifiersToMapToUniprot = new HashSet<>();
		// TAXONOMY ID is 9606
		final String humanOranismTaxonID = "9606";
		try(FileReader reader = new FileReader(bioGridFile);
				FileWriter writer = new FileWriter(BIOGRID_PPIS_FILE);
				FileWriter unmappedWriter = new FileWriter("unmapped_from_biogrid.txt"))
		{
			CSVFormat format = CSVFormat.DEFAULT
								.withDelimiter('\t')
								.withFirstRecordAsHeader();
			Set<String> unmappedIdentifiers = new HashSet<>();
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				List<CSVRecord> records = parser.getRecords();
				int totalNumRecords = records.size();
				int unMappedCount = 0;
				int recordsWithNonHumanInteractor = 0;
				for (CSVRecord record : records)
				{
					String organismA = record.get("Organism Interactor A");
					String organismB = record.get("Organism Interactor B");

					if (!organismA.trim().equals(humanOranismTaxonID) || !organismB.trim().equals(humanOranismTaxonID))
					{
//						System.out.println("Non-human organism ID: " + interactorA + " , " + interactorB);
//						System.exit(-1);
						recordsWithNonHumanInteractor++;
					}

					if (organismA != null && organismB != null
						&& organismA.trim().equals(humanOranismTaxonID)
						&& organismB.trim().equals(humanOranismTaxonID))
					{
						String score = record.get("Score");
//						if (score != null && !score.equals("-") && Double.parseDouble(score) > 0.75d)
						{
							String entrezGeneA = record.get("Entrez Gene Interactor A");
							String entrezGeneB = record.get("Entrez Gene Interactor B");
							identifiersToMapToUniprot.add(entrezGeneA);
							identifiersToMapToUniprot.add(entrezGeneB);
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
//							if (unMappedCount % 5000 == 0 || record.getRecordNumber() >= totalNumRecords)
//							{
//								unMappedCount = 0;
//								getMappingsFromUniProt(entrezGene2UniProt, UniprotDB.Entrez_Gene);
//							}
							interactors.add(entrezGeneA + "|" + entrezGeneB);
						}

					}
				}
				System.out.println(recordsWithNonHumanInteractor + " interactors pairs contained a non-human interactor.");
			}
			entrezGene2UniProt = getMappingsFromUniProt(identifiersToMapToUniprot, UniprotDB.Entrez_Gene);
			int unmappedCount = 0;
			Set<String> interactorSet = new HashSet<>();
			int mappedCount = 0;
			// Now, we go through the interactors, map them to UniProt and write them to a file.
			for (String interactorPair : interactors)
			{
				String[] parts = interactorPair.split("\\|");
				String interactorA = parts[0];
				String interactorB = parts[1];

				String mappedInteractorAId = entrezGene2UniProt.get(interactorA);
				String mappedInteractorBId = entrezGene2UniProt.get(interactorB);


				if (mappedInteractorAId!=null && mappedInteractorBId!=null &&
						!mappedInteractorAId.equals(EMPTY_TOKEN) && !mappedInteractorBId.equals(EMPTY_TOKEN))
				{
					mappedCount++;
					interactorSet.add(putProteinsInOrder(mappedInteractorAId, mappedInteractorBId));
//					writer.write( putProteinsInOrder(mappedInteractorAId, mappedInteractorBId)+"\n" );
				}
				else
				{
					unmappedCount ++;
					if (mappedInteractorAId==null || mappedInteractorAId.equals(EMPTY_TOKEN))
					{
//						unmappedWriter.write(interactorA + "\n");
						unmappedIdentifiers.add(interactorA);
//						System.out.println("Could not map interactor " + interactorA + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
					}
					if (mappedInteractorBId==null || mappedInteractorBId.equals(EMPTY_TOKEN))
					{
//						unmappedWriter.write(interactorB + "\n");
						unmappedIdentifiers.add(interactorB);
//						System.out.println("Could not map interactor " + interactorB + " to uniprot, interactions that include it will NOT be in the overlap analysis.");
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
			System.out.println("BioGrid Interactions that were mapped: "+mappedCount+"; Unmapped interactors: "+unmappedCount);

		}

	}
}
