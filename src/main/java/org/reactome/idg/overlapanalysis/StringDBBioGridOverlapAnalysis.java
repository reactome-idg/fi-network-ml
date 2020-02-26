package org.reactome.idg.overlapanalysis;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.util.StringDBUtil;

public class StringDBBioGridOverlapAnalysis
{
	private static Logger logger = LogManager.getLogger();
	private static String pathToData = "src/main/resources/data/";
	private static String pathToOutput = "output/overlaps/";

	public static void main(String... args)
	{
		logger.info("Calculating StringDB/Biogrid PPI overlap...");
		String stringDBProteinActionsFile = pathToData + "9606.protein.actions.v11.0.txt";
		String stringDBProteinLinksFile = pathToData + "9606.protein.links.full.v11.0.txt";
		String stringToEntrezGeneMapping = pathToData + "all_organisms.entrez_2_string.2018.tsv";
		String bioGridPPIFile = pathToData + "BIOGRID-ORGANISM/BIOGRID-ORGANISM-Homo_sapiens-3.5.181.tab2.txt";
		String humanOranismTaxonID = "9606";
		// Step 1: get PPIs from String
		// Step 2: get PPIs from BioGrid
		// Step 3: map BioGrid PPIs to StringDB
		// Step 4: calculate overlap
		Set<String> interactionsWithExperiments = StringDBUtil.getPPIsWithExperiments(0, stringDBProteinLinksFile);
		Set<String> bindingInteractionsWithExperiments = StringDBUtil.getBindingPPIs(stringDBProteinActionsFile, interactionsWithExperiments);
		logger.info("Number of PPIs from StringDB: {}", bindingInteractionsWithExperiments.size());
		// Get BioGrid data
		Set<String> ppisFromBioGrid = new HashSet<>();
		try(FileReader reader = new FileReader(bioGridPPIFile))
		{
			CSVFormat format = CSVFormat.DEFAULT
					.withDelimiter('\t')
					.withFirstRecordAsHeader();
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


							ppisFromBioGrid.add(StringDBUtil.putProteinsInOrder(entrezGeneA, entrezGeneB));
						}

					}
				}
			}
		}
		catch (IOException e)
		{
			logger.error("File error!", e);
		}
		logger.info("Number of PPIs from BioGrid: {}", ppisFromBioGrid.size());
		Map<String, String> entrezGeneToString = new HashMap<>();
		// Map the BioGrid data to StringDB using StringDB's EntrezGene<->StringDB mapping file.
		try(FileReader reader = new FileReader(stringToEntrezGeneMapping))
		{
			CSVFormat format = CSVFormat.DEFAULT.withCommentMarker('#').withDelimiter('\t');
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				parser.forEach( record -> {
					if (record.get(0).equals(humanOranismTaxonID))
					{
						String entrezGeneID = record.get(1);
						String stringDBID = record.get(2);
						if (entrezGeneID != null && stringDBID != null)
						{
							entrezGeneToString.put(entrezGeneID, stringDBID);
						}
					}
				});
			}
		}
		catch (IOException e)
		{
			logger.error("File error!", e);
		}
		logger.info("Number of EntrezGene -> StringDB identifier mappings: {}", entrezGeneToString.size());


		// Now need to update the BioGrid PPIs
		Set<String> mappedPPIsFromBioGrid = new HashSet<>();
		ppisFromBioGrid.forEach(ppi -> {
			String[] parts = ppi.split("\\t");
			String mapped1 = entrezGeneToString.get(parts[0]);
			String mapped2 = entrezGeneToString.get(parts[1]);

			if ((mapped1 != null && mapped2 != null) &&
				(!mapped1.equals(mapped2)))
			{
				mappedPPIsFromBioGrid.add(StringDBUtil.putProteinsInOrder(mapped1, mapped2));
			}
		});
		logger.info("Number of MAPPED (to StringDB) BioGrid PPIs: {}", mappedPPIsFromBioGrid.size());

		// Now calculate overlap.
		Set<String> overlap = new HashSet<>(bindingInteractionsWithExperiments);
		overlap.retainAll(mappedPPIsFromBioGrid);
		logger.info("Size of StringDB/BioGrid PPI Overlap: {}", overlap.size());
		try(FileWriter overlapWriter = new FileWriter(pathToOutput + "StringDB-BioGrid-PPIoverlap.tsv");
			FileWriter stringDBRemainderWriter = new FileWriter(pathToOutput + "StringDB-only-PPIs.tsv");
			FileWriter bioGridRemainderWriter = new FileWriter(pathToOutput + "BioGrid-only-PPIs.tsv"))
		{
			for (String ppi : overlap.stream().sorted().collect(Collectors.toList()))
			{
				overlapWriter.write(ppi+"\n");
			}
			// Now write the StringDB PPIs that are not in overlap
			for (String ppi : bindingInteractionsWithExperiments.stream().sorted().collect(Collectors.toList()))
			{
				if (!overlap.contains(ppi))
				{
					stringDBRemainderWriter.write(ppi + "\n");
				}
			}
			// Now write the BioGrid PPIs that are not in overlap
			for (String ppi : mappedPPIsFromBioGrid.stream().sorted().collect(Collectors.toList()))
			{
				if (!overlap.contains(ppi))
				{
					bioGridRemainderWriter.write(ppi + "\n");
				}
			}
		}
		catch (IOException e)
		{
			logger.error("File error!", e);
		}
	}

}
