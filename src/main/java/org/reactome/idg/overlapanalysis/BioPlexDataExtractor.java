package org.reactome.idg.overlapanalysis;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

public class BioPlexDataExtractor extends DataExtractor
{
	private static final String BIOPLEX_PPIS_FILE = "bioplex-ppis.txt";
	public static void extractFromBioPlexFile() throws FileNotFoundException, IOException
	{
		System.out.println("Extracting BioPlex Data...");
		String bioPlexV3File1 = "BioPlex_2.3_interactionList.tsv";
		String bioPlexV3File2 = "BioPlex_unpublishedInteractions_May_2019.tsv";
		Set<String> interactors = new HashSet<>();
		try(FileReader reader1 = new FileReader(bioPlexV3File1);
			FileReader reader2 = new FileReader(bioPlexV3File2);
				FileWriter writer = new FileWriter(BIOPLEX_PPIS_FILE))
		{
			CSVFormat format = CSVFormat.DEFAULT
								.withDelimiter('\t')
								.withFirstRecordAsHeader();

			Function<FileReader, Integer> processBioPlexFile = (reader) -> {
				int numInFile = 0;

				try(CSVParser parser = new CSVParser(reader, format);)
				{
					for (CSVRecord record : parser.getRecords())
					{
						if (Double.parseDouble(record.get("p(Interaction)")) > 0.75d
								&& (!record.get("UniprotA").equals("UNKNOWN") && !record.get("UniprotB").equals("UNKNOWN")))
						{
							String protein1 = record.get("UniprotA");
							String protein2 = record.get("UniprotB");
							numInFile++;
							interactors.add( putProteinsInOrder(protein1, protein2) );
						}
					}
					System.out.println("Number of Bioplex interactions: "+numInFile + "; out of a total of "+parser.getRecordNumber());
				}
				catch (IOException e)
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				return numInFile;
			};
			processBioPlexFile.apply(reader1);
			processBioPlexFile.apply(reader2);
			for (String interactor : interactors)
			{
				writer.write(interactor + "\n");
			}
		}
	}
}
