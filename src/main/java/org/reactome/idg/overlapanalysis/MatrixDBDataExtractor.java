package org.reactome.idg.overlapanalysis;

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

public class MatrixDBDataExtractor  extends DataExtractor
{
	private static final String HOMO_SAPIENS = "Homo sapiens";
	private static final String matrixDBFile = "src/main/resources/data/matrixdb_CORE.tab";
	private static final String MATRIXDB_PPIS_FILE = "output/overlaps/matrixdb-ppis.txt";

	public static void extractFromMatrixDBFile() throws FileNotFoundException, IOException
	{
		System.out.println("Extracting MatrixDB Data...");
		final String uniprotPrefix = "uniprotkb:";
		Set<String> interactors = new HashSet<>();
		try(FileReader reader = new FileReader(matrixDBFile);
				FileWriter writer = new FileWriter(MATRIXDB_PPIS_FILE))
		{
			CSVFormat format = CSVFormat.DEFAULT
					.withDelimiter('\t')
					.withQuote('`')
					.withFirstRecordAsHeader();
			try(CSVParser parser = new CSVParser(reader, format);)
			{
				List<CSVRecord> records = parser.getRecords();
				for (CSVRecord record : records)
				{
//					if (!record.get("Confidence value(s)").equals("-"))
					{
						String protein1 = record.get("#ID(s) interactor A");
						String protein2 = record.get("ID(s) interactor B");
						String taxon1 = record.get("Taxid interactor A");
						String taxon2 = record.get("Taxid interactor B");
						// Only include MatrixDB interactions where both identifiers are UniProt.
						if (protein1.contains(uniprotPrefix) && protein2.contains(uniprotPrefix)
							&& taxon1.contains(HOMO_SAPIENS) && taxon2.contains(HOMO_SAPIENS))
						{
//							writer.write( putProteinsInOrder(protein1, protein2)+"\n" );
							interactors.add(putProteinsInOrder(protein1.replace(uniprotPrefix, ""), protein2.replace(uniprotPrefix, "")));
						}
					}
				}
				System.out.println(records.size() +  " records parsed, " + interactors.size() + " interactors selected.");

				for (String interactor : interactors)
				{
					writer.write(interactor + "\n");
				}
			}
		}

	}

}
