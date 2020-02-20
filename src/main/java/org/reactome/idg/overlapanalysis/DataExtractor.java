package org.reactome.idg.overlapanalysis;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.idg.util.UniprotFileRetriever;
import org.reactome.idg.util.UniprotFileRetriever.UniprotDB;

abstract class DataExtractor
{
	protected static final String EMPTY_TOKEN = "<EMPTY>";
	protected static String putProteinsInOrder(String p1, String p2)
	{
		// Use String's compareTo to put the proteins in order, and put a tab between them.
		return p1.compareTo(p2) < 0
				? p1 + "\t" + p2
				: p2 + "\t" + p1;
	}

	protected static /*int */ void getMappingsFromUniProt(Map<String, String> mappings, UniprotDB targetDB) throws URISyntaxException
	{
//		int unMappedCount;
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
//		unMappedCount = 0;
//		return unMappedCount;
	}

	protected Map<String, String> getMappingsFromUniProt(Set<String> identifiersToMapToUniprot, UniprotDB targetDB) throws URISyntaxException
	{
		Map<String, String> ensemblToUniprotMap = new HashMap<>();
		int i = 0, j = 0;
		int k = 2000;
		String inputString = "";
		for (String ensemblIdentifier : identifiersToMapToUniprot)
		{
			i++;
			j++;
			inputString += ensemblIdentifier + " ";
//			if (ensemblIdentifier.equals("ENSP00000243344"))
//			{
//				System.out.println("ensemblIdentifier == ENSP00000243344");
//			}

			if (i % k == 0 || identifiersToMapToUniprot.size() == j)
			{

				UniprotFileRetriever retriever = new UniprotFileRetriever();
				retriever.setMapToDbEnum(UniprotDB.UniProtAcc);
				retriever.setMapFromDbEnum(targetDB);
				retriever.setUri(new URI("https://www.uniprot.org/uploadlists/"));
				retriever.setDataInputStream(new BufferedInputStream(new ByteArrayInputStream(inputString.getBytes())));
//				if (inputString.contains("ENSP00000243344"))
//				{
//					System.out.println("input string for UniProt mapping service contains ENSP00000243344");
//				}
				List<String> dataLines = retriever.downloadAndReturnDataLines();
				System.out.println(dataLines.size() + " " + targetDB.toString() + "-to-Uniprot mappings retrieved.");
				// Just in case...
				if (dataLines!=null)
				{
					for (String line : dataLines)
					{
//						if (line.contains("ENSP00000243344"))
//						{
//							System.out.println(line);
//						}

						if (!line.equals("From	To") && !line.equals("not mapped"))
						{
							String[] parts = line.split("\\t");
							String ensemblIDentifier = parts[0];
							String uniProtIdentifier = parts[1];

							ensemblToUniprotMap.put(ensemblIDentifier, uniProtIdentifier);
						}
					}
				}

				i = 0;
				inputString = "";
			}
		}

		return ensemblToUniprotMap;
	}
}
