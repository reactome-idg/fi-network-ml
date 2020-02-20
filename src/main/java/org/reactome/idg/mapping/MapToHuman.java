package org.reactome.idg.mapping;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.reactome.idg.util.Species;

public class MapToHuman
{
	private static final Logger logger = LogManager.getLogger(MapToHuman.class);
	private static final String PATH_TO_DATA_FILES = Paths.get("src","main","resources","data").toString() + File.separator;
	private Set<String> unmappedProteins = new TreeSet<>();


	private String outputPath = "";

	private int stringDBSpeciesCode;

	public static void main(String[] args)
	{
		MapToHuman mapper = new MapToHuman();
		mapper.generateOtherSpeciesToHumanPPIs(Species.YEAST);
		mapper.generateOtherSpeciesToHumanPPIs(Species.SCHPO);
	}

	enum ProteinIdentifierType
	{
		STRINGDB, UNIPROT_GENE, UNIPROT_ACCESSION, SGD
	}

	class Protein implements Comparable<Protein>
	{
		private String identifierValue;
		ProteinIdentifierType identifierType;
		public Protein(String value, ProteinIdentifierType type)
		{
			this.identifierValue = value;
			this.identifierType = type;
		}

		@Override
		public String toString()
		{
			// Don't want identifierType in string representation of protein. Maybe have a "detailedToString" that could include it.
			return this.identifierValue;
		}

		public String toDetailedString()
		{
			return this.identifierType.toString() + ":" + this.identifierValue;
		}

		@Override
		public int hashCode()
		{
			return this.toDetailedString().hashCode();
		}

		@Override
		public boolean equals(Object other)
		{
			if (other == null)
			{
				return false;
			}
			if (other instanceof Protein)
			{
				return (this.identifierValue == ((Protein)other).identifierValue && this.identifierType == ((Protein)other).identifierType);
			}
			return other.equals(this);
		}

		@Override
		public int compareTo(Protein protein2)
		{
			return this.identifierValue.compareTo(protein2.identifierValue);
		}

		public String getIdentifierValue()
		{
			return identifierValue;
		}

		public ProteinIdentifierType getIdentifierType()
		{
			return identifierType;
		}
	}

	class Ppi implements Comparable<Ppi>
	{
		private Protein protein1;
		private Protein protein2;

		public Ppi(Protein protein1, Protein protein2)
		{
			this.protein1 = protein1;
			this.protein2 = protein2;
		}

		@Override public String toString()
		{
			return this.protein1.compareTo(this.protein2) < 0
					? this.protein1 + "\t" + this.protein2
					: this.protein2 + "\t" + this.protein1;
		}

		@Override
		public int hashCode()
		{
			return this.toString().hashCode();
		}

		@Override
		public boolean equals(Object other)
		{
			if (other == null)
			{
				return false;
			}
			return this.toString().equals(other.toString());
		}

		public Protein getProtein1()
		{
			return protein1;
		}

		public Protein getProtein2()
		{
			return protein2;
		}

		@Override
		public int compareTo(Ppi other)
		{
			return this.toString().compareTo(other.toString());
		}
	}

	private Set<Ppi> getBindingPPIs(String stringDBActionsFile)
	{
		Set<Ppi> ppis = new HashSet<>();
		String bindingPPIs = outputPath + this.stringDBSpeciesCode + "_binding_PPIs.tsv";
		try(FileReader reader = new FileReader(stringDBActionsFile);
			CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
				.withDelimiter('\t')
				.withFirstRecordAsHeader());
			FileWriter writer = new FileWriter(bindingPPIs);)
		{
			List<CSVRecord> records = parser.getRecords();
			for (CSVRecord record : records)
			{
				// Only extract records where mode==binding
				if (record.get("mode").equals("binding"))
				{
					Protein protein1 = new Protein(record.get("item_id_a"), ProteinIdentifierType.STRINGDB);
					Protein protein2 = new Protein(record.get("item_id_b"), ProteinIdentifierType.STRINGDB);
					Ppi ppi = new Ppi(protein1, protein2);
					ppis.add(ppi);
					writer.write(ppi.toString()+"\n");
				}
			}
		}
		catch (IOException e)
		{
			logger.error("File Error!", e);
		}
		return ppis;
	}

	private Set<Ppi> getPPIsWithExperiments(String stringDBProteinLinksFile)
	{
		Set<Ppi> ppis = new HashSet<>();
		String ppisWithExperimentsScoreFile = outputPath + this.stringDBSpeciesCode + "_PPIs_with_experiments.tsv";
		// First we have to filter for interactions in protein.links.full where experiment > 0
		try (FileReader reader = new FileReader(stringDBProteinLinksFile);
				CSVParser parser = new CSVParser(reader, CSVFormat.DEFAULT
										.withDelimiter(' ')
										.withFirstRecordAsHeader());
				FileWriter writer = new FileWriter(ppisWithExperimentsScoreFile);)
		{
			List<CSVRecord> records = parser.getRecords();
			logger.info("Reading StringDB \"links\" file; {} records will be parsed.", records.size());
			for (CSVRecord  record : records)
			{
				int experiments = Integer.parseInt(record.get("experiments"));
//				int dbScore = Integer.parseInt(record.get("database"));
				// Only extract data for records where "experiments > 0"
				if (experiments > 0 /*|| dbScore > 0*/)
				{
					Protein protein1 = new Protein (record.get("protein1") , ProteinIdentifierType.STRINGDB);
					Protein protein2 = new Protein (record.get("protein2"), ProteinIdentifierType.STRINGDB);
					Ppi ppi = new Ppi(protein1, protein2);
					ppis.add(ppi);
					writer.write(ppi.toString()+"\n");
				}
			}
		}
		catch (IOException e)
		{
			logger.error("File Error!", e);
		}
		return ppis;
	}

	private Set<Ppi> getBindingPPIsWithExperiments(String stringDBActionsFile, String stringDBProteinLinksFile)
	{
		Set<Ppi> ppisWithExperiments = getPPIsWithExperiments(stringDBProteinLinksFile);
		logger.info("{} PPIs with experiments > 0", ppisWithExperiments.size());
		Set<Ppi> ppisBinding = getBindingPPIs(stringDBActionsFile);
		logger.info("{} PPIs with mode==binding", ppisBinding.size());
		Set<Ppi> ppisBindingAndExperiments = ppisBinding;
		// this will calculate the intersection of "PPIs with experiments > 0" and "Binding PPIs"
		ppisBindingAndExperiments.retainAll(ppisWithExperiments);
		logger.info("{} PPIs with experiments > 0 AND mode==binding", ppisBindingAndExperiments.size());
		return ppisBindingAndExperiments;
	}

	private String[] extractIdentifiersFromPantherLine(String[] parts, String species1)
	{
		String[] humanParts;
		String[] otherSpeciesParts;
		if (species1.equals(Species.HUMAN.toString()))
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
		String humanIdentifier = this.extractUniprotIdentifier(humanParts);
		String otherSpeciesIdentifier = this.extractUniprotIdentifier(otherSpeciesParts);
		// Return the human and non-human identifier from the Panther file.
		return new String[]{ humanIdentifier, otherSpeciesIdentifier } ;
	}

	private Map<String, Set<String>> getOrthologMappedProteins(String pathToOrthologMappingFile, boolean allowBidirectionalMapping, String speciesName)
	{
		Map<String, Set<String>> otherSpeciesMappedToHuman = new HashMap<>();
		try(BufferedReader br = new BufferedReader(new FileReader(pathToOrthologMappingFile)))
		{
			String line = br.readLine();

			while (line != null)
			{
				String[] parts = line.split("\\t");
				// Only the first two parts are useful, and only if they contain HUMAN or YEAST (or whatever speciesName is).
				String species1 = (parts[0].split("\\|"))[0];
				String species2 = (parts[1].split("\\|"))[0];
				// TODO: Generate new mappings based on gene common families. Maybe store these in a separate mapping structure so
				// that they can be marked as originating from a *generated* mapping rather than a *known* mapping.
				String geneFamily = parts[4];
				boolean lineIsUseful = !species1.equals(species2) && ((allowBidirectionalMapping && species1.equals(Species.HUMAN.toString()) && species2.equals(speciesName))
																	|| (species1.equals(speciesName) && species2.equals(Species.HUMAN.toString())));
				if (lineIsUseful)
				{
					String[] extractedIdentifiers = this.extractIdentifiersFromPantherLine(parts, species1);
					String humanIdentifier = extractedIdentifiers[0];
					String otherSpeciesIdentifier = extractedIdentifiers[1];
					Set<String> identifiers;
					if (otherSpeciesMappedToHuman.containsKey(otherSpeciesIdentifier))
					{
						identifiers = otherSpeciesMappedToHuman.get(otherSpeciesIdentifier);
					}
					else
					{
						identifiers = new TreeSet<>();
					}
					identifiers.add(humanIdentifier);
					otherSpeciesMappedToHuman.put(otherSpeciesIdentifier, identifiers);
				}
				line = br.readLine();
			}
		}
		catch (IOException e)
		{
			logger.error("File error!", e);
		}
		return otherSpeciesMappedToHuman;
	}

	private Map<String, Set<String>> getStringDBToUniProtMappings(String pathToMappings)
	{
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
				if (Integer.parseInt(parts[0]) == this.stringDBSpeciesCode)
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

	private void generateOtherSpeciesToHumanPPIs(Species species)
	{
		this.stringDBSpeciesCode = species.getSpeciesCode();
		this.outputPath = "output/"+species+"_results/";

		logger.info("Generating PPIs for {}", species);

		String stringDBProteinActionsFile = PATH_TO_DATA_FILES + this.stringDBSpeciesCode + ".protein.actions.v11.0.txt";
		String stringDBProteinLinksFile = PATH_TO_DATA_FILES + this.stringDBSpeciesCode + ".protein.links.full.v11.0.txt";

		Set<Ppi> ppisFromStringDB = this.getBindingPPIsWithExperiments(stringDBProteinActionsFile, stringDBProteinLinksFile);
		String pathToOrthologMappingFile = PATH_TO_DATA_FILES + "Orthologs_HCOP";
		Map<String, Set<String>> orthologProteins = getOrthologMappedProteins(pathToOrthologMappingFile, true, species.toString());
		logger.info("{} ortholog proteins", orthologProteins.size());
		Map<String, Set<String>> otherSpeciesStringDBToUniProt = getStringDBToUniProtMappings(PATH_TO_DATA_FILES + "all_organisms.uniprot_2_string.2018.tsv");
		logger.info("{} StringDB-to-UniProt identifiers", otherSpeciesStringDBToUniProt.size());
		String pathToOutputPPIFile = this.outputPath + species + "_MAPPED_PPIS.tsv";
		int selfInteractions = 0;
		int unmappedPPIs = 0;
		Set<String> ppiLines = new TreeSet<>();

		Map<Ppi, Set<String>> orthologMap = new HashMap<>(ppisFromStringDB.size());
		// If a PPI from StringDB (with experiments > 0 && mode==binding) has both proteins in the mapping, then it's OK! ...and should be written to output file.
		try(FileWriter writer = new FileWriter(pathToOutputPPIFile);
			FileWriter unMappedIdentifiers = new FileWriter(this.outputPath + species + "_unmapped_proteins.txt"))
		{
			for (Ppi ppi : ppisFromStringDB.stream().sorted().collect(Collectors.toList()))
			{
				Protein p1 = ppi.getProtein1();
				Protein p2 = ppi.getProtein2();
				Set<String> p1Orthologs = null;
				Set<String> p2Orthologs = null;
				String p1AsUniProt = null;
				String p2AsUniProt = null;
				if (p1.getIdentifierType().equals(ProteinIdentifierType.STRINGDB))
				{
					p1AsUniProt = getProteinAsUniProt(otherSpeciesStringDBToUniProt, p1);
				}
				if (p2.getIdentifierType().equals(ProteinIdentifierType.STRINGDB))
				{
					p2AsUniProt = getProteinAsUniProt(otherSpeciesStringDBToUniProt, p2);
				}

				if (p1AsUniProt != null && p2AsUniProt != null)
				{
					p1Orthologs = orthologProteins.get(p1AsUniProt);
					p2Orthologs = orthologProteins.get(p2AsUniProt);
					// Only continue if this is NOT a self-interaction.
					if (!p1AsUniProt.equals(p2AsUniProt))
					{
						// Output will include the proteins, and their mapping.
						StringBuilder sb = new StringBuilder();
						sb.append(p1AsUniProt).append('\t').append(p2AsUniProt).append("\t")
							.append(p1AsUniProt).append(" map from: ").append(ppi.getProtein1()).append(", ")
							.append(p2AsUniProt).append(" map from: ").append(ppi.getProtein2()).append(System.lineSeparator());
						ppiLines.add(sb.toString());
						// Now... we do the orthologs... (if any exist)
						if (p1Orthologs != null && p2Orthologs != null)
						{
							for (String p1Ortholog : p1Orthologs)
							{
								for (String p2Ortholog : p2Orthologs)
								{
									if (!p1Ortholog.equals(p2Ortholog))
									{
										sb = new StringBuilder();
										// Create a new Ortholog PPI
										Ppi orthologPPI = new Ppi(new Protein(p1Ortholog, ProteinIdentifierType.UNIPROT_ACCESSION), new Protein(p2Ortholog, ProteinIdentifierType.UNIPROT_ACCESSION));
										// Create a PPI based on the UniProt identifiers.
										Ppi ppiAsUniProt = new Ppi(new Protein(p1AsUniProt, ProteinIdentifierType.UNIPROT_ACCESSION), new Protein(p2AsUniProt, ProteinIdentifierType.UNIPROT_ACCESSION));
										// build up the list of reasons that this PPI is in the output: for orthologPPI, what are its proteins orthologs of? That goes in the output
										// to make it easier to trace why the PPI is there.
										Set<String> reasons = orthologMap.computeIfAbsent(orthologPPI, x -> new TreeSet<>());
										sb.append(orthologPPI.getProtein1()).append(" ortholog of ").append(ppiAsUniProt.getProtein1());
										reasons.add(sb.toString());
										// Now a new reason for Protein 2 (each protein should be a separate reason).
										sb = new StringBuilder();
										sb.append(orthologPPI.getProtein2()).append(" ortholog of ").append(ppiAsUniProt.getProtein2());
										reasons.add(sb.toString());
										// Use the new Ortholog PPI as the *key* in orthologMap - that way we can build up a list of "reasons" for this PPI.
										orthologMap.put(orthologPPI, reasons);
									}
									else
									{
										selfInteractions++;
									}
								}
							}
						}
					}
					else
					{
						selfInteractions++;
					}
				}
				else
				{
					// If one of the proteins-as-UniProt values is NULL, it means there was a mapping problem.
					// NOTE: this counter will be incremented for each PPI where a mapping problem occurs, but
					// a single unmapped protein could be in multiple PPIs.
					unmappedPPIs++;
				}
			}
			// Write the PPIs to output.
			logger.info("{} PPIs in output", ppiLines.size());
			for (String line : ppiLines)
			{
				writer.write(line);
			}
			// Write the PPIs from ortholog mapping to output.
			logger.info("{} PPIs (from Orthologs) in output", orthologMap.size());
			for (Entry<Ppi, Set<String>> ortholog : orthologMap.entrySet()/*.stream().sorted( (e1,e2) -> e1.getKey().compareTo(e2.getKey()) ).collect(Collectors.toList())*/)
			{
				writer.write(ortholog.getKey().toString() + "\t(");
				for (String reason : ortholog.getValue())
				{
					writer.write(reason + "; ");
				}
				writer.write(")" + System.lineSeparator());
			}
			logger.info("{} distinct unmapped proteins.", this.unmappedProteins.size());
			for (String unmappedProtein : this.unmappedProteins)
			{
				unMappedIdentifiers.write(unmappedProtein + System.lineSeparator());
			}
		}
		catch (IOException e)
		{
			logger.error("File error!", e);
		}
		logger.info("{} PPIs had mapping problems.", unmappedPPIs);
		logger.info("{} self-interactions were ommitted.", selfInteractions);
	}

	private String getProteinAsUniProt(Map<String, Set<String>> otherSpeciesStringDBToUniProt, Protein protein)
	{
		String proteinAsUniProt = null;
		if (otherSpeciesStringDBToUniProt.containsKey(protein.getIdentifierValue()))
		{
			proteinAsUniProt = otherSpeciesStringDBToUniProt.get(protein.getIdentifierValue()).stream().findFirst().get();
		}
		else
		{
			this.unmappedProteins.add(protein.toString());
		}
		return proteinAsUniProt;
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
}
