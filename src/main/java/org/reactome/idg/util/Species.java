package org.reactome.idg.util;

/**
 * Represents species. Each species has a corresponding numeric code that can be obtained with the getSpeciesCode method.
 * @author sshorser
 *
 */
public enum Species
{
	YEAST(4932, "Saccharomyces cerevisiae"), HUMAN(9606, "Homo sapiens"), SCHPO(4896, "Schizosaccharomyces pombe");

	private int speciesCode;
	private String fullName;

	Species(int code, String name)
	{
		this.speciesCode = code;
		this.fullName = name;
	}

	public int getSpeciesCode()
	{
		return this.speciesCode;
	}

	public String getFullName()
	{
		return this.fullName;
	}
}