package org.reactome.idg.model;

public class Protein implements Comparable<Protein>
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