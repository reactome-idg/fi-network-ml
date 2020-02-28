package org.reactome.idg.model;


public class ProteinProteinInteraction implements Comparable<ProteinProteinInteraction>
{
	private Protein protein1;
	private Protein protein2;

	public ProteinProteinInteraction(Protein protein1, Protein protein2)
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
	public int compareTo(ProteinProteinInteraction other)
	{
		return this.toString().compareTo(other.toString());
	}
}
