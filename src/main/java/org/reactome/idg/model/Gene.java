package org.reactome.idg.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.Table;

@Entity
@Table(name = "gene", indexes = { @Index(columnList = "symbol", unique = true, name = "idx_unq_symbol") })
public class Gene implements Comparable<Gene>
{
	@Column(name = "symbol", length = 16, nullable = false)
	private String symbol;
	
	@Id
	@GeneratedValue(strategy = GenerationType.AUTO)
	private Integer id;
	
	public Gene() {
	    
	}
	
	public String getSymbol()
	{
		return symbol;
	}
	
	public void setSymbol(String symbol)
	{
		this.symbol = symbol;
	}
	
	public Integer getId()
	{
		return id;
	}
	
	public void setId(Integer id)
	{
		this.id = id;
	}
	
	/**
	 * Genes are compared by their symbol, not their ID.
	 * @param otherGene - the other Gene to compare to.
	 */
	@Override
	public int compareTo(Gene otherGene)
	{
		return this.symbol.compareTo(otherGene.getSymbol());
	}
	
}
