package org.reactome.idg.model;

import javax.persistence.Column;
import javax.persistence.ConstraintMode;
import javax.persistence.Entity;
import javax.persistence.ForeignKey;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

//TODO: Introduce GENE class, and it is the gene IDs normalized so that this table has foreign keys to genes, rather than strings that ARE genes.
@Entity
@Table(name = "gene_pair_correlation", indexes = { @Index(columnList = "gene_1_id,gene_2_id", unique = false, name = "idx_gene_pair"),
													@Index(columnList = "gene_1_id,gene_2_id,provenance_id", unique = true, name = "idx_gene_pair_provenance")})
public class GenePairCorrelation
{
	public GenePairCorrelation()
	{
		// Hibernate requires a deafult no-arg constructor.
	}
	
	public GenePairCorrelation(Gene gene1, Gene gene2, double correlationValue, Provenance provenance)
	{
		this.gene1 = gene1;
		this.gene2 = gene2;
		this.correlationValue = correlationValue;
		this.provenance = provenance;
	}
	
	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id;

	@JoinColumn(name = "gene_1_id", foreignKey = @ForeignKey(value = ConstraintMode.CONSTRAINT), nullable = false)
	@ManyToOne
	private Gene gene1;
	
	@JoinColumn(name = "gene_2_id", foreignKey = @ForeignKey(value = ConstraintMode.CONSTRAINT), nullable = false)
	@ManyToOne
	private Gene gene2;
	
	@Column(name = "correlation_value", nullable = false, columnDefinition = "DOUBLE PRECISION(5,4)")
	private double correlationValue;
	
	@JoinColumn(name = "provenance_id", foreignKey = @ForeignKey(value = ConstraintMode.CONSTRAINT), nullable = false)
	@ManyToOne
	private Provenance provenance;

	public GenePairCorrelation(Long id)
	{
		this.id = id;
	}
 
	public Long getId()
	{
		return id;
	}

	public Gene getGene1()
	{
		return gene1;
	}

	public void setGene1(Gene gene1)
	{
		this.gene1 = gene1;
	}

	public Gene getGene2()
	{
		return gene2;
	}

	public void setGene2(Gene gene2)
	{
		this.gene2 = gene2;
	}

	public double getCorrelationValue()
	{
		return correlationValue;
	}

	public void setCorrelationValue(double correlationValue)
	{
		this.correlationValue = correlationValue;
	}

	public Provenance getProvenance()
	{
		return provenance;
	}

	public void setProvenance(Provenance provenance)
	{
		this.provenance = provenance;
	}
	
}