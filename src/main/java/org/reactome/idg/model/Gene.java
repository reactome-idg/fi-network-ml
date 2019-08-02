package org.reactome.idg.model;

import java.util.HashSet;
import java.util.Set;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.Table;

@Entity
@Table(name = "gene", indexes = { @Index(columnList = "symbol", unique = true, name = "idx_unq_symbol") })
public class Gene implements Comparable<Gene> {
	// Give it a long enough length in case there are some weird gene names that should not be there.
    @Column(name = "symbol", length = 64, nullable = false)
	private String symbol;
	
	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	@Column(name = "id", nullable = false, columnDefinition = "SMALLINT")
	private Integer id;
	
    @ManyToMany(cascade = CascadeType.MERGE, fetch = FetchType.EAGER)
    @JoinTable(name = "gene_provenances",
               joinColumns = @JoinColumn (name = "gene_id"),
               inverseJoinColumns = @JoinColumn(name = "provenance_id"))
	private Set<Provenance> provenances;
	
	public Gene() {
	}
	
	public Set<Provenance> getProvenances() {
        return provenances;
    }

    public void setProvenances(Set<Provenance> provenances) {
        this.provenances = provenances;
    }
    
    public void addProvenance(Provenance provenance) {
        if (provenances == null)
            provenances = new HashSet<>();
        provenances.add(provenance);
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
