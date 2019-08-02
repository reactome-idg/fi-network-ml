package org.reactome.idg.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Table;

/**
 * Represents the provenance of a gene-pair correlation.
 * Attributes of Provenance are name, URL, category, subcategory.
 * @author sshorser
 *
 */
@Entity
@Table(name = "provenance" /* , indexes = { @Index(columnList = "name"), @Index(columnList = "url"), @Index(columnList = "category") } */ )
public class Provenance  {
	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	@Column(name = "id", nullable = false, columnDefinition = "SMALLINT")
	private Integer id;
	// E.g. dataset in the harmonizome
	private String name;
	// Most cases url should be unique. However, it may not in cases
	// like ARCHS4, since we will do some post-process (e.g. tissue-specific calculation)
	private String url;
	// e.g. the category in the Harmonizome data
	private String category;
	// e.g. cell line or cell type. Free text for now.
	private String biologicalEntity;
	// e.g. Harmonizome or ARCHS4. Whatever the name used by the original database
	private String source;

	public Provenance() {
	}

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public Integer getId() {
        return id;
    }

    public void setId(Integer id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }

    public String getCategory() {
        return category;
    }

    public void setCategory(String category) {
        this.category = category;
    }

    public String getBiologicalEntity() {
        return biologicalEntity;
    }

    public void setBiologicalEntity(String biologicalEntity) {
        this.biologicalEntity = biologicalEntity;
    }
	
}
