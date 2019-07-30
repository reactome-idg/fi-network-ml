package org.reactome.idg.model;

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
	@GeneratedValue(strategy = GenerationType.AUTO)
	private Integer id;
	
	private String name;
	private String url;
	
	// Should these be free-form text fields or should they be foreign keys to proper classes?
	private String category;
	private String subcategory;

    public Provenance()
    {
    }
    
	public Provenance(String name, String URL, String category, String subcategory)
	{
		this.setName(name);
		this.setUrl(URL);
		this.setCategory(category);
		this.setSubcategory(subcategory);
	}

	public String getName()
	{
		return name;
	}
	
	public void setName(String name)
	{
		this.name = name;
	}
	public String getUrl()
	{
		return url;
	}
	
	public void setUrl(String url)
	{
		this.url = url;
	}
	
	public String getCategory()
	{
		return category;
	}
	
	public void setCategory(String category)
	{
		this.category = category;
	}
	
	public String getSubcategory()
	{
		return subcategory;
	}
	
	public void setSubcategory(String subcategory)
	{
		this.subcategory = subcategory;
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
	 * String representation of a Provenance. Includes the fields: ID, Name, URL. Other fields ommitted for brevity.
	 */
	@Override
	public String toString()
	{
		return "ID: " + this.getId() + "; Name: " + this.getName() + "; URL: "+this.getUrl();
	}
	
	//TODO: Add: species, tissue-type, other fields...
}
