package org.reactome.idg.util;

import java.io.InputStream;
import java.util.Properties;

import org.apache.log4j.Logger;
import org.gk.persistence.MySQLAdaptor;

public class DatabaseConfig {
	private static final Logger logger = Logger.getLogger(DatabaseConfig.class);
	
	/**
	 * Get the configured MySQLAdaptor database.
	 * @return
	 */
	public static MySQLAdaptor getMySQLDBA() {
		try {
			// Need to load all human genes from a Reactome database
	        InputStream is = DatabaseConfig.class.getClassLoader().getResourceAsStream("db.properties");
	        Properties prop = new Properties();
	        prop.load(is);
	        is.close();
	        MySQLAdaptor dba = new MySQLAdaptor("localhost",
	                                            prop.getProperty("reactome.db.name"), 
	                                            prop.getProperty("mysql.user"), 
	                                            prop.getProperty("mysql.password"));
	        return dba;
		}
		catch(Exception e) {
			logger.error(e.getMessage(), e);
			return null;
		}
	}

}
