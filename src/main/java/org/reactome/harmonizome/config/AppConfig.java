package org.reactome.harmonizome.config;

import static org.hibernate.cfg.AvailableSettings.C3P0_ACQUIRE_INCREMENT;
import static org.hibernate.cfg.AvailableSettings.C3P0_MAX_SIZE;
import static org.hibernate.cfg.AvailableSettings.C3P0_MAX_STATEMENTS;
import static org.hibernate.cfg.AvailableSettings.C3P0_MIN_SIZE;
import static org.hibernate.cfg.AvailableSettings.C3P0_TIMEOUT;
import static org.hibernate.cfg.AvailableSettings.DIALECT;
import static org.hibernate.cfg.AvailableSettings.DRIVER;
import static org.hibernate.cfg.AvailableSettings.HBM2DDL_AUTO;
import static org.hibernate.cfg.AvailableSettings.PASS;
import static org.hibernate.cfg.AvailableSettings.SHOW_SQL;
import static org.hibernate.cfg.AvailableSettings.STATEMENT_BATCH_SIZE;
import static org.hibernate.cfg.AvailableSettings.URL;
import static org.hibernate.cfg.AvailableSettings.USER;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import javax.sql.DataSource;

import org.reactome.harmonizome.CorrelationMatrixLoader;
import org.reactome.harmonizome.HarmonizomeBatch;
import org.reactome.harmonizome.HarmonizomeDataDownloader;
import org.reactome.idg.dao.GeneCorrelationDAOImpl;
import org.reactome.idg.dao.GeneDAOImpl;
import org.reactome.idg.dao.ProvenanceDAOImpl;
import org.springframework.aop.target.CommonsPool2TargetSource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.config.ConfigurableBeanFactory;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.ComponentScans;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.PropertySource;
import org.springframework.context.annotation.Scope;
import org.springframework.core.env.Environment;
import org.springframework.orm.hibernate5.HibernateTransactionManager;
import org.springframework.orm.hibernate5.LocalSessionFactoryBean;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import com.mysql.cj.jdbc.MysqlDataSource;

@Configuration
@PropertySource("file:${pathToProperties}")
@EnableTransactionManagement
@ComponentScans(value = { @ComponentScan("org.reactome.idg.dao"), @ComponentScan("org.reactome.idg.loader"), @ComponentScan("org.reactome.idg") })
public class AppConfig
{
	@Autowired
	private Environment env;

	@Bean(name="dbName")
	public String getDbName()
	{
		return env.getProperty("mysql.db_name");
	}
	
	@Bean(name="chunkSize")
	public int getChunkSize()
	{
		return Integer.parseInt(env.getProperty("chunkSize","1000000"));
	}
	
	@Bean(name="harmonizomeDownloaders")
	public List<HarmonizomeDataDownloader> harminozomeDownloaders() throws IOException, URISyntaxException
	{
		List<HarmonizomeDataDownloader> downloaders = new ArrayList<>();
		List<String> lines = Files.readAllLines(Paths.get(env.getProperty("pathToHarmonizomeFile")));
		String downloadPath = env.getProperty("harmonizomeDownloadPath");
		for (String line : lines)
		{
			String[] parts = line.split("\\t");
			HarmonizomeDataDownloader downloader = new HarmonizomeDataDownloader(new URI(parts[2]), parts[0], parts[1], downloadPath);
			downloaders.add(downloader);
		}
		return downloaders;
	}
	
	@Bean(name = "sessionFactory", autowireCandidate = true)
	public LocalSessionFactoryBean getSessionFactory()
	{
		LocalSessionFactoryBean factoryBean = new LocalSessionFactoryBean();

		Properties props = new Properties();
		// Setting JDBC properties
		props.put(DRIVER, env.getProperty("mysql.driver"));
		props.put(URL, env.getProperty("mysql.url"));
		props.put(USER, env.getProperty("mysql.user"));
		props.put(PASS, env.getProperty("mysql.password"));

		// Setting Hibernate properties
		props.put(SHOW_SQL, env.getProperty("hibernate.show_sql"));
		props.put(HBM2DDL_AUTO, env.getProperty("hibernate.hbm2ddl.auto"));
		// The following cannot work. Have to set in the running property as the
		// following:
		// -Dhibernate.dialect.storage_engine=innodb
//      props.put(STORAGE_ENGINE, env.getProperty("hibernate.dialect.storage_engine"));
		props.put(DIALECT, env.getProperty("hibernate.dialect"));
		
		// Setting C3P0 properties
		props.put(C3P0_MIN_SIZE, env.getProperty("hibernate.c3p0.min_size"));
		props.put(C3P0_MAX_SIZE, env.getProperty("hibernate.c3p0.max_size"));
		props.put(C3P0_ACQUIRE_INCREMENT, env.getProperty("hibernate.c3p0.acquire_increment"));
		props.put(C3P0_TIMEOUT, env.getProperty("hibernate.c3p0.timeout"));
		props.put(C3P0_MAX_STATEMENTS, env.getProperty("hibernate.c3p0.max_statements"));
		props.put(STATEMENT_BATCH_SIZE, env.getProperty("hibernate.jdbc.batch_size"));
		// TODO: Import proper constants for these, also maybe make them configurable in properties file.
		props.put("hibernate.connection.autocommit", false);
		props.put("hibernate.order_inserts", true);
		props.put("hibernate.generate_statistics", false);
//		props.put("hibernate.cache.use_second_level_cache", false);

		// want to see query details.
//		logging.level.org.hibernate.SQL=DEBUG
//		logging.level.org.hibernate.type=TRACE
//		props.put("logging.level.org.hibernate.SQL", "DEBUG");
//		props.put("logging.level.org.hibernate.type", "TRACE");
		
		factoryBean.setHibernateProperties(props);
		factoryBean.setDataSource(this.getDataSource());
		factoryBean.setPackagesToScan("org.reactome.idg.model");
		
		return factoryBean;
	}
	
	@Bean(name = "dataSource")
	public DataSource getDataSource()
	{
		MysqlDataSource dataSource = new MysqlDataSource();
		dataSource.setUrl(env.getProperty("mysql.url"));
		dataSource.setUser(env.getProperty("mysql.user"));
		dataSource.setPassword(env.getProperty("mysql.password"));
		dataSource.setDatabaseName(env.getProperty("mysql.db_name"));
		return dataSource;
	}

	@Bean(name = "transactionManager")
	public HibernateTransactionManager getTransactionManager()
	{
		HibernateTransactionManager transactionManager = new HibernateTransactionManager();
		transactionManager.setSessionFactory(getSessionFactory().getObject());
		return transactionManager;
	}
	
	// The "dao" bean is a protoype for the DAO-pool. Don't try to use this bean directly, request a DAO from "daoPool" via getTarget();
	@Bean(name = "dao")
	@Scope(scopeName = ConfigurableBeanFactory.SCOPE_PROTOTYPE)
	public GeneCorrelationDAOImpl getDao()
	{
		GeneCorrelationDAOImpl dao = new GeneCorrelationDAOImpl();
		dao.setBatchSize(100000);
		return dao;
	}
	
	@Bean(name = "provenanceDao")
	public ProvenanceDAOImpl getProvenanceDao()
	{
		ProvenanceDAOImpl dao = new ProvenanceDAOImpl();
		return dao;
	}
	
	@Bean(name = "geneDao")
	public GeneDAOImpl getGeneDao()
	{
		GeneDAOImpl dao = new GeneDAOImpl();
		return dao;
	}
	
	// Not sure it's necessary to pool the DAOs anymore. Was needed when running multiple connections for INSERTS, but with bulk-loading from a file
	// with a single connection, this might not be needed anymore.
	@Bean(name = "daoPool")
	public CommonsPool2TargetSource getDaoPool()
	{
		CommonsPool2TargetSource pool = new CommonsPool2TargetSource();
		
		pool.setTargetBeanName("dao");
		pool.setMinIdle(1);
		pool.setMaxSize(8);
		
		return pool;
	}
	
//	@Bean(name= "archs4Loader")
//	public Archs4Loader getArchs4Loader()
//	{
//		return new Archs4Loader(env.getProperty("files.archs4.correlation_file"));
//	}
	
	@Bean(name="correlationMatrixLoader")
	public CorrelationMatrixLoader getCorrelationMatrixLoader()
	{
		return new CorrelationMatrixLoader();
	}
	
	@Bean(name="harmonizomeBatch")
	public HarmonizomeBatch getHarmonizomeBatch()
	{
		return new HarmonizomeBatch();
	}
}