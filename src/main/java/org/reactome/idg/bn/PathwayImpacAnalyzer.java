package org.reactome.idg.bn;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.hc.client5.http.fluent.Request;
import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.booleannetwork.HillFunction;
import org.reactome.booleannetwork.FuzzyLogicSimulator.ANDGateMode;
import org.reactome.idg.util.ApplicationConfig;
import org.reactome.idg.util.DatabaseConfig;
import org.reactome.pathway.booleannetwork.BNPerturbationAnalyzer;
import org.reactome.pathway.booleannetwork.DrugToTargetsMapper;
import org.reactome.pathway.booleannetwork.PathwayImpactAnalysisResult;
import org.reactome.pathway.booleannetwork.PathwayToBooleanNetworkConverter;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * This class is used to do a comprehensive analysis of dark proteins in the context of Reactome pathways. To run methods
 * in this class, this web service app should be up and running. The RESTful API is called to get all information needed.
 */
@SuppressWarnings("unchecked")
public class PathwayImpacAnalyzer {
	private final Logger logger = Logger.getLogger(PathwayImpacAnalyzer.class);
    // The local WS API started from mvn tomcat7:run
    private final String TCRD_WS_URL = "http://localhost:8060/tcrdws";
    private final String FI_CORE_WS_URL = "http://localhost:8080/corews";
    // Cache this to avoid multiple typing
    private ObjectMapper mapper;
    // Cache load predicted FIs for quick performance
    private Map<String, Map<String, Double>> gene2gene2score;
    // Used to map to g values from score
    private Map<Double, Double> score2precision;
    private double[] scores; // Sorted array from the above map for matching prediction scores

    public PathwayImpacAnalyzer() {
    	mapper = new ObjectMapper();
    }
    
    /**
     * Conduct a systematic impact analysis using the fuzzy logic model.
     * @throws Exception
     */
    @Test
    public void performFuzzyLogicAnalysis() throws Exception {
    	MySQLAdaptor mysqlDBA = DatabaseConfig.getMySQLDBA();
    	// Get the list of pathways having ELV
    	List<GKInstance> pathways = loadPathwaysForAnalysis(mysqlDBA);
    	logger.info("Total pathways from MySQL database loaded: " + pathways.size());
    	// We will port class org.reactome.r3.fi.PathwayDrugImpactAnalyzer for doing this. 
    	BNPerturbationAnalyzer bnAnalyzer = new BNPerturbationAnalyzer();
    	// Configure
    	String bnDefaultValue = ApplicationConfig.getConfig().getAppConfig("bn.default.value");
    	bnAnalyzer.setDefaultValue(Double.parseDouble(bnDefaultValue));
    	bnAnalyzer.setTransferFunction(new HillFunction());
    	bnAnalyzer.setAndGateMode(ANDGateMode.PROD); // Use product so that we can consider better of activation.
    	
    	GeneFITargetMapper activationMapper = new ActivationGeneFITargetMapper();
    	GeneFITargetMapper inhibitionMapper = new InhibitionGeneFITargetMapper();
    	PathwayToBooleanNetworkConverter converter = new PathwayToBooleanNetworkConverter();
    	String gene = "CALHM6"; // Test case
    	for (GKInstance pathway : pathways) {
    		if (!pathway.getDisplayName().equals("Costimulation by the CD28 family"))
    			continue;
    		System.out.println("Pathway: " + pathway);
    		PathwayImpactAnalysisResult result = bnAnalyzer.performDrugImpactAnalysis(pathway, 
    																			      converter, 
    																			      gene, 
    																			      activationMapper);
    		if (result == null)
    			continue;
    		System.out.println("Activation Average: " + result.getAverage());
    		result = bnAnalyzer.performDrugImpactAnalysis(pathway,
    				converter,
    				gene, 
    				inhibitionMapper);
    		if (result == null)
    			continue;
    		System.out.println("Inhibition Average: " + result.getAverage());
    		break;
    	}
    }
    
    private List<GKInstance> loadPathwaysForAnalysis(MySQLAdaptor mysqlDBA) throws Exception {
    	String url = FI_CORE_WS_URL + "/FIService/network/pathways/logicmodels/hsa";
    	String text = Request.get(url).execute().returnContent().asString();
    	List<String> pathwayStIds = Stream.of(text.split("\n")).collect(Collectors.toList());
    	logger.info("Total pathways to be checked: " + pathwayStIds.size());
    	List<GKInstance> pathways = fetchPathways(pathwayStIds, mysqlDBA);
    	return pathways;
    }
    
    /**
     * Fetch pathways in GKInstance for a list of pathway stable ids.
     * @param stIds
     * @param dba
     * @return
     * @throws Exception
     */
    private List<GKInstance> fetchPathways(List<String> stIds, MySQLAdaptor dba) throws Exception {
    	List<GKInstance> pathways = new ArrayList<GKInstance>();
    	for (String stId : stIds) {
    		// Get the stable id first
    		Collection<GKInstance> c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.StableIdentifier,
    				ReactomeJavaConstants.identifier,
    				"=", 
    				stId);
    		if (c == null || c.size() == 0)
    			throw new IllegalStateException(stId + " cannot be found in the MySQL database.");
    		// There should be only one 
    		GKInstance stIdInst = c.stream().findAny().get();
    		// Get the pathway for this StableID
    		c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
    				ReactomeJavaConstants.stableIdentifier,
    				"=",
    				stIdInst);
    		if (c == null || c.size() == 0)
    			throw new IllegalStateException("Cannot find a pathway for " + stIdInst);
    		pathways.add(c.stream().findAny().get());
    	}
    	return pathways;
    }

    /**
     * Get a list of dark proteins. To run this method, make sure the tcrdws is running.
     * @return
     * @throws Exception
     */
    public List<String> getDarkProteins() throws Exception {
        String url = TCRD_WS_URL + "/tdark/uniprots";
        URL urlAddress = new URL(url);
        // For get-based URL, we may use this simple API
        List<String> list = mapper.readValue(urlAddress, new TypeReference<List<String>>(){});
        return list;
    }
    
    private abstract class GeneFITargetMapper implements DrugToTargetsMapper {
    	
    	public GeneFITargetMapper() {
    	}
    	
    	protected Map<String, Map<String, Double>> getGene2Gene2score() throws IOException {
    		if (gene2gene2score == null) {
				String threshold = ApplicationConfig.getConfig().getAppConfig("score.threshold.for.bn");
				gene2gene2score = new PredictionScoreAnalyzer().loadPredictionFIFile(Double.parseDouble(threshold));
			}
    		return gene2gene2score;
		}

		@Override
		public Set<String> getDrugTargets(String gene) throws Exception {
			Map<String, Map<String, Double>> gene2gene2score = getGene2Gene2score();
			Map<String, Double> gene2score = gene2gene2score.get(gene);
			if (gene2score == null)
				return Collections.EMPTY_SET;
			return new HashSet<String>(gene2score.keySet());
		}
		
		protected Map<String, Double> getGeneToStrength(Collection<String> genes, String gene) throws Exception {
			Map<String, Map<String, Double>> gene2gene2score = getGene2Gene2score();
			Map<String, Double> gene2score = gene2gene2score.get(gene);
			if (gene2score == null)
				return Collections.EMPTY_MAP;
			return convertScoreToG(gene2score);
		}
		
		/**
		 * This method is used to convert prediction scores to g values to be used in fuzzy logic simulation. In the
		 * current implementation, this is done by multiplying score by precision.
		 * @param gene2score
		 * @return
		 */
		protected Map<String, Double> convertScoreToG(Map<String, Double> gene2score) throws IOException {
			if (score2precision == null) {
				loadScore2Precision();
			}
			return gene2score.keySet().stream().collect(Collectors.toMap(Function.identity(),
					g -> {
						Double score = gene2score.get(g);
						// Get the matched scores
						int found = Arrays.binarySearch(scores, score);
						if (found >= 0) return scores[found];
						int insertionPoint = -found - 1; 
						if (insertionPoint == 0) return scores[0];
						if (insertionPoint == scores.length) return scores[scores.length - 1];
						// Do a comparison
						double v1 = scores[insertionPoint - 1];
						double v2 = scores[insertionPoint];
						double selected = score - v1 < v2 - score ? v1 : v2;
						return score2precision.get(selected) * selected;
					}));
		}

		private void loadScore2Precision() throws IOException {
			// Load scor2precision
			String fileName = ApplicationConfig.getConfig().getAppConfig("precision.recall.file");
			try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
				score2precision = lines.skip(1)
								       .map(line -> line.split(","))
								       .collect(Collectors.toMap(tokens -> new Double(tokens[3]), 
								    		                     tokens -> new Double(tokens[1])));
				scores = new double[score2precision.size()];
				int index = 0;
				for (Double d : score2precision.keySet()) {
					scores[index ++] = d;
				}
				Arrays.sort(scores);
			}
		}
    }
    
    private class ActivationGeneFITargetMapper extends GeneFITargetMapper {

		@Override
		public Map<String, Double> getGeneToInhibition(Collection<String> genes, String gene) throws Exception {
			return Collections.EMPTY_MAP;
		}

		@Override
		public Map<String, Double> getGeneToActivation(Collection<String> genes, String gene) throws Exception {
			Map<String, Double> rtn = getGeneToStrength(genes, gene);
//			rtn.keySet().retainAll(Stream.of("LCK").collect(Collectors.toSet()));
			return rtn;
		}
    	
    }
    
    private class InhibitionGeneFITargetMapper extends GeneFITargetMapper {

		@Override
		public Map<String, Double> getGeneToInhibition(Collection<String> genes, String gene) throws Exception {
			return getGeneToStrength(genes, gene);
		}

		@Override
		public Map<String, Double> getGeneToActivation(Collection<String> genes, String gene) throws Exception {
			return Collections.EMPTY_MAP;
		}
    	
    }

}
