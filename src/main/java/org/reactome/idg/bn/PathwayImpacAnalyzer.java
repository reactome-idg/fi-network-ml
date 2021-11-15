package org.reactome.idg.bn;

import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.persistence.MySQLAdaptor;
import org.gk.render.Renderable;
import org.gk.render.RenderablePathway;
import org.gk.util.FileUtilities;
import org.gk.util.GKApplicationUtilities;
import org.junit.Test;
import org.reactome.booleannetwork.FuzzyLogicSimulator.ANDGateMode;
import org.reactome.booleannetwork.HillFunction;
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
	private static final Logger logger = Logger.getLogger(PathwayImpacAnalyzer.class);
    // The local WS API started from mvn tomcat7:run
    private final String TCRD_WS_URL = "http://localhost:8060/tcrdws";
    // Cache this to avoid multiple typing
    private ObjectMapper mapper;
    // Cache load predicted FIs for quick performance
    private Map<String, Map<String, Double>> gene2gene2score;
    // Used to map to g values from score
    private Map<Double, Double> score2precision;
    private double[] scores; // Sorted array from the above map for matching prediction scores
    // For output
    private String outFileName = null;

    public PathwayImpacAnalyzer() {
    	mapper = new ObjectMapper();
    }
    
    public static void main(String[] args) {
    	if (args.length < 1) {
    		System.out.println("java -Xmx12G org.reactome.idg.bn.PathwayImpactAnalyzer {out_file_name}");
    		System.exit(0);
    	}
    	PathwayImpacAnalyzer analyzer = new PathwayImpacAnalyzer();
    	analyzer.setOutputFile(args[0]);
    	try {
    		analyzer.performFuzzyLogicAnalysis();
    	}
    	catch(Exception e) {
    		logger.error(e.getMessage(), e);
    	}
    }
    
    public void setOutputFile(String file) {
    	this.outFileName = file;
    }
    
    /**
     * Conduct a systematic impact analysis using the fuzzy logic model.
     * @throws Exception
     * TODO: Add a filter to avoid run simulations for pathways having genes annotated.
     */
    @Test
    public void performFuzzyLogicAnalysis() throws Exception {
//    	outFileName = "test.txt";
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
//    	String gene = "CALHM6"; // Test case
    	
    	PrintStream ps = null;
    	if (outFileName != null)
    		ps = new PrintStream(outFileName);
    	else
    		ps = System.out;
    	// Generate header
    	ps.println("Gene\tDBID\tPathwayName\t"
    			+ "Sum_Activation\tAverage_Activation\t"
    			+ "Sum_Inhibition\tAverage_Inhibition");
    	
    	Set<String> genes = ApplicationConfig.getConfig().getAllGenes();
    	int counter = 0;
    	long time0 = System.currentTimeMillis();
    	for (String gene : genes) {
    		if (!gene.equals("LRFN1"))
    			continue;
    		logger.warn("Working on " + gene + "...");
    		long time1 = System.currentTimeMillis();
    		for (GKInstance pathway : pathways) {
    			if (!pathway.getDBID().equals(6794361L))
    				continue;
    			//    		if (!pathway.getDisplayName().equals("Costimulation by the CD28 family"))
    			//    			continue;
//    			logger.warn("Pathway: " + pathway);
    			PathwayImpactAnalysisResult resultActivation = bnAnalyzer.performDrugImpactAnalysis(pathway, 
    					converter, 
    					gene, 
    					activationMapper);
//    			if (resultActivation == null)
//    				logger.warn("Nothing for activation simulation.");
    			PathwayImpactAnalysisResult resultInhibition = bnAnalyzer.performDrugImpactAnalysis(pathway,
    					converter,
    					gene, 
    					inhibitionMapper);
//    			if (resultInhibition == null)
//    				logger.warn("Nothing for inhibition simulation");
    			if (resultActivation == null && resultInhibition == null)
    				continue;
    			ps.println(gene + "\t" + pathway.getDBID() + "\t" + pathway.getDisplayName() + "\t" + 
    					(resultActivation == null ? "" : resultActivation.getSum()) + "\t" + 
    					(resultActivation == null ? "" : resultActivation.getAverage()) + "\t" + 
    					(resultInhibition == null ? "" : resultInhibition.getSum()) + "\t" + 
    					(resultInhibition == null ? "" : resultInhibition.getAverage()));
    		}
    		long time2 = System.currentTimeMillis();
    		logger.warn("Done gene: " + gene + " (" + (time2 - time1) / 1000.0d + " seconds)");
    		counter ++;
//    		if (counter == 100)
//    			break;
    	}
    	logger.warn("Total time for simulation: " + (System.currentTimeMillis() - time0) / 1000.0d + " seconds.");
    	ps.close();
    }
    
    /**
     * This method is used to map pathways subject to simulation to the top level pathways for visualization.
     * @throws Exception
     */
    @Test
    public void mapPathwaysToTopLevels() throws Exception {
    	MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
    	List<GKInstance> pathways = loadPathwaysForAnalysis(dba);
    	// Pull out the top level pathways
    	Long frontPageId = 9731223L;
    	GKInstance frontPage = dba.fetchInstance(frontPageId);
    	List<GKInstance> topPathways = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
    	Map<GKInstance, Set<GKInstance>> pathway2top = new HashMap<>();
    	for (GKInstance topPathway : topPathways) {
    		Set<GKInstance> comps = InstanceUtilities.getContainedEvents(topPathway);
    		comps.add(topPathway); // Add itself in case the top is small enough
    		for (GKInstance comp : comps) {
    			pathway2top.compute(comp, (key, set) -> {
    				if (set == null)
    					set = new HashSet<>();
    				set.add(topPathway);
    				return set;
    			});
    		}
    	}
    	
    	String output = ImpactResultAnalyzer.RESULT_DIR + "pathway2topic_100721.txt";
    	FileUtilities fu = new FileUtilities();
    	fu.setOutput(output);
    	fu.printLine("Pathway\tTopic");
    	for (GKInstance pathway : pathways) {
    		Set<GKInstance> tops = pathway2top.get(pathway);
    		if (tops == null)
    			throw new IllegalStateException("Cannot find a top pathway for: " + pathway);
    		String topNames = null;
    		if (tops.size() > 1) {
    			topNames = tops.stream().map(t -> t.getDisplayName().trim()).sorted().collect(Collectors.joining(","));
    		}
    		else
    			topNames = tops.stream().findAny().get().getDisplayName().trim(); // Have to do trim
    		fu.printLine(pathway.getDisplayName().trim()  + "\t" + topNames);
    	}
    	fu.close();
    }
    
    /**
     * The implementation of this method is ported from org.reactome.r3.fi.ReactomeObjectHandler.loadPathwayDBIDs()
     * in the FI core WS project. The goal is to get a list of pathways that have ELVs for impact analysis.
     * @param mysqlDBA
     * @return
     * @throws Exception
     */
    protected List<GKInstance> loadPathwaysForAnalysis(MySQLAdaptor mysqlDBA) throws Exception {
        Collection<GKInstance> diagrams = mysqlDBA.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram);
        DiagramGKBReader reader = new DiagramGKBReader();
        List<GKInstance> rtn = new ArrayList<GKInstance>();
        for (GKInstance diagram : diagrams) {
        	// Make sure it is for human pathway
        	List<GKInstance> pathways = diagram.getAttributeValuesList(ReactomeJavaConstants.representedPathway);
        	// For normal pathway only
        	boolean isHuman = false;
        	GKInstance pathwayInst = null;
        	for (int i = 0; i < pathways.size(); i++) {
        		pathwayInst = pathways.get(i);
        		GKInstance species = (GKInstance) pathwayInst.getAttributeValue(ReactomeJavaConstants.species);
        		if (species.getDBID().equals(GKApplicationUtilities.HOMO_SAPIENS_DB_ID)) {
        			isHuman = true;
        			break;
        		}
        	}
        	if (!isHuman)
        		continue;
            boolean hasELV = false;
            RenderablePathway pathway = reader.openDiagram(diagram);
            List<Renderable> components = pathway.getComponents();
            if (components != null) {
                for (Renderable r : components) {
                    if (r.getReactomeId() == null)
                        continue;
                    GKInstance inst = mysqlDBA.fetchInstance(r.getReactomeId());
                    if (inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                        hasELV = true;
                        break;
                    }
                }
            }
            if (hasELV) 
                rtn.add(pathwayInst);
        }
        return rtn;
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
    
    public DrugToTargetsMapper getInteractionMapper() {
    	return new ActivationGeneFITargetMapper();
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
