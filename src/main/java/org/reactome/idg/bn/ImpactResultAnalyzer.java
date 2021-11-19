package org.reactome.idg.bn;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.data.ReactomeAnalyzerTopicHelper;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.DatabaseConfig;
import org.reactome.pathway.booleannetwork.DrugToTargetsMapper;


@SuppressWarnings("unchecked")
public class ImpactResultAnalyzer {
	public static final String RESULT_DIR = "results/impact_analysis/";
	private static final Logger logger = Logger.getLogger(ImpactResultAnalyzer.class);
	
	public ImpactResultAnalyzer() {
	}
	
	/**
	 * Files generated in this method are used for performing some MeSH term based NLP analysis via collaborating
	 * with Augustin Luna at Dana Faber.
	 * @throws IOException
	 */
	@Test
	public void dumpFilesForNLP() throws Exception {
		String outDir = RESULT_DIR + "nlp_files/";
		File dir = new File(outDir);
		if (!dir.exists())
			dir.mkdir();
		FileUtilities fu = new FileUtilities();
		// Dump the pathways used for simulation and interaction analysis
		MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
		PathwayImpactAnalyzer impactAnalyzer = new PathwayImpactAnalyzer();
		List<GKInstance> pathways = impactAnalyzer.loadPathwaysForAnalysis(dba);
		InstanceUtilities.sortInstances(pathways);
	    dba.loadInstanceAttributeValues(pathways, new String[]{ReactomeJavaConstants.stableIdentifier});
//		String fileName = outDir + "AnalyzedPathways_111512.txt";
//		fu.setOutput(fileName);
//		fu.printLine("Pathway_Name\tStableId\tDB_ID");
//		for (GKInstance pathway : pathways) {
//			String id = getStableId(pathway);
//			fu.printLine(pathway.getDisplayName() + "\t" + 
//					     id + "\t" + 
//					     pathway.getDBID());
//		}
//		fu.close();
		// Dump the associated PMIDs for the above pathways
		String fileName = outDir + "PathwayId2PMID_111512.txt";
		fu.setOutput(fileName);
		dumpPMIDs(fu, pathways, dba);
		// Dump the dark proteins in the IDG three protein families
//		Map<String, String> target2levels = loadProtein2DevLevel();
//		Map<String, String> target2family = loadProtein2Family();
//		Set<String> selectedFamilies = Stream.of("IC", "Kinase", "GPCR").collect(Collectors.toSet());
//		List<String> targets = target2levels.keySet()
//				.stream()
//				.filter(t -> target2levels.get(t).equals("Tdark"))
//				.filter(t -> selectedFamilies.contains(target2family.get(t)))
//				.sorted()
//				.collect(Collectors.toList());
//		String fileName = outDir + "DarkProteins_3_families_111512.txt";
//		dumpTargets(targets, target2levels, target2family, fileName, fu);
//		// Randomly sample the same number of non-dark proteins 
//		List<String> sampledTargets = target2levels.keySet()
//				.stream()
//				.filter(t -> !target2levels.get(t).equals("Tdark"))
//				.collect(Collectors.toList());
//		Set<String> set = MathUtilities.randomSampling(sampledTargets, targets.size());
//		fileName = outDir + "NonDarkProteins_Sampled_111521.txt";
//		dumpTargets(set.stream().sorted().collect(Collectors.toList()),
//				target2levels,
//				target2family,
//				fileName, 
//				fu);
//		// Randomly sample the same number of dark proteins
//		sampledTargets = target2levels.keySet()
//				.stream()
//				.filter(t -> target2levels.get(t).equals("Tdark"))
//				.collect(Collectors.toList());
//		set = MathUtilities.randomSampling(sampledTargets, targets.size());
//		fileName = outDir + "DarkProteins_Sampled_111521.txt";
//		dumpTargets(set.stream().sorted().collect(Collectors.toList()),
//				target2levels,
//				target2family,
//				fileName, 
//				fu);
	}

	private void dumpPMIDs(FileUtilities fu,
	                       List<GKInstance> pathways,
	                       MySQLAdaptor dba) throws Exception {
		fu.printLine("Pathway_Stable_ID\tPMID");
		Set<Integer> totalPMIDs = new HashSet<>();
		for (GKInstance pathway : pathways) {
			List<Integer> pmids = loadPathwayPMIDs(pathway);
			totalPMIDs.addAll(pmids);
			String stableId = getStableId(pathway);
			for (Integer pmid : pmids)
				fu.printLine(stableId + "\t" + pmid);
		}
		fu.close();
		System.out.println("Total PMIDs: " + totalPMIDs.size());
		// Check PMIDs not in the pathway list
		Collection<GKInstance> references = dba.fetchInstancesByClass(ReactomeJavaConstants.LiteratureReference);
		dba.loadInstanceAttributeValues(references, new String[] {ReactomeJavaConstants.pubMedIdentifier});
		int c = 0;
		for (GKInstance reference : references) {
			Integer pmid = (Integer) reference.getAttributeValue(ReactomeJavaConstants.pubMedIdentifier);
			if (totalPMIDs.contains(pmid))
				continue;
			System.out.println(reference + ":");
			Collection<GKInstance> referrers = reference.getReferers(ReactomeJavaConstants.literatureReference);
			referrers.forEach(r -> System.out.println("\t" + r));
			c++;
			if (c == 20)
				break;
		}
	}

	private void dumpTargets(List<String> targets,
	                         Map<String, String> target2levels,
	                         Map<String, String> target2family,
	                         String fileName,
	                         FileUtilities fu)
	        throws IOException {
		fu.setOutput(fileName);
		fu.printLine("Protein\tFamily\tDevLevel");
		for (String target : targets) {
			fu.printLine(target + "\t" + 
		                 target2family.get(target) + "\t" + 
					     target2levels.get(target));
		}
		fu.close();
		System.out.println("Total targets in " + fileName + ": " + targets.size());
	}

	private String getStableId(GKInstance pathway) throws InvalidAttributeException, Exception {
		GKInstance stableId = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
		String id = (String) stableId.getAttributeValue(ReactomeJavaConstants.identifier);
		return id;
	}
	
	private List<Integer> loadPathwayPMIDs(GKInstance pathway) throws Exception {
		// Grep all event instances
		Set<GKInstance> pathwayEvents = InstanceUtilities.getContainedEvents(pathway);
		pathwayEvents.add(pathway);
		Set<GKInstance> references = new HashSet<>();
		for (GKInstance event : pathwayEvents) {
			loadPathwayPMIDs(event, references);
		}
		// References in the inferred from should be fetched
		for (GKInstance event : pathwayEvents) {
			List<GKInstance> inferredFromEvents = event.getAttributeValuesList(ReactomeJavaConstants.inferredFrom);
			if (inferredFromEvents == null || inferredFromEvents.size() == 0)
				continue;
			for (GKInstance inferrdFromEvent : inferredFromEvents) {
				loadPathwayPMIDs(inferrdFromEvent, references);
			}
		}
		List<Integer> rtn = new ArrayList<>();
		for (GKInstance reference : references) {
			// Some references may not have pubMedId
			if (!reference.getSchemClass().isValidAttribute(ReactomeJavaConstants.pubMedIdentifier))
				continue;
			Integer pmid = (Integer) reference.getAttributeValue(ReactomeJavaConstants.pubMedIdentifier);
			if (pmid == null) {
				System.err.println(reference + " doesn't have a PMID!");
				continue;
			}
			rtn.add(pmid);
		}
		return rtn.stream().sorted().collect(Collectors.toList());
	}

	private void loadPathwayPMIDs(GKInstance event, Set<GKInstance> references)
	        throws InvalidAttributeException, Exception {
		if (event.getSchemClass().isValidAttribute(ReactomeJavaConstants.literatureReference)) {
			List<GKInstance> litRefs = event.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
			if (litRefs != null && litRefs.size() > 0)
				references.addAll(litRefs);
		}
		grepReferences(event, "catalystActivityReference", references);
		grepReferences(event, "regulationReference", references);
		grepReferences(event, "summation", references);
	}
	
	private void grepReferences(GKInstance event, 
	                            String attributeName,
	                            Set<GKInstance> references) throws Exception {
		if (event.getSchemClass().isValidAttribute(attributeName)) {
			List<GKInstance> values = event.getAttributeValuesList(attributeName);
			if (values != null) {
				for (GKInstance value : values) {
					if (!value.getSchemClass().isValidAttribute(ReactomeJavaConstants.literatureReference))
						continue;
					List<GKInstance> litRefs = value.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
					if (litRefs != null && litRefs.size() > 0)
						references.addAll(litRefs);
				}
			}
		}
	}
	
	/**
	 * Generate the data frames for some R-based visualization.
	 * @throws IOException
	 */
	@Test
	public void createDataFrame() throws IOException {
		// Used for filtering
		Map<String, String> gene2level = loadProtein2DevLevel();
		Map<String, String> gene2family = loadProtein2Family();
		// Selected proteins
		Set<String> selectedFamilies = Stream.of("IC", "Kinase", "GPCR").collect(Collectors.toSet());
		
		String src = RESULT_DIR + "impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt";
		
		// Output the data frame for FDR
//		String out = RESULT_DIR + "impact_analysis_092121_df_fdr_no_filter_three_families_all_100821.txt";
//		int index = 10;
//		boolean needNegLog = true;
//		double threshold = 1.0d;
//		boolean isFDR = true;
		// Output for Average_Activation
//		String out = RESULT_DIR + "impact_analysis_092121_df_activation_three_families_all_100821.txt";
//		int index = 4;
//		boolean needNegLog = false;
//		double threshold = 0.0d;
//		boolean isFDR = false;
		// Output for average_inhibition
		String out = RESULT_DIR + "impact_analysis_092121_df_inhibition_three_families_all_100821.txt";
		int index = 6;
		boolean needNegLog = false;
		double threshold = 0.0d;
		boolean isFDR = false;
		
				
		Set<String> allPathways = new HashSet<>();
		Map<String, Map<String, Double>> gene2pathway2score = new HashMap<>();
		FileUtilities fu = new FileUtilities();
		fu.setInput(src);
		String line = fu.readLine();
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			String gene = tokens[0];
//			if (!gene2level.get(gene).equals("Tdark"))
//				continue;
			if (!selectedFamilies.contains(gene2family.get(gene)))
				continue;
			Double score = new Double(tokens[index]);
			if (isFDR && score > threshold)
				continue;
			if (!isFDR && score < threshold)
				continue; // Other cases
			Map<String, Double> pathway2score = gene2pathway2score.get(gene);
			if (pathway2score == null) {
				pathway2score = new HashMap<>();
				gene2pathway2score.put(gene, pathway2score);
			}
			if (needNegLog)
				score = -Math.log10(score);
			pathway2score.put(tokens[2], score);
			allPathways.add(tokens[2]);
		}
		fu.close();
		
		exportDataFrame(gene2pathway2score, allPathways, out, fu);
	}

	private void exportDataFrame(Map<String, Map<String, Double>> gene2pathway2score,
	                             Set<String> allPathways,
	                             String out,
	                             FileUtilities fu) throws IOException {
		System.out.println("Total selected genes: " + gene2pathway2score.size());
		System.out.println("Total selected patwhays: " + allPathways.size());
		List<String> pathwayList = allPathways.stream().sorted().collect(Collectors.toList());
		fu.setOutput(out);
		StringBuilder builder = new StringBuilder();
		builder.append("Gene");
		pathwayList.forEach(p -> builder.append("\t").append(p));
		fu.printLine(builder.toString());
		List<String> geneList = gene2pathway2score.keySet().stream().sorted().collect(Collectors.toList());
		for (String gene : geneList) {
			Map<String, Double> pathway2score = gene2pathway2score.get(gene);
			builder.setLength(0);
			builder.append(gene);
			for (String pathway : pathwayList) {
				Double score = pathway2score.get(pathway);
				if (score == null)
					score = 0.0d;
				builder.append("\t").append(score);
			}
			fu.printLine(builder.toString());
		}
		fu.close();
	}
	
	private Map<String, String> loadProtein2Family() throws IOException {
		InputStream is = DatabaseConfig.class.getClassLoader().getResourceAsStream("UniProtGeneDevLevelsTypes_100721.txt");
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line = br.readLine();
		Map<String, String> gene2family = new HashMap<>();
		while ((line = br.readLine()) != null) {
			String[] tokens = line.split("\t");
			gene2family.put(tokens[1], tokens[3]);
		}
		br.close();
		return gene2family;
	}
	
	private Map<String, String> loadProtein2DevLevel() throws IOException {
		InputStream is = DatabaseConfig.class.getClassLoader().getResourceAsStream("UniProtGeneDevLevelsTypes_100721.txt");
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line = br.readLine();
		Map<String, String> gene2level = new HashMap<>();
		while ((line = br.readLine()) != null) {
			String[] tokens = line.split("\t");
			gene2level.put(tokens[1], tokens[2]);
		}
		br.close();
		return gene2level;
	}
	
	@Test
	public void attachDevLevels() throws Exception {
		Map<String, String> gene2level = loadProtein2DevLevel();
		System.out.println("Total gen2level: " + gene2level.size());
//		String srcFile = RESULT_DIR + "neuronal_system_impact_analysis_092121_with_enrichment_092921_100521.txt";
//		String destFile = RESULT_DIR + "neuronal_system_impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt";
		
		String srcFile = RESULT_DIR + "impact_analysis_092121_with_enrichment_092921.txt";
		String destFile = RESULT_DIR + "impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt";
		
		FileUtilities fu = new FileUtilities();
		fu.setInput(srcFile);
		fu.setOutput(destFile);
		String line = fu.readLine();
		fu.printLine(line + "\tDevLevel");
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			String level = gene2level.get(tokens[0]);
			fu.printLine(line + "\t" + level);
		}
		fu.close();
	}
	
	@Test
	public void filterToTopPathway() throws Exception {
		MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
		// Top pathway
		Long dbId = 112316L; // Neuronal System
		dbId = 6794362L; // Protein-protein interactions at synapses
		GKInstance pathway = dba.fetchInstance(dbId);
		Set<GKInstance> containedPathways = InstanceUtilities.getContainedEvents(pathway);
		containedPathways.add(pathway);
		logger.info("Total pathways: " + containedPathways.size());
		Set<Long> pathwayIds = containedPathways.stream().map(p -> p.getDBID()).collect(Collectors.toSet());
		
//		String srcFile = RESULT_DIR + "impact_analysis_092121_with_enrichment_092921.txt";
//		String destFile = RESULT_DIR + "neuronal_system_impact_analysis_092121_with_enrichment_092921_100521.txt";
		
		String srcFile = RESULT_DIR + "impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt";
		String destFile = RESULT_DIR + "protein_interactions_at_synapses_impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt";
		
		FileUtilities fu = new FileUtilities();
		fu.setInput(srcFile);
		fu.setOutput(destFile);
		String line = fu.readLine();
		fu.printLine(line);
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			if (pathwayIds.contains(new Long(tokens[1])))
				fu.printLine(line);
		}
		fu.close();
	}
	
	@Test
	public void attachPathwayEnrichmentResults() throws Exception {
		// Set up tools for enrichment analysis
		Map<String, Set<String>> pathwayId2Genes = loadPathwayIdToGenes();
		Map<String, Set<String>> gene2PathwayIds = InteractionUtilities.switchKeyValues(pathwayId2Genes);
		PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
		DrugToTargetsMapper interactionMapper = new PathwayImpactAnalyzer().getInteractionMapper();
		// Perform enrichment analysis
		String srcFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121.txt";
		String destFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121_with_enrichment_092921.txt";
		destFile = ImpactResultAnalyzer.RESULT_DIR + "test_delete_me.txt";
		FileUtilities fu = new FileUtilities();
		fu.setInput(srcFile);
		fu.setOutput(destFile);
		String line = fu.readLine();
		// These results are from pathway enrichment analysis
		line += "\tIsAnnotated\tMapped_Genes\tpValue\tFDR"; // FDR is corrected for one gene here
		fu.printLine(line);
		Map<String, List<GeneSetAnnotation>> gene2annotation = new HashMap<>();
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			
			// Test code
			if (!tokens[0].equals("LRFN1"))
				continue;
			
			List<GeneSetAnnotation> annotations = gene2annotation.get(tokens[0]);
			if (annotations == null) {
				Set<String> interactors = interactionMapper.getDrugTargets(tokens[0]);
				logger.info("Interacting pathways for " + tokens[0] + "...");
				annotations = annotator.annotateGeneSet(interactors, gene2PathwayIds);
				gene2annotation.put(tokens[0], annotations);
			}
			// Find the annotation
			GeneSetAnnotation annotation = null;
			for (GeneSetAnnotation tmp : annotations) {
				if (tmp.getTopic().equals(tokens[1])) {
					annotation = tmp;
					break;
				}
			}
			if (annotation == null)
				throw new IllegalStateException("No pathway annotation for: " + line);
			Set<String> pathwayGenes = pathwayId2Genes.get(tokens[1]);
			boolean isAnnotated = pathwayGenes.contains(tokens[0]);
			line += "\t" + isAnnotated + "\t" + annotation.getHitNumber() + "\t" + 
			        annotation.getPValue() + "\t" + annotation.getFdr() + "\t" + annotation.getHitIds();
			fu.printLine(line);
		}
		fu.close();
		logger.info("Done");
	}
	
	@Test
	public void checkResults() throws Exception {
		String resultFile = RESULT_DIR + "impact_analysis_092121_filtered_092321.txt";
		Set<String> totalGenes = new HashSet<String>();
		Set<String> totalPathways = new HashSet<String>();
		FileUtilities fu = new FileUtilities();
		fu.setInput(resultFile);
		String line = fu.readLine();
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			totalGenes.add(tokens[0]);
			totalPathways.add(tokens[2]);
		}
		fu.close();
		System.out.println("Total genes: " + totalGenes.size());
		System.out.println("Total pathways: " + totalPathways.size());
	}
	
	/**
	 * During simulation, pathways that have genes annotated are subject to simulation too. In our final analysis,
	 * we will assume genes play significant roles in these annotated pathways without simulation. This method is 
	 * used to filter pathways in the report to make the report simpler.
	 * @throws Exception
	 */
	@Test
	public void filterOutAnnotatedPathways() throws Exception {
		Map<String, Set<String>> pathwayId2Genes = loadPathwayIdToGenes();
		System.out.println("Size of pathways2genes: " + pathwayId2Genes.size());
		
		String src = RESULT_DIR + "impact_analysis_092121.txt";
		String dest = RESULT_DIR + "impact_analysis_092121_filtered_092321.txt";
		dest = RESULT_DIR + "impact_analysis_092121_filtered_092921.txt";
		FileUtilities fu = new FileUtilities();
		fu.setInput(src);
		fu.setOutput(dest);
		String line = fu.readLine();
		fu.printLine(line);
		int totalFilteredLines = 0;
		int totalLines = 0;
		while ((line = fu.readLine()) != null) {
			totalLines ++;
			String[] tokens = line.split("\t");
			Set<String> pathwayGenes = pathwayId2Genes.get(tokens[1]);
			if (pathwayGenes.contains(tokens[0])) {
				totalFilteredLines ++;
				continue;
			}
			fu.printLine(line);
		}
		fu.close();
		System.out.println("Total filtered lines: " + totalFilteredLines);
		System.out.println("Total lines: " + totalLines);
		double ratio = (double) totalFilteredLines / totalLines;
		System.out.println("Ratio of filtered to total: " + ratio);
	}
	
	/**
	 * A simple way to grep genes for individual pathways from the MySQL database used for simulation.
	 * @return
	 * @throws Exception
	 */
	private Map<String, Set<String>> loadPathwayIdToGenes() throws Exception {
		MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
		PathwayImpactAnalyzer impactAnalyzer = new PathwayImpactAnalyzer();
		List<GKInstance> pathways = impactAnalyzer.loadPathwaysForAnalysis(dba);
		ReactomeAnalyzerTopicHelper topicHelper = new ReactomeAnalyzerTopicHelper();
		Map<String, Set<String>> pathwayId2Genes = new HashMap<>();
		for (GKInstance pathway : pathways) {
			Set<GKInstance> participants = topicHelper.grepPathwayParticipants(pathway);
			Set<String> genes = new HashSet<String>();
			for (GKInstance participant : participants) {
				Set<GKInstance> refGenes = grepRefPepSeqs(participant);
				for (GKInstance refGene : refGenes) {
					if (refGene.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
						String gene = (String) refGene.getAttributeValue(ReactomeJavaConstants.geneName);
						if (gene != null)
							genes.add(gene);
					}
				}
			}
			if (genes.size() == 0)
				continue; // Just in case there is no genes in a pathway
			pathwayId2Genes.put(pathway.getDBID().toString(), genes);
		}
		return pathwayId2Genes;
	}
	
	/**
	 * The following method basic is copied from ReactomeAnalyzerTopicHelper. However, in this implementation,
	 * hasCandidate is added back. The reason is to make it consistent with impact analysis since hasCandidate
	 * is converted there into Boolean networks. Also the web site uses the mapping file which includes hasCandidate.
	 * However, the FIs extracted from Reactome reactions and complexes, which are used as the positive training dataset,
	 * have not considered these candidates still. This may cause the inconsistence and confusion. We may need to change
	 * later on in the future.
	 */
    private Set<GKInstance> grepRefPepSeqs(GKInstance pe) throws Exception {
    	Set<GKInstance> refSeqs = new HashSet<>();
        Set<GKInstance> ewases = InstanceUtilities.getContainedInstances(pe,
                                                                         ReactomeJavaConstants.hasComponent,
                                                                         ReactomeJavaConstants.hasMember,
                                                                         ReactomeJavaConstants.hasCandidate);
        ewases.add(pe);
        for (GKInstance ewas : ewases) {
            if (!ewas.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                continue;
            GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refEntity == null)
                continue;
            if (refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct) ||
                refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence)) {
                refSeqs.add(refEntity);
            }
        }
        return refSeqs;
    }

}
