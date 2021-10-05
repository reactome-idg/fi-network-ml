package org.reactome.idg.bn;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.data.ReactomeAnalyzerTopicHelper;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.DatabaseConfig;
import org.reactome.pathway.booleannetwork.DrugToTargetsMapper;

public class ImpactResultAnalyzer {
	private static final String RESULT_DIR = "results/impact_analysis/";
	private static final Logger logger = Logger.getLogger(ImpactResultAnalyzer.class);
	
	public ImpactResultAnalyzer() {
	}
	
	@Test
	public void attachPathwayEnrichmentResults() throws Exception {
		// Set up tools for enrichment analysis
		Map<String, Set<String>> pathwayId2Genes = loadPathwayIdToGenes();
		Map<String, Set<String>> gene2PathwayIds = InteractionUtilities.switchKeyValues(pathwayId2Genes);
		PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
		DrugToTargetsMapper interactionMapper = new PathwayImpacAnalyzer().getInteractionMapper();
		// Perform enrichment analysis
		String srcFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121.txt";
		String destFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121_with_enrichment_092921.txt";
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
			        annotation.getPValue() + "\t" + annotation.getFdr();
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
		PathwayImpacAnalyzer impactAnalyzer = new PathwayImpacAnalyzer();
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
