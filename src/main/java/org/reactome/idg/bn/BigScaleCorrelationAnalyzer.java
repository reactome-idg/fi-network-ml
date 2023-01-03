package org.reactome.idg.bn;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.coexpression.CoExpressionLoader;

/**
 * This class is used to perform support evidence analysis using BigScale2 scRNA-seq correlation results. The original
 * correlation file was generated by BigScale2 R script using human Blood scRNA-seq data. For details, see: TODO: add URL
 * to the R script.
 * @author wug
 *
 */
public class BigScaleCorrelationAnalyzer {
	private static final Logger logger = Logger.getLogger(BigScaleCorrelationAnalyzer.class);
	
	private final String DIR_NAME = "results/impact_analysis/scRNASeq/";
	private final String COR_FILE_NAME = DIR_NAME + "tabula_blood_corr_recursive.csv";
	// Use the same ratio as in the original BigScale paper
	private final double PERCENTAGE_CUTOFF = 0.001d; // 0.1%
	
	public BigScaleCorrelationAnalyzer() {
	}
	
	private Set<String> loadCoExpNetwork() throws IOException {
		CoExpressionLoader loader = new CoExpressionLoader();
		loader.setUseAbsoluteValue(false); // Only positive correlations are used. Most likely negative correlations are not
		                                   // reliable considering the strong positive biased distribution. See the Jupyte notebook
		                                   // for the actual distribution.
		return loader.loadCoExpressionViaPercentile(new File(COR_FILE_NAME),
									                PERCENTAGE_CUTOFF);
	}
	
	/**
	 * Perform a pathway enrichment for genes in the network.
	 * @throws Exception
	 */
	@Test
	public void performPathwayEnrichmentAnalysis() throws Exception {
		Set<String> network = loadCoExpNetwork();
		// Need to remove quotation markers in the set
		network = network.stream().map(cor -> cor.replaceAll("\"", "")).collect(Collectors.toSet());
		logger.info("Total size of the network: " + network.size());
		// Need to split it into gene to partners
		Map<String, Set<String>> gene2partners = InteractionUtilities.generateProteinToPartners(network);
		// The following is modified from method attachPathwayEnrichmentResults() in ImpactResultAnalyzer
		// Set up tools for enrichment analysis
		logger.info("Starting pathway enrichment analysis...");
		Map<String, Set<String>> pathwayId2Genes = new ImpactResultAnalyzer().loadPathwayIdToGenes();
		Map<String, Set<String>> gene2PathwayIds = InteractionUtilities.switchKeyValues(pathwayId2Genes);
		PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
		// Perform enrichment analysis
		String srcFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121_with_enrichment_092921.txt";
		String destFile = ImpactResultAnalyzer.RESULT_DIR + "impact_analysis_092121_with_enrichment_092921_with_scRNASeq_082322.txt";
		FileUtilities fu = new FileUtilities();
		fu.setInput(srcFile);
		fu.setOutput(destFile);
		String line = fu.readLine();
		// These results are from pathway enrichment analysis
		line += "\tBigScale_pValue\tBigScale_FDR\tBigScale_HitIds"; // FDR is corrected for one gene here
		fu.printLine(line);
		Map<String, List<GeneSetAnnotation>> gene2annotation = new HashMap<>();
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			String gene = tokens[0];
			if (!gene2partners.keySet().contains(gene)) {
				fu.printLine(line + "\tnan\tnan\tnan"); 
				continue;
			}
			// Test code
//			if (!gene.equals("DICER1"))
//				continue;
			List<GeneSetAnnotation> annotations = gene2annotation.get(gene);
			if (annotations == null) {
				Set<String> interactors = gene2partners.get(gene);
				logger.info("Interacting pathways for " + gene + "...");
				annotations = annotator.annotateGeneSet(interactors, gene2PathwayIds);
				gene2annotation.put(gene, annotations);
			}
			// Find the annotation
			GeneSetAnnotation annotation = null;
			for (GeneSetAnnotation tmp : annotations) {
				if (tmp.getTopic().equals(tokens[1])) {
					annotation = tmp;
					break;
				}
			}
			// In this context, we may see pathways that cannot be mapped into this limited gene set
			if (annotation == null)
				fu.printLine(line + "\tnan\tnan\tnan");
			else {
				line += "\t" + annotation.getPValue() + "\t" + annotation.getFdr() + "\t" + annotation.getHitIds();
				fu.printLine(line);
			}
		}
		fu.close();
		logger.info("Done");
	}
	
}