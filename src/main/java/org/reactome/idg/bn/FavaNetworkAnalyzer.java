package org.reactome.idg.bn;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.ApplicationConfig;
import org.reactome.r3.util.FileUtility;

/**
 * We hope to use the functional interaction network generated by FAVA, a variational autoencoder-based approach,
 * to find support evidence for our interacting pathway scores. The original document for FAVA is here:
 * https://github.com/mikelkou/fava. 
 * @author wug
 *
 */
public class FavaNetworkAnalyzer {
	
	public FavaNetworkAnalyzer() {
	}
	
	/**
	 * Use this method to find a threshold for the FAVA network. The idea is that we want to get almost the same
	 * size of network as we are using for our interacting pathway analysis.
	 * @throws IOException
	 */
	@Test
	public void selectFAVAThreshold() throws IOException {
		Map<String, Map<String, Double>> gene2gene2score = loadFAVANetwork();
		System.out.println("Total genes: " + gene2gene2score.size());
		System.out.println("Threshold\tTotal_FIs");
		for (int i = 0; i < 20; i++) {
			double threshold = i / 20.0d;
			Set<String> fis = generateNetwork(gene2gene2score, threshold);
			System.out.println(threshold + "\t" + fis.size());
		}
		// Check with Reactome's random forest trained FIs we have used in our impact analysis
		String threshold = ApplicationConfig.getConfig().getAppConfig("score.threshold.for.bn");
		Map<String, Map<String, Double>> idgGene2gene2score = new PredictionScoreAnalyzer().loadPredictionFIFile(Double.parseDouble(threshold));
		Set<String> idgFIs = generateNetwork(idgGene2gene2score, Double.parseDouble(threshold));
		System.out.println("Random forest predicted fis with threshold " + threshold + ": " + idgFIs.size());
		System.out.println("The above FIs are used in our interacting pathway analysis.");
	}
	
	private Set<String> generateNetwork(Map<String, Map<String, Double>> gene2gene2score, double threshold) {
		Set<String> fis = new HashSet<>();
		gene2gene2score.keySet().stream().forEach(gene -> {
			Map<String, Double> gene2score = gene2gene2score.get(gene);
			for (String gene1 : gene2score.keySet()) {
				Double score = gene2score.get(gene1);
				if (score > threshold)
					fis.add(InteractionUtilities.generateFIFromGene(gene, gene1));
			}
		});
		return fis;
	}
	
	/**
	 * Load the original FAVA file.
	 * @return
	 * @throws IOException
	 */
	public Map<String, Map<String, Double>> loadFAVANetwork() throws IOException {
		String fileName = ApplicationConfig.getConfig().getAppConfig("fava.network.file");
		FileUtility fu = new FileUtility();
		fu.setInput(fileName);
		String line = fu.readLine();
		Map<String, Map<String, Double>> gene2gene2score = new HashMap<>();
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split(",");
			Map<String, Double> gene2score = gene2gene2score.get(tokens[0]);
			if (gene2score == null) {
				gene2score = new HashMap<>();
				gene2gene2score.put(tokens[0], gene2score);
			}
			gene2score.put(tokens[1], Double.parseDouble(tokens[2]));
		}
		fu.close();
		return gene2gene2score;
	}

	
}
