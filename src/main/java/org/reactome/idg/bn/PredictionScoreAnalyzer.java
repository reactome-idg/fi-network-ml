package org.reactome.idg.bn;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.ApplicationConfig;

import tech.tablesaw.api.FloatColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;
import tech.tablesaw.api.TextColumn;
import tech.tablesaw.io.DataFrameWriter;
import tech.tablesaw.plotly.Plot;
import tech.tablesaw.plotly.api.Histogram;
import tech.tablesaw.plotly.api.HorizontalBarPlot;
import tech.tablesaw.plotly.api.LinePlot;
import tech.tablesaw.plotly.components.Layout;

/**
 * Some methods that are used to analyze the prediction results from trained NBC are in
 * class org.reactome.idg.fi.FeatureFileGeneratorTests in the test directory. Some refactoring
 * will be needed to group these methods together here. 
 * @author wug
 *
 */
public class PredictionScoreAnalyzer {

	public PredictionScoreAnalyzer() {
	}
	
	/**
	 * Plot the score and precision relationships.
	 * @throws IOException
	 */
	@Test
	public void plotScoreAndPrecision() throws IOException {
		String fileName = ApplicationConfig.getConfig().getAppConfig("precision.recall.file");
		System.out.println(fileName);
		Table table = Table.read().csv(fileName);
		System.out.println(table.columnNames());
		// Need to add a new column between thresholds and Precision
		FloatColumn preScoreCol = FloatColumn.create("Precision*Score", table.rowCount());
		for (int i = 0; i < table.rowCount(); i++) {
			Row row = table.row(i);
			Double preScore = row.getDouble("Precision") * row.getDouble("Thresholds");
			preScoreCol.set(i, preScore.floatValue());
		}
		// Have to add it later on after filling the data
		table.addColumns(preScoreCol);
		Plot.show(LinePlot.create("Precesion and Thresholds",
				table,
				"Thresholds", 
				"Precision*Score"));
		Plot.show(LinePlot.create("Precesion and Thresholds",
				table,
				"Thresholds", 
				"Precision"));
		// Output for better plot
		File file = new File(fileName);
		fileName = file.getName().split("\\.")[0];
		File output = new File(file.getParent(), fileName + "_090621.csv");
		table.write().csv(output);
	}
	
	@Test
	public void plotFeatures() throws IOException {
		String fileName = ApplicationConfig.getConfig().getAppConfig("feature.file");
		Map<String, Integer> posFeature2Counts = new HashMap<String, Integer>();
		Map<String, Integer> negFeature2Counts = new HashMap<String, Integer>();
		FileUtilities fu = new FileUtilities();
		fu.setInput(fileName);
		String line = fu.readLine();
		String[] headers = line.split(",");
		int totalLines = 0;
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split(",");
			Map<String, Integer> feature2Count;
			if (tokens[1].equals("0"))
				feature2Count = negFeature2Counts;
			else
				feature2Count = posFeature2Counts;
			for (int i = 2; i < headers.length; i++) {
				if (tokens[i].equals("1")) {
					feature2Count.compute(headers[i], (k, v) -> {
						if (v == null)
							v = 1;
						else
							v ++;
						return v;
					});
				}
			}
			totalLines ++;
//			if (totalLines == 10000)
//				break;
		}
		fu.close();
		System.out.println("Total lines: " + totalLines);
		System.out.println(posFeature2Counts);
		System.out.println(negFeature2Counts);
		// Use tablesaw to create plot
		Table table = Table.create();
		// Use IntColumn so that we can do transpose
		IntColumn typeCol = IntColumn.create("FI");
		table.addColumns(typeCol);
		for (int i = 2; i < headers.length; i++) {
			IntColumn featureCol = IntColumn.create(headers[i]);
			table.addColumns(featureCol);
		}
		Row posRow = table.appendRow();
		posRow.setInt("FI", 1);
		posFeature2Counts.forEach((f, c) -> posRow.setInt(f, c));
		Row negRow = table.appendRow();
		negRow.setInt("FI", 0);
		negFeature2Counts.forEach((f, c) -> negRow.setInt(f, c));
		Set<String> features = new HashSet<String>(posFeature2Counts.keySet());
		features.addAll(negFeature2Counts.keySet());
		String[] cols = features.stream().sorted().toArray(String[]::new);
		Plot.show(HorizontalBarPlot.create(
				"Feature Counts",
				table,
				"FI",
				Layout.BarMode.GROUP,
				cols));
		// Output the file
		DataFrameWriter writer = table.write();
		fileName = "results/features_check/feature_counts_090321.csv";
		System.out.println("Output the table to file: " + fileName);
		writer.csv(fileName);
	}
	
	@Test
	public void testSimplePlot() {
		Table table = Table.create();
		TextColumn typeCol = TextColumn.create("FI");
		table.addColumns(typeCol);
		Row posRow = table.appendRow();
		posRow.setText("FI", "1");
	}

	/**
	 * Load the prediction scores calculated based on the trained random forest.
	 * @return
	 * @throws IOException
	 */
	public Map<String, Map<String, Double>> loadPredictionFIFile() throws IOException {
		return loadPredictionFIFile(0.0d); // Load all
	}
	
	public Map<String, Map<String, Double>> loadPredictionFIFile(double threshold) throws IOException {
		String fileName = ApplicationConfig.getConfig().getAppConfig("prediction.result.file");
		try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
			Map<String, Map<String, Double>> gene2gene2score = new HashMap<String, Map<String,Double>>();
			lines.skip(1)
			.map(line ->  line.split(","))
			.forEach(tokens -> {
				// Gene pair
				String[] genes = tokens[0].split("_");
				Double score = new Double(tokens[4]);
				if (score < threshold)
					return; // Do nothing
				insertGene2Score(genes[0], genes[1], score, gene2gene2score);
				insertGene2Score(genes[1], genes[0], score, gene2gene2score);
			});
			return gene2gene2score;
		}
	}

	@Test
	public void analyzePredictionFile() throws IOException {
		Map<String, Map<String, Double>> gene2gene2score = loadPredictionFIFile();
		System.out.println("Size: " + gene2gene2score.size());
		DescriptiveStatistics stat = new DescriptiveStatistics();
		gene2gene2score.forEach((gene, gene2score) -> {
			stat.addValue(gene2score.size());
		});
		System.out.println(stat.toString());
		// How about use a threshold = 0.50
		stat.clear();
		double threshold = 0.50d;
		gene2gene2score.forEach((gene, gene2score) -> {
			int counter = gene2score.keySet()
					.stream()
					.filter(g -> gene2score.get(g) > threshold)
					.collect(Collectors.counting()).intValue();
			stat.addValue(counter);
		});
		System.out.println("\nUse a threshold = " + threshold);
		System.out.println(stat.toString());
		// We want to see the distribution of these scores for Reactome FIs only
		stat.clear();
		Set<String> fis = ApplicationConfig.getConfig().loadReactomeFIsInGenes();
		Map<String, Double> fi2score = new HashMap<String, Double>();
		gene2gene2score.forEach((gene, gene2score) -> {
			gene2score.keySet()
			.stream()
			.forEach(target -> {
				String fi = InteractionUtilities.generateFIFromGene(gene, target);
				if (fis.contains(fi)) {
					fi2score.put(fi, gene2score.get(target));
				}
			});
		});
		// for plot
		double[] values = new double[fi2score.size()];
		int i = 0;
		for (Double v : fi2score.values()) {
			values[i] = v;
			stat.addValue(v);
			i ++;
		}
		System.out.println("\nPredicted score distribution for Reactome FIs:");
		System.out.println(stat.toString());
		System.out.println("median: " + stat.getPercentile(50.0d));
		System.out.println("5 percentile: " + stat.getPercentile(5.0d));
		// Plot the score for visualization 
		Plot.show(Histogram.create("Distribution of Predicted Score of Reactome FIs", values));
		// Output the score file
		String scoreFileName = "results/feature_files/prediction/prd_probabilities_prediction_061820_reactome_fi.txt";
		PrintWriter pr = new PrintWriter(scoreFileName);
		pr.println("FI\tScore");
		fi2score.keySet().stream().sorted().forEach(fi -> {
			pr.println(fi + "\t" + fi2score.get(fi));
		});
		pr.close();
	}

	private void insertGene2Score(String gene,
			String target,
			Double score,
			Map<String, Map<String, Double>> gene2gene2score) {
		Map<String, Double> gene2score = gene2gene2score.get(gene);
		if (gene2score == null) {
			gene2score = new HashMap<String, Double>();
			gene2gene2score.put(gene, gene2score);
		}
		gene2score.put(target, score);
	}
}
