package org.reactome.idg.harmonizome;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.annotations.FeatureDesc;
import org.reactome.idg.model.FeatureSource;
import org.reactome.idg.model.FeatureType;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to load pre-processed pairwise relationships. The files are quite large and
 * it may take time and space to load these files.
 * @author wug
 *
 */
public class HarmonizomePairwiseLoader {
    private final static Logger logger = Logger.getLogger(HarmonizomePairwiseLoader.class);
    private final String DIR = ApplicationConfig.getConfig().getAppConfig("harmonizome.filtered.dir");
    
    public HarmonizomePairwiseLoader() {
    }
    
    public List<File> getPairwiseFiles() {
        return getPairwiseFiles(new File(DIR));
    }
    
    public List<File> getPairwiseFiles(File dir) {
        return getFiles(dir, "_filtered.txt");
    }
    
    public Set<String> loadPairwises(File file) throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(file.getAbsolutePath()))) {
            return lines.map(line -> line.split("\t"))
                        .map(tokens -> tokens[0] + "\t" + tokens[1])
                        .collect(Collectors.toSet());
        }
    }
    
    public List<File> getProcessedFiles() {
        String dirName = ApplicationConfig.getConfig().getAppConfig("harmonizome.processed.dir");
        File file = new File(dirName);
        return getFiles(file, "_processed.txt");
    }
    
    public List<File> getDownloadedFiles() {
        String dirName = ApplicationConfig.getConfig().getAppConfig("harmonizome.downloaded.dir");
        File dir = new File(dirName);
        return getFiles(dir, ".txt.tgz");
    }
    
    public Map<File, Double> getSelectedDownloadFiles() throws IOException {
        String dirName = ApplicationConfig.getConfig().getAppConfig("harmonizome.downloaded.dir");
        String fileName = ApplicationConfig.getConfig().getAppConfig("harmonizome.selected.download.file");
        InputStream is = ApplicationConfig.getConfig().getInputStream(fileName);
        try (Stream<String> lines = new BufferedReader(new InputStreamReader(is)).lines()) {
            Map<File, Double> file2Percentile = lines.skip(1)
                                    .filter(line -> !line.startsWith("#")) // Filter these files for easy editing
                                    .map(line -> line.split("\t"))
                                    .collect(Collectors.toMap(tokens -> new File(dirName, tokens[0]),
                                                              tokens -> new Double(tokens[1])));
            return file2Percentile;
        }
    }
    
    private List<File> getFiles(File dir, String ext) {
        return Arrays.asList(dir.listFiles())
                .stream()
                .filter(file -> file.getName().endsWith(ext))
                .collect(Collectors.toList());
    }
    
    /**
     * Load a set of pairwise relationships directly from a downloaded harmonizome file.
     * @param file
     * @param percentile
     * @return
     * @throws IOException
     */
    @FeatureDesc(sources = {FeatureSource.Harmonizome},
                 type = FeatureType.Gene_Similarity)
    public Set<String> loadPairwisesFromDownload(File file, double percentile) throws Exception {
        File dir = file.getParentFile();
        logger.info("Handling " + file.getName() + "...");
        // Unzip the file first
        String fileName = file.getName();
        int index = fileName.lastIndexOf(".");
        fileName = fileName.substring(0, index);
        File unzipped = new File(dir, fileName);
        logger.info("Unzipping to " + fileName);
        new DataDownloader().unzipDownload(file, unzipped);
        logger.info("Processing the file ...");
        DataProcessor processor = new DataProcessor();
        // To avoid code duplication, we will generate two files: processed and filtered
        // and then delete them after the use
        // By resetting the following two values, we will process all files
        DataProcessor.GENE_NUMBER_THRESHOLD = 0;
        DataProcessor.VALUE_THRESHOLD = 0.0d;
        File processed = new File(dir, fileName + ".processed");
        File filtered = new File(dir, fileName + ".filtered");
        processor.processCorrelations(unzipped, processed, filtered);
        Set<String> rtn = loadPairwisesFromProcessed(processed, percentile);
        // Clean up all temporary files
        unzipped.delete();
        processed.delete();
        filtered.delete();
        return rtn;
    }
    
    /**
     * Load the pairwise relationships from the top percentile stored in a processed file.
     * @param file
     * @param percentile
     * @return
     * @throws IOException
     */
    public Set<String> loadPairwisesFromProcessed(File file, double percentile) throws IOException {
        logger.info("Loading all values for " + file.getName() + "...");
        // Need to find the threshold first
        List<Float> values = new ArrayList<>();
        FileUtility fu = new FileUtility();
        fu.setInput(file.getAbsolutePath());
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].length() == 0)
                    continue;
                Float value = new Float(tokens[i]);
                values.add(Math.abs(value));
            }
        }
        fu.close();
        values.sort(Comparator.reverseOrder());
        int index = (int) (values.size() * percentile);
        double cutoff = values.get(index);
        logger.info("Total values: " + values.size());
        logger.info("Cutoff value: " + cutoff);
        // If the cutoff value is 1.0, we should adjust it to a lower value
        // Otherwise, nothing will be got
        if (cutoff > ApplicationConfig.getConfig().getMaximumCutoff()) {
            cutoff = ApplicationConfig.getConfig().getMaximumCutoff();
            logger.info("Cutoff adjusted to: " + cutoff);
        }
        // Start to load the expression values based on the cutoff
        Set<String> rels = new HashSet<>();
        fu.setInput(file.getAbsolutePath());
        line = fu.readLine();
        String[] genes = line.split("\t");
        int lineIndex = 0;
        while ((line = fu.readLine()) != null) {
            lineIndex ++;
            String[] tokens = line.split("\t");
            for (int i = lineIndex + 1; i < tokens.length; i++) {
                Double value = new Double(tokens[i]);
                if (Math.abs(value) <= cutoff)
                    continue;
                String rel = InteractionUtilities.generateFIFromGene(tokens[0], genes[i]);
                if (rel != null)
                    rels.add(rel);
            }
        }
        fu.close();
        logger.info("Total collected relationships: " + rels.size());
        return rels;
    }

}
