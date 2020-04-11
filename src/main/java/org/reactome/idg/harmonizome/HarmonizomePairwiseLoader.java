package org.reactome.idg.harmonizome;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to load pre-processed pairwise relationships. The files are quite large and
 * it may take time and space to load these files.
 * @author wug
 *
 */
public class HarmonizomePairwiseLoader {
    private final String DIR = ApplicationConfig.getConfig().getAppConfig("harmonizome.filtered.dir");
    
    public HarmonizomePairwiseLoader() {
    }
    
    public List<File> getPairwiseFiles() {
        return getPairwiseFiles(new File(DIR));
    }
    
    public List<File> getPairwiseFiles(File dir) {
        File[] files = dir.listFiles();
        return Arrays.asList(files).stream()
                     .filter(file -> file.getName().endsWith("_filtered.txt"))
                     .collect(Collectors.toList());
    }
    
    public Set<String> loadPairwises(File file) throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(file.getAbsolutePath()))) {
            return lines.map(line -> line.split("\t"))
                        .map(tokens -> tokens[0] + "\t" + tokens[1])
                        .collect(Collectors.toSet());
        }
    }

}
