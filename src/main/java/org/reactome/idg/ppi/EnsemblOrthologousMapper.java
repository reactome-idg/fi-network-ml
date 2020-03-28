package org.reactome.idg.ppi;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is copied directly from EnsemlAnalyzer class in the FINetworkBuild project.
 * In order to make methods here work, a protein family should be created first as described
 * in the original class.
 * TODO: Copying is to reduce the dependency. Most likely, we should try to refactor these two
 * projects (this project and the FINetwork project) to avoid copying.
 * @author wug
 *
 */
public class EnsemblOrthologousMapper implements OrthologousMapper {

    public EnsemblOrthologousMapper() {
    }

    public Map<String, Set<String>> loadYeastIdToHumanUniProtMap() throws IOException {
        return loadToHumanMapInUniProt("559292"); // need to use this yeast type!!!
        //        return loadToHumanMapInUniProt("4932");
    }
    
    public Map<String, Set<String>> loadFissionYeastIdToHumanUniProtMap() throws IOException {
        return null;
    }

    public Map<String, Set<String>> loadWormIdToHumanUniProtMap() throws IOException {
        return loadToHumanMapInUniProt("6239");
    }

    public Map<String, Set<String>> loadFlyIdToHumanUniProtMap() throws IOException {
        return loadToHumanMapInUniProt("7227");
    }

    public Map<String, Set<String>> loadMouseIdToHumanUniProtMap() throws IOException {
        return loadToHumanMapInUniProt("10090");
    }
    
    private Map<String, Set<String>> loadProteinFamilies() throws IOException {
        String fileName = ApplicationConfig.getConfig().getAppConfig("ensebml.protein.family.file");
        Map<String, Set<String>> idToFamilies = new HashMap<>();
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            lines.map(line -> line.split("\t"))
                 .forEach(tokens -> {
                     idToFamilies.compute(tokens[0], (key, set) -> {
                         if (set == null)
                             set = new HashSet<>();
                         set.add(tokens[1]);
                         return set;
                     });
                 });
            return idToFamilies;
        }
    }

    private Map<String, Set<String>> loadToHumanMapInUniProt(String taxonId) throws IOException {
        Map<String, Set<String>> familyToProteins = loadProteinFamilies();
        Set<String> humanIds = new HashSet<String>();
        Set<String> otherIds = new HashSet<String>();
        // To be returned
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        for (String family : familyToProteins.keySet()) {
            Set<String> proteins = familyToProteins.get(family);
            humanIds.clear();
            otherIds.clear();
            splitIds(proteins, humanIds, otherIds, taxonId);
            for (String otherId : otherIds) {
                Set<String> humanSet = map.get(otherId);
                if (humanSet == null) {
                    humanSet = new HashSet<String>();
                    map.put(otherId, humanSet);
                }
                humanSet.addAll(humanIds);
            }
        }
        return map;
    }

    private void splitIds(Set<String> proteins,
                          Set<String> humanIds,
                          Set<String> otherIds,
                          String taxonId) {
        for (String protein : proteins) {
            if (protein.startsWith("9606:")) {
                // This is a human protein
                humanIds.add(protein.substring(5));
            }
            else if (protein.startsWith(taxonId)) {
                otherIds.add(protein.substring(taxonId.length() + 1)); // 1 for ":".
            }
        }
    }

}
