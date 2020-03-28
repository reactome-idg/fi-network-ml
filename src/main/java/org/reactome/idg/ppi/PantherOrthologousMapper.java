package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.gk.util.FileUtilities;
import org.reactome.idg.util.ApplicationConfig;

public class PantherOrthologousMapper implements OrthologousMapper {
    private final String MAP_FILE = ApplicationConfig.getConfig().getAppConfig("panther.orthologous.map");
    private final FileUtilities fu = new FileUtilities();
    
    public PantherOrthologousMapper() {
    }
    
    public Map<String, Set<String>> loadYeastIdToHumanUniProtMap() throws IOException {
        String speciesName = ApplicationConfig.getConfig().getAppConfig("panther.yeast.name");
        return loadMODToHumanProtMap(speciesName, 2);
    }
    
    public Map<String, Set<String>> loadFissionYeastIdToHumanUniProtMap() throws IOException {
        String speciesName = ApplicationConfig.getConfig().getAppConfig("panther.fission.yeast.name");
        return loadMODToHumanProtMap(speciesName, 2);
    }
    
    public Map<String, Set<String>> loadMouseIdToHumanUniProtMap() throws IOException {
        String speciesName = ApplicationConfig.getConfig().getAppConfig("panther.mouse.name");
        return loadMODToHumanProtMap(speciesName, 2);
    }
    
    public Map<String, Set<String>> loadFlyIdToHumanUniProtMap() throws IOException {
        String speciesName = ApplicationConfig.getConfig().getAppConfig("panther.fly.name");
        return loadMODToHumanProtMap(speciesName, 2);
    }
    
    public Map<String, Set<String>> loadWormIdToHumanUniProtMap() throws IOException {
        String speciesName = ApplicationConfig.getConfig().getAppConfig("panther.worm.name");
//        String text = ApplicationConfig.getConfig().getAppConfig("panther.worm.index");
        return loadMODToHumanProtMap(speciesName, 2);
    }

    private Map<String, Set<String>> loadMODToHumanProtMap(String speciesName,
                                                           int modIdIndex) throws IOException {
        fu.setInput(MAP_FILE);
        String line = null;
        Map<String, Set<String>> yeastToHuman = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!tokens[1].startsWith(speciesName))
                continue;
            String humanId = extractId(tokens[0], 2);
            String yeastId = extractId(tokens[1], modIdIndex);
            yeastToHuman.compute(yeastId, (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(humanId);
                return set;
            });
        }
        fu.close();
        return yeastToHuman;
    }
    
    private String extractId(String token, int tokenIndex) {
        String[] tokens = token.split("\\|");
        if (tokens.length != 3)
            throw new IllegalStateException(token + " doesn't have three tokens!");
        int index = tokens[tokenIndex].indexOf("=");
        return tokens[tokenIndex].substring(index + 1);
    }

}
