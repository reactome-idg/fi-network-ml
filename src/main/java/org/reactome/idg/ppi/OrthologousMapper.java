package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

public interface OrthologousMapper {
    
    public Map<String, Set<String>> loadYeastIdToHumanUniProtMap() throws IOException;
    
    public Map<String, Set<String>> loadFissionYeastIdToHumanUniProtMap() throws IOException;
    
    public Map<String, Set<String>> loadMouseIdToHumanUniProtMap() throws IOException;
    
    public Map<String, Set<String>> loadFlyIdToHumanUniProtMap() throws IOException;
    
    public Map<String, Set<String>> loadWormIdToHumanUniProtMap() throws IOException;


}
