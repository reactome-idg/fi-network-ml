package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Set;

public abstract class PPIDataHandler {

    protected String getPPI(String protein1, String protein2) {
        if (protein1 == null || protein2 == null)
            return null;
        int compare = protein1.compareTo(protein2);
        if (compare > 0)
            return protein2 + "\t" + protein1;
        else if (compare < 0)
            return protein1 + "\t" + protein2;
        return null; // Avoid self-interaction
    }
    
    /**
     * Load mouse PPIs should be in UniProt identifiers.
     * @return
     * @throws IOException
     */
    public abstract Set<String> loadMousePPIs() throws IOException;
    
    public abstract Set<String> loadFlyPPIs() throws IOException;
    
    public abstract Set<String> loadYeastPPIs() throws IOException;
    
    public abstract Set<String> loadWormPPIs() throws IOException;
    
    /**
     * Human PPIs are loaded in Human gene names.
     * @return
     * @throws IOException
     */
    public abstract Set<String> loadHumanPPIs() throws IOException;
    
}
