package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Set;

import org.reactome.fi.util.InteractionUtilities;

public abstract class PPIDataHandler {

    protected String getPPI(String gene1, String gene2) {
        // In case either of genes is null
        if (gene1 == null || gene2 == null || gene1.equals(gene2))
            return null;
        return InteractionUtilities.generateFIFromGene(gene1, gene2);
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
