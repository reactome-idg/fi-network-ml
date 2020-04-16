package org.reactome.idg.fi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.idg.util.ApplicationConfig;

public class FeaturesHandlerTests {
    private FeaturesCheckers handler;
    
    public FeaturesHandlerTests() {
        handler = new FeaturesCheckers();
    }
    
    @Test
    public void testCheckPPIFeatures() throws Exception {
//        handler.checkPPIFeatures();
//        handler.checkMiscFeatures();
//        handler.checkHarmonizomeFeatures();
        handler.checkGeneExpressionFeatures();
    }
    
    @Test
    public void checkMemoryUsage() throws Exception {
        Set<String> allGenes = ApplicationConfig.getConfig().getAllGenes();
        List<String> geneList = new ArrayList<>(allGenes);
        System.out.println("Total genes: " + allGenes.size());
        int[] sizes = new int[] {300000, 3000000, 30000000, 50000000};
        for (int size : sizes) {
            System.out.println("The total size of the list: " + size);
//            List<String> pairs = new ArrayList<>();
            Set<String> pairs = new HashSet<>();
            long memory = Runtime.getRuntime().totalMemory();
            System.out.println("Total memory before filling: " + memory / (1024 * 1024.0d) + " M");
            Random random = new Random();
            for (int i = 0; i < size; i++) {
                int index1 = random.nextInt(geneList.size());
                int index2 = random.nextInt(geneList.size());
                if (index1 == index2)
                    continue;
                pairs.add(InteractionUtilities.generateFIFromGene(geneList.get(index1),
                                                                  geneList.get(index2)));
            }
            memory = Runtime.getRuntime().totalMemory();
            System.out.println("Total memory after filling: " + memory / (1024 * 1024.0d) + " M");
            System.out.println("Size of pairs: " + pairs.size());
        }
    }

}
