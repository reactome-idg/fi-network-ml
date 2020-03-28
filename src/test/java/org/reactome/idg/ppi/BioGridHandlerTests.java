package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Set;

import org.junit.Test;

public class BioGridHandlerTests extends PPIDataHandlerTests {
    private BioGridHandler handler;
    
    public BioGridHandlerTests() {
        handler = new BioGridHandler();
    }
    
    @Test
    public void testGetMousePPIs() throws IOException {
        Set<String> ppis = handler.loadMousePPIs();
        System.out.println("Total mouse ppis: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void getGetFlyPPIs() throws IOException {
        Set<String> ppis = handler.loadFlyPPIs();
        System.out.println("Total fly PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void testGetWormPPIs() throws IOException {
        Set<String> ppis = handler.loadWormPPIs();
        System.out.println("Total worm PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void getGetYeastPPIs() throws IOException {
        Set<String> ppis = handler.loadYeastPPIs();
        System.out.println("Total yeast PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void testGetHumanPPIs() throws IOException {
        Set<String> ppis = handler.loadHumanPPIs();
        System.out.println("Total human PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void testGetHumanPPIsFromFissionYeast() throws IOException {
        Set<String> ppis = handler.loadFissionYeastPPIs();
        System.out.println("Total fission yeast PPIs: " + ppis.size());
        peekPPIs(ppis);
    }

}
