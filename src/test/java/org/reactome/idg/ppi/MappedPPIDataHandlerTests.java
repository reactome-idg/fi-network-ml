package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Set;

import org.apache.log4j.Logger;
import org.junit.Test;

public class MappedPPIDataHandlerTests extends PPIDataHandlerTests {
    private final static Logger logger = Logger.getLogger(MappedPPIDataHandlerTests.class);
    private MappedPPIDataHandler handler = new MappedPPIDataHandler();
    
    public MappedPPIDataHandlerTests() {
    }
    
    @Test
    public void testGetMousePPIs() throws IOException {
        Set<String> ppis = handler.loadMousePPIs();
        logger.info("Total mouse ppis: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void getGetFlyPPIs() throws IOException {
        Set<String> ppis = handler.loadFlyPPIs();
        logger.info("Total fly PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void testGetWormPPIs() throws IOException {
        Set<String> ppis = handler.loadWormPPIs();
        logger.info("Total worm PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void getGetYeastPPIs() throws IOException {
        Set<String> ppis = handler.loadYeastPPIs();
        logger.info("Total yeast PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
    @Test
    public void testGetHumanPPIs() throws IOException {
        Set<String> ppis = handler.loadHumanPPIs();
        logger.info("Total human PPIs: " + ppis.size());
        peekPPIs(ppis);
    }
    
}
