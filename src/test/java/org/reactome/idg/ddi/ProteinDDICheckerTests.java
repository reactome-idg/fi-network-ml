package org.reactome.idg.ddi;

import java.io.IOException;
import java.util.Set;

import org.junit.Test;
import org.reactome.idg.misc.ProteinDDIChecker;
import org.reactome.idg.ppi.PPIDataHandlerTests;

public class ProteinDDICheckerTests extends PPIDataHandlerTests {
    
    private ProteinDDIChecker checker;
    
    public ProteinDDICheckerTests() {
        checker = new ProteinDDIChecker();
    }
    
    @Test
    public void testLoadGenePairsViaDDIs() throws IOException {
        Set<String> totalPairs = checker.loadGenePairsViaDDIs();
        peekPPIs(totalPairs);
    }

}
