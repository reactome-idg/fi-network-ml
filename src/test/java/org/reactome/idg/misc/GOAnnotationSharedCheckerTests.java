package org.reactome.idg.misc;

import java.io.IOException;
import java.util.Set;

import org.junit.Test;
import org.reactome.idg.misc.GOAnnotationShareChecker;
import org.reactome.idg.ppi.PPIDataHandlerTests;

public class GOAnnotationSharedCheckerTests extends PPIDataHandlerTests {
    
    private GOAnnotationShareChecker checker;
    
    public GOAnnotationSharedCheckerTests() {
        checker = new GOAnnotationShareChecker();
    }
    
    @Test
    public void testLoadGenePairsViaGOBPShare() throws IOException {
        Set<String> pairs = checker.loadGenePairsViaGOBPShare();
        peekPPIs(pairs);
    }

}
