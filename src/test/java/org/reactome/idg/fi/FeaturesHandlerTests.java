package org.reactome.idg.fi;

import org.junit.Test;

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

}
