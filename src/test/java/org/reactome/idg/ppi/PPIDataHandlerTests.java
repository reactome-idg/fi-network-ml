package org.reactome.idg.ppi;

import java.util.Set;

import org.apache.log4j.Logger;

public abstract class PPIDataHandlerTests {
    private static final Logger logger = Logger.getLogger(PPIDataHandlerTests.class);
    
    protected void peekPPIs(Set<String> ppis) {
//        if (true)
//            return;
        logger.info("Printing out 10 PPIs:");
        int c = 0;
        for (String ppi : ppis) {
            logger.info(ppi);
            c ++;
            if (c == 10)
                break;
        }
    }

}
