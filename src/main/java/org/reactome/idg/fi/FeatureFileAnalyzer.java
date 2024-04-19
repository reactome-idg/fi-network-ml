package org.reactome.idg.fi;

import java.io.IOException;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to check the generated feature files to get some statistics. 
 * @author wug
 *
 */
public class FeatureFileAnalyzer {
    
    public FeatureFileAnalyzer() {
        
    }
    
    @Test
    public void checkTrainFeatureFile() throws IOException {
        String fileName = "results/feature_files/training/feature_matrix_041720.csv";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        int lines = 0;
        int fis = 0;
        int goInFI = 0;
        int humanPPIInFI = 0;
        int yeastPPIInFI = 0;
        while ((line = fu.readLine()) != null) {
            lines ++;
            String[] tokens = line.split(",");
            if (tokens[1].equals("1")) {
                fis ++;
                if (tokens[2].equals("1"))
                    humanPPIInFI ++;
                if (tokens[6].equals("1"))
                    yeastPPIInFI ++;
                if (tokens[8].equals("1"))
                    goInFI ++;
            }
        }
        fu.close();
        System.out.printf("Total lines: %,d\n"
                + "Total FIs: %,d\n"
                + "Total GO sharing in FIs: %,d\n" 
                + "Total HumanPPI in FIs: %,d\n"
                + "Total YeastPPI in FIs: %,d", 
                lines, fis, goInFI, humanPPIInFI, yeastPPIInFI);
    }

}
