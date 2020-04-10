package org.reactome.idg.ppi;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.reactome.fi.util.FileUtility;
import org.reactome.idg.util.ApplicationConfig;

public class BioPlexHandler extends PPIDataHandler {
    private final String DIR = ApplicationConfig.getConfig().getAppConfig("bioplex.dir") + File.separator;

    public BioPlexHandler() {
    }

    @Override
    public Set<String> loadMousePPIs() throws IOException {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Set<String> loadFlyPPIs() throws IOException {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Set<String> loadYeastPPIs() throws IOException {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Set<String> loadWormPPIs() throws IOException {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Set<String> loadHumanPPIs() throws IOException {
        Set<String> ppis = new HashSet<>();
        String fileName = DIR + ApplicationConfig.getConfig().getAppConfig("bioplex.human.293t.file");
        ppis.addAll(loadHumanPPIs(fileName));
        fileName = DIR + ApplicationConfig.getConfig().getAppConfig("bioplex.human.HTC116.file");
        ppis.addAll(loadHumanPPIs(fileName));
        return ppis;
    }
    
    private Set<String> loadHumanPPIs(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        Set<String> ppis = new HashSet<>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // We want to have gene names for PPIs
            String gene1 = stripQuotations(tokens[4]);
            String gene2 = stripQuotations(tokens[5]);
            String ppi = getPPI(gene1, gene2);
            ppis.add(ppi);
        }
        fu.close();
        return ppis;
    }
    
    private String stripQuotations(String token) {
        if (token.startsWith("\""))
            token = token.substring(1);
        if (token.endsWith("\""))
            token = token.substring(0, token.length() - 1);
        return token;
    }
    
}
