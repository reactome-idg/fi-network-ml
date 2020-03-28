package org.reactome.idg.ppi;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.gk.util.FileUtilities;
import org.reactome.idg.util.ApplicationConfig;

public class StringDBHandler extends PPIDataHandler {
    private final String DIR = ApplicationConfig.getConfig().getAppConfig("stringdb.dir");
    private final String EXP_CHANNEL_NAME = "experiments";
    private final ApplicationConfig config = ApplicationConfig.getConfig();

    public StringDBHandler() {
    }
    
    public Set<String> loadHumanPPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.human.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.human.map");
        return loadModPPIs(ppiFile, mapFile);
    }
    
    public Set<String> loadFlyPPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.fly.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.fly.map");
        return loadModPPIs(ppiFile, mapFile);
    }
    
    public Set<String> loadYeastPPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.yeast.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.yeast.map");
        return loadModPPIs(ppiFile, mapFile);
    }
    
    public Set<String> loadWormPPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.worm.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.worm.map");
        return loadModPPIs(ppiFile, mapFile);
    }
    
    public Set<String> loadMousePPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.mouse.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.mouse.map");
        return loadModPPIs(ppiFile, mapFile);
    }

    private Set<String> loadModPPIs(String ppiFileName, String mapFileName) throws IOException {
        Set<String> mousePPIs = grepPPIsBasedOnChannel(ppiFileName, EXP_CHANNEL_NAME);
        Map<String, String> stringDBToUniprot = loadStringIdToUniProt(mapFileName);
        Set<String> rtn = new HashSet<>();
        for (String ppi : mousePPIs) {
            String[] tokens = ppi.split("\t");
            String id1 = stringDBToUniprot.get(tokens[0]);
            String id2 = stringDBToUniprot.get(tokens[1]);
            String mapped = getPPI(id1, id2);
            if (mapped != null)
                rtn.add(mapped);
        }
        return rtn;
    }
    
    private Map<String, String> loadStringIdToUniProt(String fileName) throws IOException {
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        Map<String, String> stringDBIDToUniProt = new HashMap<String, String>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue; // Escape comment line
            String[] tokens = line.split("\t");
//            System.out.println(line);
            String uniprot = tokens[1].split("\\|")[0];
            if (stringDBIDToUniProt.containsKey(tokens[2]))
                throw new IllegalStateException("Same StringDB id mapped to more than one UniProt: " + tokens[2]);
            stringDBIDToUniProt.put(tokens[2], uniprot);
        }
        fu.close();
        return stringDBIDToUniProt;
    }
    
    private Set<String> grepPPIsBasedOnChannel(String fileName, 
                                               String channel) throws IOException {
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get the index of the channel
        String[] tokens = line.split(" ");
        int channelIndex = -1;
        for (int i = 0; i < tokens.length; i++) {
            if (channel.equalsIgnoreCase(tokens[i])) {
                channelIndex = i;
                break;
            }
        }
        if (channelIndex == -1) 
            throw new IllegalArgumentException(channel + " is not defined in the data file: " + fileName);
        Set<String> ppis = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split(" ");
            int channelScore = new Integer(tokens[channelIndex]);
            if (channelScore > 0) {
                String ppi = getPPI(tokens[0], tokens[1]);
                if (ppi == null)
                    continue; 
                ppis.add(ppi);
            }
        }
        fu.close();
        return ppis;
    }
    
}
