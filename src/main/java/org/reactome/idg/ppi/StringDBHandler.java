package org.reactome.idg.ppi;

import static org.hamcrest.CoreMatchers.nullValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.util.FileUtilities;
import org.reactome.idg.util.ApplicationConfig;

public class StringDBHandler extends PPIDataHandler {
    private final String EXP_CHANNEL_NAME = "experiments";
    private final ApplicationConfig config = ApplicationConfig.getConfig();
    private final String DIR = config.getAppConfig("stringdb.dir");
    private final String UNIPROT_AC = "UniProt_AC";
    
    public StringDBHandler() {
    }
    
    /**
     * Gene names are directly loaded into the set. So there is no need for mapping. However,
     * UniProt ids are used for non-human species.
     * The mapping file used should be 9606.protein.info.v12.0.txt (note the version in the name).
     */
    public Set<String> loadHumanPPIs() throws IOException {
        String ppiFile = DIR + File.separator + config.getAppConfig("stringdb.human.file");
        String mapFile = DIR + File.separator + config.getAppConfig("stringdb.human.map");
        return loadHumanPPs(ppiFile, mapFile);
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
        Map<String, Set<String>> stringDBToUniprot = loadStringIdToUniProt(mapFileName);
        Set<String> rtn = new HashSet<>();
        for (String ppi : mousePPIs) {
            String[] tokens = ppi.split("\t");
            Set<String> ids1 = stringDBToUniprot.get(tokens[0]);
            if (ids1 == null || ids1.size() == 0)
                continue;
            Set<String> ids2 = stringDBToUniprot.get(tokens[1]);
            if (ids2 == null || ids2.size() == 0)
                continue;
            // Permutate these two sets
            for (String id1 : ids1) {
                for (String id2 : ids2) {
                    String mapped = getPPI(id1, id2);
                    if (mapped != null) // There is a chance!
                        rtn.add(mapped);
                }
            }
        }
        return rtn;
    }
    
    private Set<String> loadHumanPPs(String ppiFileName, 
                                     String mapFileName) throws IOException {
        Set<String> humanPPIs = grepPPIsBasedOnChannel(ppiFileName, EXP_CHANNEL_NAME);
        // Create the mapping file
        Map<String, String> stringIdToHumanGene = Files.lines(Paths.get(mapFileName))
                                                       .skip(1) // title line
                                                       .map(line -> line.split("\t"))
                                                       .collect(Collectors.toMap(tokens -> tokens[0], 
                                                                                 tokens -> tokens[1]));
        Set<String> rtn = new HashSet<>();
        for (String ppi : humanPPIs) {
            String[] tokens = ppi.split("\t");
            String id1 = stringIdToHumanGene.get(tokens[0]);
            String id2 = stringIdToHumanGene.get(tokens[1]);
            String mapped = getPPI(id1, id2);
            if (mapped != null)
                rtn.add(mapped);
        }
        return rtn;
    }
    
    /**
     * The mapping file used is 4932.protein.aliases.v12.0.txt. The database key is UniProt_AC. One StringDB id
     * may be mapped to more than one UniProt accession number. However, some of them may be unreviewed. In this
     * implementation, we just load all and let the orthologous mapping to do filtering!
     * @param fileName
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadStringIdToUniProt(String fileName) throws IOException {
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        Map<String, Set<String>> stringDBIDToUniProt = new HashMap<>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue; // Escape comment line
            String[] tokens = line.split("\t");
            if (!tokens[2].equals(UNIPROT_AC))
                continue;
//            System.out.println(line);
            Set<String> set = stringDBIDToUniProt.get(tokens[0]);
            if (set == null) {
                set = new HashSet<>();
                stringDBIDToUniProt.put(tokens[0], set);
            }
            set.add(tokens[1]);
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
            int channelScore = Integer.parseInt(tokens[channelIndex]);
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
