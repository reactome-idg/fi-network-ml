package org.reactome.idg.ppi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to handle import of BioGrid PPI data.
 * @author wug
 *
 */
public class BioGridHandler extends PPIDataHandler {
    private final String BIOGRID_DIR = ApplicationConfig.getConfig().getAppConfig("biogrid.dir");
    private ApplicationConfig config = ApplicationConfig.getConfig();
    
    public BioGridHandler() {
    }
    
    public Set<String> loadHumanPPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + 
                          config.getAppConfig("biogrid.human.name");
        String speciesId = config.getAppConfig("biogrid.human.id");
        return loadBioGridPPIs(fileName, speciesId, 7, 8);
    }
    
    public Set<String> loadMousePPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + config.getAppConfig("biogrid.mouse.name");
        String speciesId = config.getAppConfig("biogrid.mouse.id");
        String speciesDB = config.getAppConfig("biogrid.mouse.db");
        return loadMODPPIs(fileName, speciesId, speciesDB);
    }
    
    public Set<String> loadWormPPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + config.getAppConfig("biogrid.worm.name");
        String speciesId = config.getAppConfig("biogrid.worm.id");
        String speciesDB = config.getAppConfig("biogrid.worm.db");
        return loadMODPPIs(fileName, speciesId, speciesDB);
    }
    
    public Set<String> loadFlyPPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + config.getAppConfig("biogrid.fly.name");
        String speciesId = config.getAppConfig("biogrid.fly.id");
        String speciesDB = config.getAppConfig("biogrid.fly.db");
        return loadMODPPIs(fileName, speciesId, speciesDB);
    }
    
    public Set<String> loadYeastPPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + config.getAppConfig("biogrid.yeast.name");
        String speciesId = config.getAppConfig("biogrid.yeast.id");
        String speciesDB = config.getAppConfig("biogrid.yeast.db");
        return loadMODPPIs(fileName, speciesId, speciesDB);
    }
    
    public Set<String> loadFissionYeastPPIs() throws IOException {
        String fileName = BIOGRID_DIR + File.separator + config.getAppConfig("biogrid.fission.yeast.name");
        String speciesId = config.getAppConfig("biogrid.fission.yeast.id");
        String speciesDB = config.getAppConfig("biogrid.fission.yeast.db");
        return loadMODPPIs(fileName, speciesId, speciesDB);
    }

    private Set<String> loadMODPPIs(String fileName, String speciesId, String speciesDB) throws IOException {
        Set<String> ppis = loadBioGridPPIs(fileName,
                                           speciesId,
                                           3, // For yeast, we want to use BioGrid ids for easy mapping to human
                                           4);
        Map<String, String> biogridIdToOtherId = loadBioGridIdToOther(speciesDB);
        // Map PPIs to SGD ids
        Set<String> ppisInSGD = ppis.stream().map(ppi -> {
            String[] tokens = ppi.split("\t");
            String id1 = biogridIdToOtherId.get(tokens[0]);
            String id2 = biogridIdToOtherId.get(tokens[1]);
            return getPPI(id1, id2);
        }).filter(ppi -> ppi !=null).collect(Collectors.toSet());
        return ppisInSGD;
    }

    private Map<String, String> loadBioGridIdToOther(String otherType) throws IOException {
        String fileName = config.getAppConfig("biogrid.id.file.selected");
        boolean isInData = false;
        Map<String, String> bIdToOther = new HashMap<>();
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (isInData) {
                String[] tokens = line.split("\t");
                if (tokens[2].equals(otherType)) {
                    bIdToOther.put(tokens[0], tokens[1]);
                }
            }
            else if (line.startsWith("BIOGRID_ID")) {
                isInData = true;
            }
        }
        fu.close();
        return bIdToOther;
    }

    private Set<String> loadBioGridPPIs(String fileName, 
                                        String speciesId,
                                        int interactorIndex1,
                                        int interactorIndex2) throws IOException {
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> rtn = new HashSet<>();
        List<String> list = new ArrayList<>();
        int lines = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (speciesId != null && !tokens[15].equals(speciesId) || !tokens[16].equals(speciesId))
                continue;
            lines ++;
            String ppi = getPPI(tokens[interactorIndex1], 
                                tokens[interactorIndex2]);
            if (ppi != null) {
                //                if (rtn.contains(ppi)) {
                //                    System.out.println("Duplicated: " + ppi);
                //                }
                rtn.add(ppi);
                list.add(ppi);
            }
        }
        fu.close();
        //        System.out.println("Total reported human PPIs: " + lines);
        //        System.out.println("PPIs in the returned set: " + rtn.size());
        //        System.out.println("PPIs in the list: " + list.size());
        return rtn;
    }

    /**
     * This file is used to shorten the identifier files to a list of species to be used
     * in the IDG project.
     * @throws Exception
     */
    @Test
    public void processIdentifierFiles() throws Exception {
        String src = ApplicationConfig.getConfig().getAppConfig("biogrid.id.file");
        String dest = ApplicationConfig.getConfig().getAppConfig("biogrid.id.file.selected");
        String selectedSpecies = ApplicationConfig.getConfig().getAppConfig("biogrid.species");
        Set<String> speciesSet = Arrays.asList(selectedSpecies.split(",")).stream().collect(Collectors.toSet());
//        System.out.println("Selected species: " + speciesSet);
        FileUtilities fu = new FileUtilities();
        fu.setInput(src);
        fu.setOutput(dest);
        boolean isInData = false;
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (!isInData) {
                fu.printLine(line);
                if (line.startsWith("BIOGRID_ID"))
                    isInData = true;
            }
            else {
                String[] tokens = line.split("\t");
                if (speciesSet.contains(tokens[3]))
                    fu.printLine(line);
            }
        }
        fu.close();
    }

}
