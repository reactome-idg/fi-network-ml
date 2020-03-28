package org.reactome.idg.ppi;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.reactome.idg.util.ApplicationConfig;

/**
 * This class is used to map MOD PPIs to human.
 * @author wug
 *
 */
public class MappedPPIDataHandler extends PPIDataHandler {
    private static Logger logger = Logger.getLogger(MappedPPIDataHandler.class);
    private PPIDataHandler biogridHandler = new BioGridHandler();
    private PPIDataHandler stringDBHandler = new StringDBHandler();
//    private OrthologousMapper mapper = new EnsemblOrthologousMapper();
    private OrthologousMapper mapper = new PantherOrthologousMapper();
    private Map<String, String> uniprotToGene;
    private boolean useUniProt = false;
    
    public MappedPPIDataHandler() {
    }
    
    @Override
    public Set<String> loadHumanPPIs() throws IOException {
        Set<String> bPPIs = biogridHandler.loadHumanPPIs();
        logger.info("Total human PPIs from BioGrid: " + bPPIs.size());
        Set<String> sPPIs = stringDBHandler.loadHumanPPIs();
        logger.info("Total human PPIs from StringDB: " + sPPIs.size());
        return mergePPIs(bPPIs, sPPIs);
    }
    
    @Override
    public Set<String> loadFlyPPIs() throws IOException {
        Set<String> bPPIs = biogridHandler.loadFlyPPIs();
        logger.info("Total fly PPIs from BioGrid: " + bPPIs.size());
        Set<String> sPPIs = stringDBHandler.loadFlyPPIs();
        logger.info("Total fly PPIs from StringDB: " + sPPIs.size());
        Set<String> merged = mergePPIs(bPPIs, sPPIs);
        return mapMODPPIsToHuman(merged, mapper.loadFlyIdToHumanUniProtMap());
    }

    private Set<String> mergePPIs(Set<String> bPPIs, 
                                  Set<String> sPPIs) {
        Set<String> merged = new HashSet<>(bPPIs);
        merged.addAll(sPPIs);
        logger.info("Merged: " + merged.size());
        // For Debug purposes
        Set<String> shared = new HashSet<>(bPPIs);
        shared.retainAll(sPPIs);
        logger.info("Shared: " + shared.size());
        return merged;
    }

    @Override
    public Set<String> loadYeastPPIs() throws IOException {
        Set<String> bPPIs = biogridHandler.loadYeastPPIs();
        logger.info("Total yeast PPIs from BioGrid: " + bPPIs.size());
        Set<String> sPPIs = stringDBHandler.loadYeastPPIs();
        logger.info("Total yeast PPIs from StringDB: " + sPPIs.size());
        Set<String> merged = mergePPIs(bPPIs, sPPIs);
        return mapMODPPIsToHuman(merged, mapper.loadYeastIdToHumanUniProtMap());
    }

    @Override
    public Set<String> loadWormPPIs() throws IOException {
        Set<String> bPPIs = biogridHandler.loadWormPPIs();
        logger.info("Total worm PPIs from BioGrid: " + bPPIs.size());
        Set<String> sPPIs = stringDBHandler.loadWormPPIs();
        logger.info("Total worm PPIs from StringDB: " + sPPIs.size());
        Set<String> merged = mergePPIs(bPPIs, sPPIs);
        return mapMODPPIsToHuman(merged, mapper.loadWormIdToHumanUniProtMap());
    }

    public Set<String> loadMousePPIs() throws IOException {
        Set<String> bPPIs = biogridHandler.loadMousePPIs();
        logger.info("Total mouse PPIs from BioGrid: " + bPPIs.size());
        Set<String> sPPIs = stringDBHandler.loadMousePPIs();
        logger.info("Total mouse PPIs from StringDB: " + sPPIs.size());
        Set<String> merged = mergePPIs(bPPIs, sPPIs);
        // Special case for the mouse mapping, which generates a much better coverage
        // than the panther
        OrthologousMapper mapper = new EnsemblOrthologousMapper();
        return mapMODPPIsToHuman(merged,
                                 mapper.loadMouseIdToHumanUniProtMap());
    }
    
    private Set<String> mapMODPPIsToHuman(Set<String> modPPIs,
                                          Map<String, Set<String>> modIdToHumanUniProtMap) throws IOException {
        Set<String> humanPPIs = new HashSet<>();
        // We want to map to human genes directly
        Map<String, String> uniprotToGene = getUniProtToGene();
        for (String yeastPPI : modPPIs) {
            String[] sgdIds = yeastPPI.split("\t");
            Set<String> humanProt1 = modIdToHumanUniProtMap.get(sgdIds[0]);
            if (humanProt1 == null || humanProt1.size() == 0)
                continue;
            Set<String> humanProt2 = modIdToHumanUniProtMap.get(sgdIds[1]);
            if (humanProt2 == null || humanProt2.size() == 0)
                continue;
            for (String hId1 : humanProt1) {
                String gene1 = uniprotToGene.get(hId1);
                for (String hId2 : humanProt2) {
                    String gene2 = uniprotToGene.get(hId2);
                    String humanPPI = null;
                    if (useUniProt)
                        humanPPI = getPPI(hId1, hId2);
                    else
                        humanPPI = getPPI(gene1, gene2);
                    if (humanPPI != null)
                        humanPPIs.add(humanPPI);
                }
            }
        }
        return humanPPIs;
    }
    
    private Map<String, String> getUniProtToGene() throws IOException {
        if (uniprotToGene == null) {
            URL url = getClass().getClassLoader().getResource(ApplicationConfig.getConfig().getAppConfig("reactome.uniprot.to.gene"));
            uniprotToGene = Files.lines(Paths.get(url.getFile()))
                                 .skip(1) // The header line
                                 .map(line -> line.split("\t"))
                                 .collect(Collectors.toMap(tokens -> tokens[0], tokens -> tokens[1]));
        }
        return uniprotToGene;
    }


}
