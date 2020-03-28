package org.reactome.idg.ppi;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.util.FileUtilities;
import org.junit.Test;

public class OrthologousMapperChecker {
    
    @Test
    public void checkPatherYeastToHumanMap() throws Exception {
        String fileName = "/Users/wug/datasets/Panther/orthologs_14.1/HUMAN_YEAST.txt";
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        Set<String> yeastProteins = new HashSet<>();
        String line = null;
        String db = "UniProtKB";
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int index = tokens[1].indexOf(db);
            if (index < 0)
                throw new IllegalStateException(line + " doesn't have UniProt!");
            String id = tokens[1].substring(index + db.length() + 1);
            yeastProteins.add(id);
//            System.out.println(id);
        }
        fu.close();
        System.out.println("Total ids: " + yeastProteins.size());
        
        Map<String, Set<String>> yeastToHuman = loadYeastGeneToHumanMap();
        System.out.println("Total mapped yeast ids: " + yeastToHuman.size());
        Map<String, String> sgdIdToSysId = loadSGDIdToSysId();
        System.out.println("Total SGD id to Sys id: " + sgdIdToSysId.size());
    }
    
    @Test
    public void mapYeastPPIToHumanPPI() throws IOException {
        String dir = "/Users/wug/datasets/StringDB/Yeast/";
        String src = dir + "4932_PPIs_with_experiments.tsv";
        String target = dir + "human_4932_PPIs_with_experiments.tsv";
        
        Map<String, Set<String>> yeastSysIdToHumanIds = getYeastSysIdToHumanIds();
        
        FileUtilities fu = new FileUtilities();
        fu.setInput(src);
        // Push into a set to avoid any duplication
        Set<String> mappedHumanPPIs = new HashSet<>();
        String line = null;
        int totalUnMapped = 0;
        Set<String> unmappedYeastIds = new HashSet<>();
        int totalYeastPPIs = 0;
        while ((line = fu.readLine()) != null) {
            totalYeastPPIs ++;
            String[] tokens = line.split("\t");
            Set<String> humanIds1 = yeastSysIdToHumanIds.get(tokens[0]);
            if (humanIds1 == null) {
                totalUnMapped ++;
                unmappedYeastIds.add(tokens[0]);
                continue;
            }
            Set<String> humanIds2 = yeastSysIdToHumanIds.get(tokens[1]);
            if (humanIds2 == null) {
                totalUnMapped ++;
                unmappedYeastIds.add(tokens[1]);
                continue;
            }
            generatePPIs(humanIds1, humanIds2, mappedHumanPPIs);
        }
        fu.close();
        System.out.println("Total yeast PPIs: " + totalYeastPPIs);
        System.out.println("Total unmapped yeast PPIs: " + totalUnMapped);
        System.out.println("Total unmapped yeast ids: " + unmappedYeastIds.size());
        System.out.println("Total mapped human PPIs: " + mappedHumanPPIs.size());
        
        fu.setOutput(target);
        for (String humanPPI : mappedHumanPPIs.stream().sorted().collect(Collectors.toList())) {
            fu.printLine(humanPPI);
        }
        fu.close();
    }
    
    private void generatePPIs(Set<String> ids1, Set<String> ids2, Set<String> ppis) {
        for (String id1 : ids1) {
            for (String id2 : ids2) {
                int compare = id1.compareTo(id2);
                if (compare < 0)
                    ppis.add(id1 + "\t" + id2);
                else if (compare > 0)
                    ppis.add(id2 + "\t" + id1);
                // We don't want to have self interaction
            }
        }
    }
    
    
    private Map<String, Set<String>> getYeastSysIdToHumanIds() throws IOException {
        Map<String, String> sgdIdToSysId = loadSGDIdToSysId();
        Map<String, Set<String>> sgdIdToHumanIds = loadYeastGeneToHumanMap();
        Map<String, Set<String>> yeastSysIdToHumanIds = new HashMap<>();
        sgdIdToHumanIds.forEach((sgdId, humanSet) -> {
            String yeastSysId = sgdIdToSysId.get(sgdId);
            yeastSysIdToHumanIds.put(yeastSysId, humanSet);
        });
        return yeastSysIdToHumanIds;
    }
    
    private Map<String, String> loadSGDIdToSysId() throws IOException {
        String fileName = "/Users/wug/datasets/SGD/dbxref.tab";
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = null;
        Map<String, String> sgdToSys = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            sgdToSys.put(tokens[4], tokens[3]);
        }
        fu.close();
        return sgdToSys;
    }
    
    private Map<String, Set<String>> loadYeastGeneToHumanMap() throws IOException {
        String fileName = "/Users/wug/datasets/Panther/orthologs_14.1/HUMAN_YEAST.txt";
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = null;
        Map<String, Set<String>> yeastToHuman = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String humanId = extractId(tokens[0], 2);
            String yeastId = extractId(tokens[1], 1);
            yeastToHuman.compute(yeastId, (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(humanId);
                return set;
            });
        }
        fu.close();
        return yeastToHuman;
    }
    
    private String extractId(String token, int tokenIndex) {
        String[] tokens = token.split("\\|");
        if (tokens.length != 3)
            throw new IllegalStateException(token + " doesn't have three tokens!");
        int index = tokens[tokenIndex].indexOf("=");
        return tokens[tokenIndex].substring(index + 1);
    }
    
    @Test
    public void testLoadSGDIdToHumanUniProtMap() throws IOException {
        PantherOrthologousMapper mapper = new PantherOrthologousMapper();
        Map<String, Set<String>> yeastToHuman = mapper.loadYeastIdToHumanUniProtMap();
        System.out.println("Size: " + yeastToHuman.size());
        String yeastId = yeastToHuman.keySet().stream().findAny().get();
        Set<String> humanIds = yeastToHuman.get(yeastId);
        System.out.println("An example: " + yeastId + " -> " + humanIds);
    }
    
    @Test
    public void testParse() {
        String text = "YEAST|SGD=S000005930|UniProtKB=Q12532";
        String[] tokens = text.split("\\|");
        Arrays.asList(tokens).forEach(System.out::println);
    }

}
