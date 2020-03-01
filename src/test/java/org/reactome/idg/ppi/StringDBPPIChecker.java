package org.reactome.idg.ppi;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.gk.util.FileUtilities;
import org.junit.Test;

public class StringDBPPIChecker {
    private final String DIR_NAME = "/Users/wug/datasets/StringDB/";
    
    public StringDBPPIChecker() {
    }
    
    @Test
    public void checkNumbers() throws IOException {
        String fileName = DIR_NAME + "Human/9606.protein.links.full.v11.0.txt";
        int expCount = 0;
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = fu.readLine();
        int dbCountInExp = 0;
        int dbCount = 0;
        int lines = 0;
        Set<String> ppis = new HashSet<>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(" ");
            Integer expScore = Integer.parseInt(tokens[9]);
            Integer dbScore = Integer.parseInt(tokens[11]);
            if (expScore > 0) {
//                System.out.println(expScore);
                expCount ++;
                if (dbScore > 0)
                    dbCountInExp ++;
            }
            if (dbScore > 0) {
                dbCount ++;
            }
            lines ++;
            String ppi = getPPI(tokens[0], tokens[1]);
            if (ppi != null)
                ppis.add(ppi);
        }
        fu.close();
        System.out.println("Total PPIs having experimental score: " + expCount);
        System.out.println("\tIncluding DB source: " + dbCountInExp);
        System.out.println("Total PPIs having db score: " + dbCount);
        System.out.println("Total PPIs in the file: " + ppis.size());
        System.out.println("Total lines: " + lines);
    }
    
    @Test
    public void checkActionFile() throws IOException {
        String fileName = DIR_NAME + "Human/9606.protein.actions.v11.0.txt";
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> ppis = new HashSet<>();
        int lines = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String ppi = getPPI(tokens[0], tokens[1]);
            ppis.add(ppi);
            lines ++;
        }
        fu.close();
        System.out.println("Total PPIs: " + ppis.size());
        System.out.println("Total lines: " + lines);
    }
    
    @Test
    public void checkExperimentalBindingPPIs() throws IOException {
        Set<String> expPPIs = grepPPIsBasedOnChannel("database");
//        Set<String> expTransPPIs = grepPPIsBasedOnChannel("database_transferred");
//        expPPIs.addAll(expTransPPIs);
        System.out.println("Total experimental PPIs: " + expPPIs.size());
        Set<String> bindingPPIs = grepPPIsBasedOnAction("binding");
        System.out.println("Total binding PPIs: " + bindingPPIs.size());
        Set<String> shared = getShared(expPPIs, bindingPPIs);
        System.out.println("Experimental binding PPIs: " + shared.size());
        // Check the overlapping between exp PPIs and all annotated PPIs
        Set<String> allActionPPIs = grepPPIsBasedOnAction(null);
        System.out.println("Total action annotated PPIs: " + allActionPPIs.size());
        shared = getShared(expPPIs, allActionPPIs);
        System.out.println("Experimental action annotated PPIs: " + shared.size());
        double percentage = shared.size() / (double) expPPIs.size();
        System.out.println("Percentage: " + percentage);
    }
    
    private Set<String> getShared(Set<String> set1, Set<String> set2) {
        Set<String> shared = new HashSet<>(set1);
        shared.retainAll(set2);
        return shared;
    }
    
    private Set<String> grepPPIsBasedOnChannel(String channel) throws IOException {
        String fileName = DIR_NAME + "Human/9606.protein.links.full.v11.0.txt";
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
    
    private Set<String> grepPPIsBasedOnAction(String mode) throws IOException {
        String fileName = DIR_NAME + "Human/9606.protein.actions.v11.0.txt";
        Set<String> rtn = new HashSet<>();
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            rtn = lines.skip(1)
                       .map(line -> line.split("\t"))
                       .filter(tokens -> (mode == null ? true : tokens[2].equalsIgnoreCase(mode)))
                       .map(tokens -> getPPI(tokens[0], tokens[1]))
                       .filter(ppi -> ppi != null)
                       .collect(Collectors.toSet());
            return rtn;
        }
    }
    
    private String getPPI(String protein1, String protein2) {
        int compare = protein1.compareTo(protein2);
        if (compare > 0)
            return protein2 + "\t" + protein1;
        else if (compare < 0)
            return protein1 + "\t" + protein2;
        return null; // Avoid self-interaction
    }
    
    @Test
    public void checkOverlappingBetweenBioGridAndStringDB() throws IOException {
        Set<String> stringDbExpPPIsInGenes = loadExpPPIsInGenes();
        System.out.println("Total StringDB PPIs in genes: " + stringDbExpPPIsInGenes.size());
        Set<String> biogridPPIs = loadBioGridMVPPIs();
        System.out.println("Total BioGrid MV PPIs: " + biogridPPIs.size());
        Set<String> shared = getShared(stringDbExpPPIsInGenes, biogridPPIs);
        System.out.println("Total shared: " + shared.size());
        Set<String> biogridAllHumanPPIs = loadBioGridHumanPPIs();
        System.out.println("Total BioGrid human PPIs: " + biogridAllHumanPPIs.size());
        shared = getShared(stringDbExpPPIsInGenes, biogridAllHumanPPIs);
        System.out.println("Total shared: " + shared.size());
    }
    
    private Set<String> loadExpPPIsInGenes() throws IOException {
        Set<String> expPPIs = grepPPIsBasedOnChannel("experiments");
        System.out.println("Total PPIs annotated from experiments: " + expPPIs.size());
        Map<String, String> idToGene = loadStringDBIdToGene();
        // Map to genes
        Set<String> rtn = new HashSet<>();
        Set<String> notMapped = new HashSet<>();
        expPPIs.forEach(ppi -> {
            String[] tokens = ppi.split("\t");
            String gene1 = idToGene.get(tokens[0]);
            if (gene1 == null) {
                notMapped.add(tokens[0]);
                return;
//                throw new IllegalStateException(tokens[0] + " cannot be mapped to a gene!");
            }
            String gene2 = idToGene.get(tokens[1]);
            if (gene2 == null) {
                notMapped.add(tokens[1]);
                return;
            }
//                throw new IllegalStateException(tokens[1] + " cannot be mapped to a gene!");
            String geneInt = getPPI(gene1, gene2);
            if (geneInt != null)
                rtn.add(geneInt);
        });
        System.out.println("Total not mapped StringDB ids: " + notMapped.size());
        notMapped.stream().sorted().forEach(System.out::println);
        return rtn;
    }
    
    private Map<String, String> loadStringDBIdToGene() throws IOException {
        String fileName = DIR_NAME + "Human/human.name_2_string.tsv";
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            return lines.skip(1)
                 .map(line -> line.split("\t"))
                 .collect(Collectors.toMap(tokens -> tokens[2], tokens -> tokens[1]));
        }
    }
    
    @Test
    public void testLoadBioGridMVPPIs() throws IOException {
        Set<String> ppis = loadBioGridMVPPIs();
        System.out.println("Total ppis: " + ppis.size());
    }
    
    private Set<String> loadBioGridMVPPIs() throws IOException {
        String fileName = "datasets/BioGrid/BIOGRID-MV-Physical-3.5.181.tab2.txt";
        return loadBioGridPPIs(fileName, "9606");
    }
    
    private Set<String> loadBioGridHumanPPIs() throws IOException {
        String fileName = "datasets/BioGrid/BIOGRID-ORGANISM-3/BIOGRID-ORGANISM-Homo_sapiens-3.5.181.tab2.txt";
        return loadBioGridPPIs(fileName, "9606");
    }

    private Set<String> loadBioGridPPIs(String fileName, String species) throws IOException {
        FileUtilities fu = new FileUtilities();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> rtn = new HashSet<>();
        List<String> list = new ArrayList<>();
        int lines = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (species != null && !tokens[15].equals(species) || !tokens[16].equals(species))
                continue;
            lines ++;
            String ppi = getPPI(tokens[7], tokens[8]);
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

}
