package org.reactome.idg.corr;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.persistence.MySQLAdaptor;
import org.gk.render.Node;
import org.gk.render.RenderablePathway;
import org.gk.schema.InvalidAttributeException;
import org.junit.Test;
import org.reactome.idg.coexpression.CoExpressionLoader;

@SuppressWarnings("unchecked")
public class TCGAGeneExpressionTest {
    private static final Logger logger = Logger.getLogger(TCGAGeneExpressionTest.class);
    
    public TCGAGeneExpressionTest() {
    }
    
    @Test
    public void testGetCutoffValueForRatio() throws IOException {
        CoExpressionLoader loader = new CoExpressionLoader();
        List<File> files = loader.getGTExCoExpressionFiles();
        double percentile = 0.001d;
        long time1 = System.currentTimeMillis();
        long memory = Runtime.getRuntime().totalMemory();
        logger.info("Total memory before loading: " + memory / (1024 * 1024 * 1024.0d) + " G.");
        for (File file : files) {
            logger.info("Loading file: " + file.getName());
            double cutoff = loader.getCutoffValueForRatio(file, percentile);
            logger.info("cutoff: " + cutoff);
            break;
        }
        long time2 = System.currentTimeMillis();
        logger.info("Total time for loading: " + (time2 - time1) / 1000.0d + " seconds.");
        memory = Runtime.getRuntime().totalMemory();
        logger.info("Total memory after loading: " + memory / (1024 * 1024 * 1024.0d) + " G.");
    }
    
    @Test
    public void testLoader() throws IOException {
        CoExpressionLoader loader = new CoExpressionLoader();
        List<File> files = loader.getGTExCoExpressionFiles();
        long time1 = System.currentTimeMillis();
        for (File file : files) {
            logger.info("Loading file: " + file.getName());
            Set<String> rels = loader.loadCoExpression(file, 0.50d);
            logger.info("Total rels: " + rels.size());
            break;
        }
        long time2 = System.currentTimeMillis();
        logger.info("Total time for loading: " + (time2 - time1) / 1000.0d + " seconds.");
    }
    
    public static void main(String[] args) throws Exception {
        new TCGAGeneExpressionTest().testGetCutoffValueForRatio();
    }
    
    @Test
    public void filterMCLClusterToELVPathways() throws Exception {
//        String dir = "/Users/wug/Desktop/";
//        
//        String input = "low_path_EA_pearson.csv";
//        String output = "low_path_EA_pearson_filtered.csv";
//        
//        input = "low_path_EA_bicor.csv";
//        output = "low_path_EA_bicor_filtered.csv";
//        
//        input = "low_path_EA_spearman.csv";
//        output = "low_path_EA_spearman_filtered.csv";
//        
//        MySQLAdaptor dba = getDBA();
//        Set<GKInstance> normalPathways = grepNormalPathways(dba);
//        Set<GKInstance> elvPathways = grepELVPathways(normalPathways, dba);
//        Set<String> stableIds = grepStableIds(elvPathways);
//        
//        FileUtilities fu = new FileUtilities();
//        fu.setInput(dir + input);
//        fu.setOutput(dir + output);
//        String line = fu.readLine();
//        fu.printLine(line); // This is the header
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split(",");
//            if (!stableIds.contains(tokens[5])) // Column 5 is for stable ids
//                continue;
//            fu.printLine(line);
//        }
//        fu.close();
    }
    
    private Set<String> grepStableIds(Set<GKInstance> pathways) throws Exception {
        Set<String> stableIds = new HashSet<>();
        for (GKInstance pathway : pathways) {
            GKInstance stableId = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
            String value = (String) stableId.getAttributeValue(ReactomeJavaConstants.identifier);
            if (value != null)
                stableIds.add(value);
        }
        System.out.println("Total stable ids: " + stableIds.size());
        stableIds.stream().sorted().forEach(System.out::println);
        return stableIds;
    }
    
    private Set<GKInstance> grepELVPathways(Set<GKInstance> pathways, MySQLAdaptor dba) throws Exception {
        Set<GKInstance> elvPathways = new HashSet<>();
        DiagramGKBReader reader = new DiagramGKBReader();
        for (GKInstance pathway : pathways) {
            Collection<GKInstance> diagrams = pathway.getReferers(ReactomeJavaConstants.representedPathway);
            if (diagrams == null || diagrams.size() == 0)
                continue;
            GKInstance diagram = diagrams.iterator().next();
            RenderablePathway rDiagram = reader.openDiagram(diagram);
            boolean hasELV = false;
            for (Object r : rDiagram.getComponents()) {
                if (r instanceof Node) {
                    Node node = (Node) r;
                    if (node.getReactomeId() != null) {
                        GKInstance inst = dba.fetchInstance(node.getReactomeId());
                        if (inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                            hasELV = true;
                            break;
                        }
                    }
                }
            }
            if (hasELV)
                elvPathways.add(pathway);
        }
        System.out.println("Total ELV pathways: " + elvPathways.size());
        return elvPathways;
    }

    private Set<GKInstance> grepNormalPathways(MySQLAdaptor dba) throws Exception, InvalidAttributeException {
        GKInstance frontPage = (GKInstance) dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage).iterator().next();
        List<GKInstance> frontPageItems = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        Set<GKInstance> pathways = new HashSet<>();
        for (GKInstance frontPageItem : frontPageItems) {
            if (frontPageItem.getDBID().equals(1643685L))
                continue;
            GKInstance species = (GKInstance) frontPageItem.getAttributeValue(ReactomeJavaConstants.species);
            if (!species.getDisplayName().equals("Homo sapiens"))
                continue;
            Set<GKInstance> subPathways = InstanceUtilities.getContainedInstances(frontPageItem, ReactomeJavaConstants.hasEvent);
            for (GKInstance subPathway : subPathways) {
                if (subPathway.getSchemClass().isa(ReactomeJavaConstants.Pathway))
                    pathways.add(subPathway);
            }
            pathways.add(frontPageItem);
        }
        System.out.println("Total normal pathways: " + pathways.size());
        return pathways;
    }
    
    private MySQLAdaptor getDBA() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "gk_current_ver70",
                                            "root",
                                            "macmysql01");
        return dba;
    }

}
