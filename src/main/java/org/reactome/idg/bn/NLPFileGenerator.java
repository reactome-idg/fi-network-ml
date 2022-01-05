package org.reactome.idg.bn;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.util.FileUtilities;
import org.junit.Test;
import org.reactome.idg.util.DatabaseConfig;

/**
 * This class is used to generate some files used for NLP related projects.
 * @author wug
 *
 */
@SuppressWarnings("unchecked")
public class NLPFileGenerator {
	
	public NLPFileGenerator() {
	}
	
	/**
	 * This method is used to dump summary texts for each analyzed pathway so that we can perform
	 * topic modeling.
	 * @throws Exception
	 */
	@Test
	public void dumpSummaryTextForNLP() throws Exception {
		String outDir = ImpactResultAnalyzer.RESULT_DIR + "nlp_files/";
		File dir = new File(outDir);
		if (!dir.exists())
			dir.mkdir();
		// Dump the pathways used for simulation and interaction analysis
		MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
		PathwayImpactAnalyzer impactAnalyzer = new PathwayImpactAnalyzer();
		List<GKInstance> pathways = impactAnalyzer.loadPathwaysForAnalysis(dba);
		InstanceUtilities.sortInstances(pathways);
		FileUtilities fu = new FileUtilities();
//		String fileName = outDir + "PathwayText_111921.txt";
		String fileName = outDir + "PathwayText_120721.txt";
		fu.setOutput(fileName);
		for (GKInstance pathway : pathways) {
			List<String> textList = grepSummationTexts(pathway);
			// Pathway name starts with ###
			fu.printLine("###" + pathway.getDBID() + "||" + pathway.getDisplayName());
			// Seperate each summation with "//"
			for (String text : textList) {
				fu.printLine(text);
				fu.printLine("//");
			}
		}
		fu.close();
	}
	
	private List<String> grepSummationTexts(GKInstance pathway) throws Exception {
		Set<GKInstance> pathwayEvents = InstanceUtilities.getContainedEvents(pathway);
		pathwayEvents.add(pathway);
		List<String> textList = new ArrayList<>();
		for (GKInstance event : pathwayEvents) {
			List<GKInstance> list = event.getAttributeValuesList(ReactomeJavaConstants.summation);
			if (list == null || list.size() == 0)
				continue;
			for (GKInstance summation : list) {
				String text = (String) summation.getAttributeValue(ReactomeJavaConstants.text);
				if (text == null || text.trim().length() == 0)
					continue;
				// See https://colab.research.google.com/drive/12hfBveGHRsxhPIUMmJYrll2lFU4fOX06 for
				// using [SEP], most likely a sentence marker
				text = event.getDisplayName() + "[SEP]" + text;
				textList.add(text);
			}
		}
		return textList;
	}

	/**
	 * Files generated in this method are used for performing some MeSH term based NLP analysis via collaborating
	 * with Augustin Luna at Dana Faber.
	 * @throws IOException
	 */
	@Test
	public void dumpFilesForNLP() throws Exception {
		String outDir = ImpactResultAnalyzer.RESULT_DIR + "nlp_files/";
		File dir = new File(outDir);
		if (!dir.exists())
			dir.mkdir();
		FileUtilities fu = new FileUtilities();
		// Dump the pathways used for simulation and interaction analysis
		MySQLAdaptor dba = DatabaseConfig.getMySQLDBA();
		PathwayImpactAnalyzer impactAnalyzer = new PathwayImpactAnalyzer();
		List<GKInstance> pathways = impactAnalyzer.loadPathwaysForAnalysis(dba);
		InstanceUtilities.sortInstances(pathways);
	    dba.loadInstanceAttributeValues(pathways, new String[]{ReactomeJavaConstants.stableIdentifier});
//		String fileName = outDir + "AnalyzedPathways_111512.txt";
//		fu.setOutput(fileName);
//		fu.printLine("Pathway_Name\tStableId\tDB_ID");
//		for (GKInstance pathway : pathways) {
//			String id = getStableId(pathway);
//			fu.printLine(pathway.getDisplayName() + "\t" + 
//					     id + "\t" + 
//					     pathway.getDBID());
//		}
//		fu.close();
		// Dump the associated PMIDs for the above pathways
		String fileName = outDir + "PathwayId2PMID_111512.txt";
		fu.setOutput(fileName);
		dumpPMIDs(fu, pathways, dba);
		// Dump the dark proteins in the IDG three protein families
//		Map<String, String> target2levels = loadProtein2DevLevel();
//		Map<String, String> target2family = loadProtein2Family();
//		Set<String> selectedFamilies = Stream.of("IC", "Kinase", "GPCR").collect(Collectors.toSet());
//		List<String> targets = target2levels.keySet()
//				.stream()
//				.filter(t -> target2levels.get(t).equals("Tdark"))
//				.filter(t -> selectedFamilies.contains(target2family.get(t)))
//				.sorted()
//				.collect(Collectors.toList());
//		String fileName = outDir + "DarkProteins_3_families_111512.txt";
//		dumpTargets(targets, target2levels, target2family, fileName, fu);
//		// Randomly sample the same number of non-dark proteins 
//		List<String> sampledTargets = target2levels.keySet()
//				.stream()
//				.filter(t -> !target2levels.get(t).equals("Tdark"))
//				.collect(Collectors.toList());
//		Set<String> set = MathUtilities.randomSampling(sampledTargets, targets.size());
//		fileName = outDir + "NonDarkProteins_Sampled_111521.txt";
//		dumpTargets(set.stream().sorted().collect(Collectors.toList()),
//				target2levels,
//				target2family,
//				fileName, 
//				fu);
//		// Randomly sample the same number of dark proteins
//		sampledTargets = target2levels.keySet()
//				.stream()
//				.filter(t -> target2levels.get(t).equals("Tdark"))
//				.collect(Collectors.toList());
//		set = MathUtilities.randomSampling(sampledTargets, targets.size());
//		fileName = outDir + "DarkProteins_Sampled_111521.txt";
//		dumpTargets(set.stream().sorted().collect(Collectors.toList()),
//				target2levels,
//				target2family,
//				fileName, 
//				fu);
	}

	private void dumpPMIDs(FileUtilities fu,
	                       List<GKInstance> pathways,
	                       MySQLAdaptor dba) throws Exception {
		fu.printLine("Pathway_Stable_ID\tPMID");
		Set<Integer> totalPMIDs = new HashSet<>();
		for (GKInstance pathway : pathways) {
			List<Integer> pmids = loadPathwayPMIDs(pathway);
			totalPMIDs.addAll(pmids);
			String stableId = getStableId(pathway);
			for (Integer pmid : pmids)
				fu.printLine(stableId + "\t" + pmid);
		}
		fu.close();
		System.out.println("Total PMIDs: " + totalPMIDs.size());
		// Check PMIDs not in the pathway list
		Collection<GKInstance> references = dba.fetchInstancesByClass(ReactomeJavaConstants.LiteratureReference);
		dba.loadInstanceAttributeValues(references, new String[] {ReactomeJavaConstants.pubMedIdentifier});
		int c = 0;
		for (GKInstance reference : references) {
			Integer pmid = (Integer) reference.getAttributeValue(ReactomeJavaConstants.pubMedIdentifier);
			if (totalPMIDs.contains(pmid))
				continue;
			System.out.println(reference + ":");
			Collection<GKInstance> referrers = reference.getReferers(ReactomeJavaConstants.literatureReference);
			referrers.forEach(r -> System.out.println("\t" + r));
			c++;
			if (c == 20)
				break;
		}
	}

	private void dumpTargets(List<String> targets,
	                         Map<String, String> target2levels,
	                         Map<String, String> target2family,
	                         String fileName,
	                         FileUtilities fu)
	        throws IOException {
		fu.setOutput(fileName);
		fu.printLine("Protein\tFamily\tDevLevel");
		for (String target : targets) {
			fu.printLine(target + "\t" + 
		                 target2family.get(target) + "\t" + 
					     target2levels.get(target));
		}
		fu.close();
		System.out.println("Total targets in " + fileName + ": " + targets.size());
	}

	private String getStableId(GKInstance pathway) throws InvalidAttributeException, Exception {
		GKInstance stableId = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
		String id = (String) stableId.getAttributeValue(ReactomeJavaConstants.identifier);
		return id;
	}
	
	private List<Integer> loadPathwayPMIDs(GKInstance pathway) throws Exception {
		// Grep all event instances
		Set<GKInstance> pathwayEvents = InstanceUtilities.getContainedEvents(pathway);
		pathwayEvents.add(pathway);
		Set<GKInstance> references = new HashSet<>();
		for (GKInstance event : pathwayEvents) {
			loadPathwayPMIDs(event, references);
		}
		// References in the inferred from should be fetched
		for (GKInstance event : pathwayEvents) {
			List<GKInstance> inferredFromEvents = event.getAttributeValuesList(ReactomeJavaConstants.inferredFrom);
			if (inferredFromEvents == null || inferredFromEvents.size() == 0)
				continue;
			for (GKInstance inferrdFromEvent : inferredFromEvents) {
				loadPathwayPMIDs(inferrdFromEvent, references);
			}
		}
		List<Integer> rtn = new ArrayList<>();
		for (GKInstance reference : references) {
			// Some references may not have pubMedId
			if (!reference.getSchemClass().isValidAttribute(ReactomeJavaConstants.pubMedIdentifier))
				continue;
			Integer pmid = (Integer) reference.getAttributeValue(ReactomeJavaConstants.pubMedIdentifier);
			if (pmid == null) {
				System.err.println(reference + " doesn't have a PMID!");
				continue;
			}
			rtn.add(pmid);
		}
		return rtn.stream().sorted().collect(Collectors.toList());
	}

	private void loadPathwayPMIDs(GKInstance event, Set<GKInstance> references)
	        throws InvalidAttributeException, Exception {
		if (event.getSchemClass().isValidAttribute(ReactomeJavaConstants.literatureReference)) {
			List<GKInstance> litRefs = event.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
			if (litRefs != null && litRefs.size() > 0)
				references.addAll(litRefs);
		}
		grepReferences(event, "catalystActivityReference", references);
		grepReferences(event, "regulationReference", references);
		grepReferences(event, "summation", references);
	}
	
	private void grepReferences(GKInstance event, 
	                            String attributeName,
	                            Set<GKInstance> references) throws Exception {
		if (event.getSchemClass().isValidAttribute(attributeName)) {
			List<GKInstance> values = event.getAttributeValuesList(attributeName);
			if (values != null) {
				for (GKInstance value : values) {
					if (!value.getSchemClass().isValidAttribute(ReactomeJavaConstants.literatureReference))
						continue;
					List<GKInstance> litRefs = value.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
					if (litRefs != null && litRefs.size() > 0)
						references.addAll(litRefs);
				}
			}
		}
	}
	
}
