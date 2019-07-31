package org.reactome.idg.harmonizome;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

/**
 * This class is modified from the Python script provided by the Harmonizome project:
 * https://amp.pharm.mssm.edu/Harmonizome/static/harmonizomedownloader.py
 * @author wug
 *
 */
public class DataDownloader {
    private final Logger logger =LogManager.getLogger(DataDownloader.class);
    private final String URL = "https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/%s/%s";
    private final String SELECTED_DATA_TYPE = "gene_similarity_matrix_cosine.txt.gz";

    public DataDownloader() {
    }
    
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.err.println("Usage java org.reactome.harmonizome.DataDownloader dirName");
        }
        DataDownloader loader = new DataDownloader();
        loader.downloadDatasets(args[0]);
    }
    
    @Test
    public void checkDownloadedFiles() throws IOException {
        String dirName = "datasets/Harmonizome/download/gene_similarity_matrix_cosine";
        File dir = new File(dirName);
        for (File file : dir.listFiles()) {
            Scanner scanner = new Scanner(file);
            String header = scanner.nextLine();
            System.out.println("Name: " + file.getName());
            System.out.println("Size: " + file.length() / (1024 * 1024) + " M");
            String[] tokens = header.split("\t");
            System.out.println("Tokens in the header: " + tokens.length);
            scanner.close();
            System.out.println();
        }
    }
    
    /**
     * The method is used to perform downloading. Currently only gene_similarity_matrix_cosine.txt.gz
     * files are downloaded. The downloaded files will be unzipped automatically.
     * @param dirName the directory to hold all downloaded files.
     * @throws Exception
     */
    public void downloadDatasets(String dirName) throws Exception {
        File dir = new File(dirName);
        if (!dir.exists()) {
            throw new IllegalArgumentException(dirName + " doesn't exist! Create the directory first!");
        }
        logger.info("Download directory: " + dirName);
        Map<String, String> datasetToPath = getDataSet2Path();
        List<String> datasets = loadReactomeIDGDatasets();
        int c = 0;
        for (String dataset : datasets) {
            String path = datasetToPath.get(dataset);
            if (path == null)
                throw new IllegalStateException("Cannot find path for " + dataset);
            if (path.equals("epigenomicsdnaaccessibility")) { // Nothing in this dataset
                logger.info("Escape epigenomicsdnaaccessibility. Nothing is there!");
                continue; 
            }
            c ++;
            download(path, dirName);
            // download 5 for test. For production, we will parse these files directly
            // into a database
            if (c == 5)
                break;
        }
        System.out.println("Total download: " + c); // We should get 66 data sets
    }
    
    public void download(String path, 
                         String dirName) throws Exception {
        String url = String.format(URL, path, SELECTED_DATA_TYPE);
        logger.info("Starting download: " + url);
        java.net.URL urlObj = new java.net.URL(url);
        InputStream is = urlObj.openStream();
//        if (true) {
//            is.close();
//            return;
//        }
        GZIPInputStream zis = new GZIPInputStream(is);
        // There should be only one entry in the file
        Scanner scanner = new Scanner(zis);
        PrintWriter pw = new PrintWriter(dirName + File.separator + path + ".txt");
        String line = null;
        while (scanner.hasNext()) {
           line = scanner.nextLine();
           pw.println(line);
        }
        pw.close();
        scanner.close();
        zis.close();
        is.close();
        logger.info("Done!");
    }
    
    public List<String> loadReactomeIDGDatasets() throws IOException {
        String fileName = "Harmonizome_datasets_annotations_062819.txt";
        InputStream is = getClass().getClassLoader().getResourceAsStream(fileName);
        Scanner scanner = new Scanner(is);
        List<String> rtn = new ArrayList<>();
        String line = scanner.nextLine(); // Escape the first head line
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();
            String[] tokens = line.split("\t");
            if (tokens[13].equalsIgnoreCase("false"))
                continue;
            // The actual data set name is the combination of the first two tokens
            if (tokens[1].length() > 0) // We need to use shortname is it is available
                rtn.add(tokens[1] + " " + tokens[2]);
            else {
                String db = tokens[0];
                if (db.startsWith("\"")) {
                    db = db.substring(1, db.length() - 1);
                }
                rtn.add(db + " " + tokens[2]);
            }
        }
        scanner.close();
        is.close();
        return rtn;
    }
    
    public Map<String, String> getDataSet2Path() throws Exception {
        // This file was generated by copying the python code:
        // https://amp.pharm.mssm.edu/Harmonizome/static/harmonizomedownloader.py
        // We need to get the path name for each data set
        String fileName = "Harmonizome_data_mapping.txt";
        // The file will be copied together with classes in the top-most level.
        // We need to use ClassLoader to avoid using the sub package
        InputStream is = getClass().getClassLoader().getResourceAsStream(fileName);
        Map<String, String> dataset2Path = new HashMap<>();
        Scanner scanner = new Scanner(is);
        String line = null;
        Pattern pattern = Pattern.compile("'([^']+)'");
        int index = 0;
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();
            if (line.startsWith("##"))
                continue; // Escape some comments
            Matcher matcher = pattern.matcher(line);
            index = 0;
            String dataset = null;
            String path = null;
            while (matcher.find(index)) {
                // Group 1 gives us the actual value without ''
                String token = matcher.group(1);
                if (index == 0)
                    dataset = token;
                else
                    path = token;
                index = matcher.end();
            }
//            System.out.println(line);
//            System.out.println(dataset + " -> " + path);
            dataset2Path.put(dataset, path);
        }
        scanner.close();
        is.close();
        // Some datasets are not convered by Python script
        dataset2Path.put("LINCS Kinativ Kinase Inhibitor Bioactivity Profiles", "kinativ");
        dataset2Path.put("LINCS KinomeScan Kinase Inhibitor Targets", "kinomescan");
        return dataset2Path;
    }
    
    @Test
    public void testGetDataSet2Path() throws Exception {
        Map<String, String> dataset2Path = getDataSet2Path();
        System.out.println("DataSet2Path: " + dataset2Path.size());
        dataset2Path.forEach((data, path) -> System.out.println(data + " -> " + path));
    }

}
