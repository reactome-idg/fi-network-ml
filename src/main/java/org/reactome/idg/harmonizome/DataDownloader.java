package org.reactome.idg.harmonizome;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
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
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;
import org.reactome.idg.model.Provenance;

/**
 * This class is modified from the Python script provided by the Harmonizome project:
 * https://amp.pharm.mssm.edu/Harmonizome/static/harmonizomedownloader.py
 * @author wug
 *
 */
public class DataDownloader {
    private final Logger logger = LogManager.getLogger(DataDownloader.class);
    public final String SOURCE = "Harmonizome";
    public final String URL = "https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/%s/%s";
    public final String SELECTED_DATA_TYPE = "gene_similarity_matrix_cosine.txt.gz";

    public DataDownloader() {
    }
    
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.err.println("Usage java org.reactome.harmonizome.DataDownloader dirName cleanupFolder{true|false}");
        }
        DataDownloader loader = new DataDownloader();
        loader.downloadDatasets(args[0]);
        if (args.length > 1 && args[1].equals("true"))
            loader.cleanupDir(args[0]);
    }
    
    /**
     * All original download files ending with .tgz will be deleted. All processed
     * files will be zipped and then deleted. 
     * @param dir
     * @throws Exception
     */
    private void cleanupDir(String dirName) throws Exception {
        logger.info("Cleaning up the download directory...");
        long time1 = System.currentTimeMillis();
        File dir = new File(dirName);
        File[] list = dir.listFiles();
        List<File> processedFiles = new ArrayList<>();
        // Delete all downloaded tgz file
        for (File file : list) {
            if (file.getName().endsWith("_processed.txt"))
                processedFiles.add(file);
            else if (file.getName().endsWith(".txt.tgz"))
                file.delete();
        }
        // Zip all processed files and then delete them.
        File dest = new File(dirName, "harmonizome_processed.zip");
        FileOutputStream fos = new FileOutputStream(dest);
        ZipOutputStream zos = new ZipOutputStream(new BufferedOutputStream(fos));
        int buffer = 100 * 1024;
        byte[] data = new byte[buffer];
        int read;
        for (File file : processedFiles) {
            FileInputStream fis = new FileInputStream(file);
            BufferedInputStream bis = new BufferedInputStream(fis, buffer);
            ZipEntry entry = new ZipEntry(file.getName());
            zos.putNextEntry(entry);
            while ((read = bis.read(data, 0, buffer)) > 0) {
                zos.write(data, 0, read);
            }
            bis.close();
            fis.close();
            file.delete();
        }
        zos.close();
        long time2 = System.currentTimeMillis();
        logger.info("Cleaning up is done: " + (time2 - time1) / (60.0d * 1000) + " minutes.");
    }
    
    /**
     * Get a list of Provenance objects for the Harmonizome project.
     * @return
     * @throws IOException
     */
    public List<Provenance> getProvenaces() throws Exception {
        List<Provenance> list = new ArrayList<>();
        Map<String, String> datasetToPath = getDataSet2Path();
        InputStream is = getDataResource();
        Scanner scanner = new Scanner(is);
        String line = scanner.nextLine(); // Escape the first head line
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();
            String[] tokens = line.split("\t");
            if (tokens[13].equalsIgnoreCase("false"))
                continue;
            String name = null;
            // The actual data set name is the combination of the first two tokens
            if (tokens[1].length() > 0) // We need to use shortname if it is available
                name = tokens[1] + " " + tokens[2];
            else {
                String db = tokens[0];
                if (db.startsWith("\"")) {
                    db = db.substring(1, db.length() - 1);
                }
                name = db + " " + tokens[2];
            }
            Provenance provenance = new Provenance();
            provenance.setName(name);
            provenance.setSource(SOURCE);
            provenance.setCategory(tokens[10]);
            provenance.setBiologicalEntity(tokens[8]);
            String path = datasetToPath.get(name);
            String url = String.format(URL, path, SELECTED_DATA_TYPE);
            provenance.setUrl(url);
            list.add(provenance);
        }
        scanner.close();
        is.close();
        return list;
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
        DataProcessor processor = new DataProcessor();
        for (String dataset : datasets) {
            String path = datasetToPath.get(dataset);
            if (path == null)
                throw new IllegalStateException("Cannot find path for " + dataset);
            if (path.equals("epigenomicsdnaaccessibility")) { // Nothing in this dataset
                logger.info("Escape epigenomicsdnaaccessibility. Nothing is there!");
                continue; 
            }
            c ++;
            String fileName = download(path, dirName);
            process(new File(fileName), dirName, processor);
            // download 2 for test. For production, we will parse these files directly
            // into a database
//            if (c == 2)
//                break;
        }
        System.out.println("Total download: " + c); // We should get 66 data sets
    }
    
    public String download(String path, 
                           String dirName) throws Exception {
        String url = String.format(URL, path, SELECTED_DATA_TYPE);
        logger.info("Starting download: " + url);
        long time1 = System.currentTimeMillis();
        java.net.URL urlObj = new java.net.URL(url);
        InputStream is = urlObj.openStream();
        // Try to download the zipped file first
        int size = 100 * 1024;  // Assign 100k for buffer
        BufferedInputStream bis = new BufferedInputStream(is, size);
        File gzFile = new File(dirName + "/" + path + ".txt.tgz");
        FileOutputStream fos = new FileOutputStream(gzFile);
        BufferedOutputStream bos = new BufferedOutputStream(fos);
        byte[] content = new byte[size];
        int read = bis.read(content, 0, size);
        while (read > 0) {
            bos.write(content, 0, read);
            read = bis.read(content, 0, size);
        }
        bos.close();
        bis.close();
        is.close();
        
//        if (true) {
//            is.close();
//            return;
//        }
        // Assign 100K as the buffer to increase the performance
        FileInputStream fis = new FileInputStream(gzFile);
        GZIPInputStream zis = new GZIPInputStream(fis, 100 * 1024);
        // There should be only one entry in the file
        Scanner scanner = new Scanner(zis);
        String fileName = dirName + "/" + path + ".txt";
        PrintWriter pw = new PrintWriter(fileName);
        String line = null;
        while (scanner.hasNext()) {
           line = scanner.nextLine();
           pw.println(line);
        }
        pw.close();
        scanner.close();
        zis.close();
        fis.close();
        long time2 = System.currentTimeMillis();
        logger.info("Total time for downloading: " + (time2 - time1) / (60.0d * 1000) + " minutes.");
        return fileName;
    }
    
    @Test
    public void testProcess() throws Exception {
        String dir = "datasets/Harmonizome/download/gene_similarity_matrix_cosine";
        File srcFile = new File(dir, "chea.txt");
        DataProcessor processor = new DataProcessor();
        process(srcFile, dir, processor);
    }
    
    /**
     * Process the downloaded data to get the filtered correlations into another file.
     * @param file
     * @throws Exception
     */
    private void process(File srcFile,
                        String resultDir,
                        DataProcessor processor) throws Exception {
        logger.info("Starting processing " + srcFile);
        String fileName = srcFile.getName().split("\\.")[0];
        File processed = new File(resultDir, fileName + "_processed.txt");
        File filtered = new File(resultDir, fileName + "_filtered.txt");
        if (processor.processCorrelations(srcFile, processed, filtered)) {
            logger.info(srcFile.getName() + " is processed successfully!");
            // We can just delete the original downloaded file
            srcFile.delete();
        }
        else {
            logger.info(srcFile.getName() + " is not processed successfully!");
        }
    }
    
    public List<String> loadReactomeIDGDatasets(boolean filterToReactome) throws IOException {
        InputStream is = getDataResource();
        Scanner scanner = new Scanner(is);
        List<String> rtn = new ArrayList<>();
        String line = scanner.nextLine(); // Escape the first head line
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();
            String[] tokens = line.split("\t");
            if (filterToReactome && tokens[13].equalsIgnoreCase("false"))
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
    
    /**
     * Default loading is implemented for loading data that can be used by Reactome only.
     * @return
     * @throws IOException
     */
    public List<String> loadReactomeIDGDatasets() throws IOException {
        return loadReactomeIDGDatasets(true);
    }

    private InputStream getDataResource() {
        String fileName = "Harmonizome_datasets_annotations_062819.txt";
        InputStream is = getClass().getClassLoader().getResourceAsStream(fileName);
        return is;
    }
    
    private Map<String, String> getDataSet2Path() throws Exception {
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
