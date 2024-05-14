package proteomics;

import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.SpecProcessor;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.Peptide;
import proteomics.Types.PeptideInfo;
import proteomics.Types.VarPtm;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "3";
    static final boolean useXcorr = false;
    public static final int minTagLenToReduceProtDb = 5;
    public static final boolean isDebugMode = java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0;


    static int MIN_TAG_LEN;
    public static double MIN_AA_INTENSITY;
    private int NUM_SIG_PTM;
    public static HashSet<Integer> lszDebugScanNum = new HashSet<>(Arrays.asList(4127));

    public static void main(String[] args) {
        long startTime = System.nanoTime();
        if (args.length != 1) {
            help();
        }
        String parameterPath = args[0].trim();
        logger.info("Running PIPI version {}.", versionStr);

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }
        try {
            dbName = String.format(Locale.US, "PIPI.temp.db");
            new PIPI(parameterPath, dbName, hostName);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String dbName, String hostName) throws Exception {
        ch.qos.logback.classic.Logger root = (ch.qos.logback.classic.Logger) org.slf4j.LoggerFactory.getLogger(ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
        root.setLevel(ch.qos.logback.classic.Level.DEBUG);
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tol = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tol = Double.valueOf(parameterMap.get("ms1_tolerance"));
        boolean addContaminants = Boolean.valueOf(parameterMap.get("add_contaminant"));
        double proteinCovThres = Double.valueOf(parameterMap.get("min_prot_coverage"));
        String spectraPath = parameterMap.get("spectra_path");
        String outputDir = parameterMap.get("output_dir");
        MIN_TAG_LEN = Integer.valueOf(parameterMap.get("min_tag_len"));
        MIN_AA_INTENSITY = (MIN_TAG_LEN == 3) ? 0.0 : 0.4;
        NUM_SIG_PTM = 18;

        logger.info("Parameter: {}.", parameterPath);

        Set<Integer> msLevelSet = new HashSet<>();
        msLevelSet.add(2);

        logger.info("Loading parameters and build fmIndex...");
        BuildIndex buildIndex = new BuildIndex(parameterMap);
        MassTool massTool = buildIndex.returnMassTool();
        InferPTM inferPTM = buildIndex.getInferPTM();

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        DatasetReader datasetReader;
        JMzReader[] spectraParserArray;
        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();
        Map<Integer, String> fileIdNameMap = new HashMap<>();
        Map<String, Integer> fileNameIdMap = new HashMap<>();
        if ((!spectraFile.exists())) {
            throw new FileNotFoundException("The spectra file not found.");
        }

        if (!spectraFile.isDirectory()) {
            spectraParserArray = new JMzReader[1];
            JMzReader spectraParser;
            String ext = spectraPath.substring(spectraPath.lastIndexOf(".") + 1);
            if (ext.contentEquals("mzXML")) {
                spectraParser = new MzXMLFile(spectraFile);
            } else if (ext.toLowerCase().contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            } else {
                throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
            }
            spectraParserArray[0] = spectraParser;
            fileIdNameMap.put(0, spectraPath.substring(spectraPath.lastIndexOf("/") + 1).split("\\.")[0].replaceAll("\\.", "_"));
            fileNameIdMap.put(spectraPath.substring(spectraPath.lastIndexOf("/") + 1).split("\\.")[0].replaceAll("\\.", "_"), 0);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tol, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        } else {
            String[] fileList = spectraFile.list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".mgf");
                }
            });
            spectraParserArray = new JMzReader[fileList.length];
            for (int i = 0; i < fileList.length; i++) {
                spectraParserArray[i] = new MgfFile(new File(spectraPath + fileList[i]));
                fileIdNameMap.put(i, fileList[i].split("\\.")[0].replaceAll("\\.", "_"));
                fileNameIdMap.put(fileList[i].split("\\.")[0].replaceAll("\\.", "_"), i);
            }
            String ext = fileList[0].substring(fileList[0].lastIndexOf(".") + 1);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tol, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        }

        SpecProcessor specProcessor = new SpecProcessor(massTool);
        Map<String, Integer> precursorChargeMap = new HashMap<>();
        Map<String, Double> precursorMassMap = new HashMap<>();

        logger.info("Preprocessing protein database...");
        int p_thread_num = Runtime.getRuntime().availableProcessors();

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Integer> precursorScanNoMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        ExecutorService threadPoolGetLongTag = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<GetLongTag.Entry>> taskListGetLongTag = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoderX = DriverManager.getConnection(sqlPath);
        Statement sqlStateGetLongTag = sqlConSpecCoderX.createStatement();
        ResultSet sqlResSetGetLongTag = sqlStateGetLongTag.executeQuery("SELECT scanName, scanNum, precursorCharge" +
                ", precursorMass, precursorScanNo, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient FROM spectraTable");

        ReentrantLock lockGetLongTag = new ReentrantLock();
        int submitNum = 0;
        Set<String> validScanSet = new HashSet<>();
        while (sqlResSetGetLongTag.next()) {
            String scanName = sqlResSetGetLongTag.getString("scanName");
            int scanNum = sqlResSetGetLongTag.getInt("scanNum");
            int precursorCharge = sqlResSetGetLongTag.getInt("precursorCharge");
            double precursorMass = sqlResSetGetLongTag.getDouble("precursorMass");
            mgfTitleMap.put(scanName, sqlResSetGetLongTag.getString("mgfTitle"));
            isotopeCorrectionNumMap.put(scanName, sqlResSetGetLongTag.getInt("isotopeCorrectionNum"));
            precursorScanNoMap.put(scanName, sqlResSetGetLongTag.getInt("precursorScanNo"));
            ms1PearsonCorrelationCoefficientMap.put(scanName, sqlResSetGetLongTag.getDouble("ms1PearsonCorrelationCoefficient"));
            if (isDebugMode) {
                boolean shouldRun = false;
                for (int debugScanNum : lszDebugScanNum) {
                    if (Math.abs(scanNum - debugScanNum) < 2) {
                        shouldRun = true;
                    }
                }
                if (!shouldRun) continue;
            }
            int fileId = fileNameIdMap.get(scanName.split("\\.")[0]);
            precursorChargeMap.put(scanName, precursorCharge);
            precursorMassMap.put(scanName, precursorMass);
            submitNum++;
            validScanSet.add(scanName);
            taskListGetLongTag.add(threadPoolGetLongTag.submit(new GetLongTag(scanNum, buildIndex, massTool, spectraParserArray[fileId], lockGetLongTag, scanName, precursorCharge
                    , precursorMass, specProcessor, ms2Tol)));
        }
        sqlResSetGetLongTag.close();
        sqlStateGetLongTag.close();

        int lastProgressGetLongTag = 0;
        int totalCountGetLongTag = taskListGetLongTag.size();
        int countGetLongTag = 0;
        Map<String, List<GetLongTag.TagRes>> prot_TagResList_Map = new HashMap<>();
        while (countGetLongTag < totalCountGetLongTag) {
            List<Future<GetLongTag.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountGetLongTag - countGetLongTag);
            for (Future<GetLongTag.Entry> task : taskListGetLongTag) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        GetLongTag.Entry entry = task.get();
                        for (String prot : entry.prot_TagResList_Map.keySet()) {
                            List<GetLongTag.TagRes> tagResList = entry.prot_TagResList_Map.get(prot);
                            if (prot_TagResList_Map.containsKey(prot)) {
                                prot_TagResList_Map.get(prot).addAll(tagResList);
                            } else {
                                List<GetLongTag.TagRes> tmpTagResList = new LinkedList<>();
                                tmpTagResList.addAll(tagResList);
                                prot_TagResList_Map.put(prot, tagResList);
                            }
                        }
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countGetLongTag += toBeDeleteTaskList.size();
            taskListGetLongTag.removeAll(toBeDeleteTaskList);
            taskListGetLongTag.trimToSize();

            int progress = countGetLongTag * 20 / totalCountGetLongTag;
            if (progress != lastProgressGetLongTag) {
                logger.info("Getting long tags for prot {}%...", progress * 5);
                lastProgressGetLongTag = progress;
            }

            if (countGetLongTag == totalCountGetLongTag) {
                break;
            }
            Thread.sleep(6000);
        }
        threadPoolGetLongTag.shutdown();
        if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolGetLongTag.shutdownNow();
            if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        if (lockGetLongTag.isLocked()) {
            lockGetLongTag.unlock();
        }

        Map<String, Integer> protLengthMap = buildIndex.protLengthMap;
        Map<String, Double> protCoverageMap = new HashMap<>();
        for (String prot : prot_TagResList_Map.keySet()) {
            int protLen = buildIndex.protSeqMap.get(prot).length();
            List<GetLongTag.TagRes> tagResList = prot_TagResList_Map.get(prot);
            Map<Integer, Double> posProductMap = new HashMap<>(protLen);
            for (int pos = 0; pos < protLen; pos++) {
                posProductMap.put(pos, 1.0);
            }
            for (GetLongTag.TagRes tagRes : tagResList) {
                int relPos = tagRes.relPos;
                List<Double> normedIaaList = tagRes.normedIaaList;
                for (int aaPos = relPos; aaPos < relPos + tagRes.tagSeq.length(); aaPos++) {
                    posProductMap.put(aaPos, posProductMap.get(aaPos) * (1 - normedIaaList.get(aaPos - relPos)));
                }
            }
            double tempCoverage = 0;
            for (int pos : posProductMap.keySet()) {
                tempCoverage += 1 - posProductMap.get(pos);
            }
            protCoverageMap.put(prot, tempCoverage / protLen);
        }

        List<Pair<String, Double>> protScoreLongList = new ArrayList<>();
        for (String prot : protCoverageMap.keySet()) {
            protScoreLongList.add(new Pair<>(prot, protCoverageMap.get(prot)));
        }
        Collections.sort(protScoreLongList, Comparator.comparing(o -> o.getSecond(), Comparator.reverseOrder()));

        Set<String> reducedProtIdSet = new HashSet<>();

        for (Pair<String, Double> pair : protScoreLongList) {
            if (pair.getSecond() < proteinCovThres) break;
            reducedProtIdSet.add(pair.getFirst());
        }

        if (isDebugMode) reducedProtIdSet = protLengthMap.keySet();

        logger.info("Building decoy...");
        String dbPath = parameterMap.get("database");
        String decoyFullDBPath = dbPath + ".decoy.fasta";
        File decoyFullFile = new File(decoyFullDBPath);
        if (decoyFullFile.exists()) {
            DbTool decoyDbTool = new DbTool(decoyFullDBPath, "others");
            Map<String, String> decoyProtSeqMap = decoyDbTool.getProtSeqMap();
            for (String decoyProtId : decoyProtSeqMap.keySet()) {
                String targetProtId = decoyProtId.split("DECOY_")[1];
                if (reducedProtIdSet.contains(targetProtId)) {
                    buildIndex.protSeqMap.put(decoyProtId, decoyProtSeqMap.get(decoyProtId));
                }
            }
        } else {
            BufferedWriter decoyWriter = new BufferedWriter(new FileWriter(decoyFullDBPath));
            ExecutorService threadPoolBuildDecoyProts = Executors.newFixedThreadPool(p_thread_num);
            ArrayList<Future<BuildDecoyProts.Entry>> taskListBuildDecoyProts = new ArrayList<>(reducedProtIdSet.size() + 10);
            ReentrantLock lockSpecCoder = new ReentrantLock();
            for (String protId : protLengthMap.keySet()) {
                taskListBuildDecoyProts.add(threadPoolBuildDecoyProts.submit(new BuildDecoyProts(parameterMap, buildIndex, protId)));
            }
            int lastProgressBuildDecoyProts = 0;
            int totalCountBuildDecoyProts = taskListBuildDecoyProts.size();
            int countBuildDecoyProts = 0;
            while (countBuildDecoyProts < totalCountBuildDecoyProts) {
                List<Future<BuildDecoyProts.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBuildDecoyProts - countBuildDecoyProts);
                for (Future<BuildDecoyProts.Entry> task : taskListBuildDecoyProts) {
                    if (task.isDone()) {
                        if (task.get() != null) {
                            BuildDecoyProts.Entry entry = task.get();
                            String protId = entry.protId;
                            String decoyProtSeq = entry.decoyProtSeq;
                            String decoyProtId = "DECOY_" + protId;
                            if (reducedProtIdSet.contains(protId)) buildIndex.protSeqMap.put(decoyProtId, decoyProtSeq);
                            decoyWriter.write(String.format(Locale.US, ">%s\n", decoyProtId));
                            decoyWriter.write(decoyProtSeq + "\n");
                        }
                        toBeDeleteTaskList.add(task);
                    } else if (task.isCancelled()) {
                        toBeDeleteTaskList.add(task);
                    }
                }
                countBuildDecoyProts += toBeDeleteTaskList.size();
                taskListBuildDecoyProts.removeAll(toBeDeleteTaskList);
                taskListBuildDecoyProts.trimToSize();

                int progress = countBuildDecoyProts * 20 / totalCountBuildDecoyProts;
                if (progress != lastProgressBuildDecoyProts) {
                    logger.info("Build Decoy Prots {}%...", progress * 5);
                    lastProgressBuildDecoyProts = progress;
                }

                if (countBuildDecoyProts == totalCountBuildDecoyProts) {
                    break;
                }
                Thread.sleep(6000);
            }
            threadPoolBuildDecoyProts.shutdown();
            if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS)) {
                threadPoolBuildDecoyProts.shutdownNow();
                if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS))
                    throw new Exception("Pool did not terminate");
            }
            if (lockSpecCoder.isLocked()) {
                lockSpecCoder.unlock();
            }
            decoyWriter.close();
        }

        Iterator<String> iter = buildIndex.protSeqMap.keySet().iterator();
        while (iter.hasNext()) {
            String pureId = iter.next().replace("DECOY_", "");
            if (!reducedProtIdSet.contains(pureId)) {
                iter.remove();
            }
        }
        buildBiDirectionFMIndex(buildIndex);
        buildIndex.minPeptideMass = 0;
        buildIndex.maxPeptideMass = 9999;
        Map<String, String> proteinAnnotationMap;

        DbTool dbTool = buildIndex.dbTool;
        DbTool contaminantsDb;
        if (addContaminants) {
            contaminantsDb = new DbTool(null, "contaminants");
            proteinAnnotationMap = contaminantsDb.getProteinAnnotateMap();
            proteinAnnotationMap.putAll(dbTool.getProteinAnnotateMap());
        } else {
            proteinAnnotationMap = dbTool.getProteinAnnotateMap();
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
        for (String protId : buildIndex.protSeqMap.keySet()) {
            writer.write(String.format(Locale.US, ">%s %s\n", protId, proteinAnnotationMap.getOrDefault(protId, "")));
            writer.write(buildIndex.protSeqMap.get(protId) + "\n");
        }
        writer.close();

        logger.info("Main searching...");

        p_thread_num = Integer.valueOf(parameterMap.get("PIPI_thread_num"));
        if (Integer.valueOf(parameterMap.get("PIPI_thread_num")) == 0) {
            p_thread_num = (int) (Runtime.getRuntime().availableProcessors() / 2);
        }
        if (isDebugMode) p_thread_num = 1;
        logger.debug("availabel processors for Main search" + Runtime.getRuntime().availableProcessors());
        logger.debug("thread num for Main search" + p_thread_num);
        ExecutorService threadPoolBone = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<MainSearch.Entry>> taskListBone = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockBone = new ReentrantLock();
        int submitNumBone = 0;
        for (String scanName : validScanSet) {
            String[] scanNameStr = scanName.split("\\.");
            int precursorCharge = precursorChargeMap.get(scanName);

            double precursorMass = precursorMassMap.get(scanName);
            int scanNum = Integer.valueOf(scanNameStr[2]);
            boolean shouldRun = false;
            if (isDebugMode) {
                for (int debugScanNum : lszDebugScanNum) {
                    if (Math.abs(scanNum - debugScanNum) < 1) {
                        shouldRun = true;
                    }
                }
                if (!shouldRun) continue;
            }
            submitNumBone++;
            int fileId = fileNameIdMap.get(scanNameStr[0]);
            taskListBone.add(threadPoolBone.submit(new MainSearch(scanNum, buildIndex, massTool, ms2Tol, ms1Tol, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass()
                    , spectraParserArray[fileId], lockBone, scanName, precursorCharge, precursorMass, specProcessor)));
        }
        logger.info(submitNumBone + " MS2 submitted to Main Search.");

        Map<String, PeptideInfo> allPeptideInfoMap = new HashMap<>();
        int lastProgressBone = 0;
        int resultCountBone = 0;
        int totalCountBone = taskListBone.size();
        int countBone = 0;

        Map<Integer, TreeMap<Double, Set<String>>> fileId_pcMassScanNameMap = new HashMap<>();
        for (int fileId : fileIdNameMap.keySet()) {
            fileId_pcMassScanNameMap.put(fileId, new TreeMap<>());
        }
        Map<String, Peptide> scanName_TopPeptide_Map = new HashMap<>();
        Map<String, String> scanName_PepString_Map = new HashMap<>();
        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatement = sqlConSpecCoder.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient,  precursorScanNo, labelling, peptide, theoMass, score,peptideSet) VALUES (?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSpecCoder.setAutoCommit(false);
        Map<VarPtm, Integer> varPtmCountMap = new HashMap<>();
        Map<String, Integer> ptmWithPos_countMap = new HashMap<>();
        while (countBone < totalCountBone) {
            List<Future<MainSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBone - countBone);
            for (Future<MainSearch.Entry> task : taskListBone) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        MainSearch.Entry entry = task.get();
                        for (VarPtm varPtm : entry.varPtmList) {
                            if (varPtmCountMap.containsKey(varPtm)) {
                                varPtmCountMap.put(varPtm, varPtmCountMap.get(varPtm) + 1);
                            } else {
                                varPtmCountMap.put(varPtm, 1);
                            }
                            String ptmWithPosStr;
                            if (varPtm.position == 4) {
                                ptmWithPosStr = varPtm.site + "(" + InferPTM.df3.format(varPtm.mass) + ")" + varPtm.position;
                            } else if (varPtm.position == 3) {
                                ptmWithPosStr = "pepC" + "(" + InferPTM.df3.format(varPtm.mass) + ")" + varPtm.position;
                            } else if (varPtm.position == 2) {
                                ptmWithPosStr = "pepN" + "(" + InferPTM.df3.format(varPtm.mass) + ")" + varPtm.position;
                            } else if (varPtm.position == 1) {
                                ptmWithPosStr = "protC" + "(" + InferPTM.df3.format(varPtm.mass) + ")" + varPtm.position;
                            } else {
                                ptmWithPosStr = "protN" + "(" + InferPTM.df3.format(varPtm.mass) + ")" + varPtm.position;
                            }
                            if (ptmWithPos_countMap.containsKey(ptmWithPosStr)) {
                                ptmWithPos_countMap.put(ptmWithPosStr, ptmWithPos_countMap.get(ptmWithPosStr) + 1);
                            } else {
                                ptmWithPos_countMap.put(ptmWithPosStr, 1);
                            }
                        }
                        allPeptideInfoMap.putAll(entry.peptideInfoMapForRef);
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanName);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
                        sqlPreparedStatement.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatement.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(8, precursorScanNoMap.get(entry.scanName));
                        sqlPreparedStatement.setString(9, entry.labelling);
                        sqlPreparedStatement.setString(10, entry.peptide);
                        sqlPreparedStatement.setDouble(11, entry.theoMass);
                        sqlPreparedStatement.setDouble(12, entry.score);
                        sqlPreparedStatement.setString(13, entry.peptideSet);

                        sqlPreparedStatement.executeUpdate();
                        int fileId = fileNameIdMap.get(entry.scanName.split("\\.")[0]);
                        TreeMap<Double, Set<String>> local_pcMassScanNameMap = fileId_pcMassScanNameMap.get(fileId);
                        if (local_pcMassScanNameMap.containsKey(entry.precursorMass)) {
                            local_pcMassScanNameMap.get(entry.precursorMass).add(entry.scanName);
                        } else {
                            Set<String> scanNumSet = new HashSet<>();
                            scanNumSet.add(entry.scanName);
                            local_pcMassScanNameMap.put(entry.precursorMass, scanNumSet);
                        }
                        scanName_TopPeptide_Map.put(entry.scanName, entry.topPeptide);
                        scanName_PepString_Map.put(entry.scanName, entry.peptideSet);

                        ++resultCountBone;
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countBone += toBeDeleteTaskList.size();
            taskListBone.removeAll(toBeDeleteTaskList);
            taskListBone.trimToSize();

            sqlConSpecCoder.commit();
            int progress = countBone * 20 / totalCountBone;
            if (progress != lastProgressBone) {
                logger.info("Main Search {}%...", progress * 5);
                lastProgressBone = progress;
            }

            if (countBone == totalCountBone) {
                break;
            }
            Thread.sleep(6000);
        }

        threadPoolBone.shutdown();
        if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolBone.shutdownNow();
            if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        sqlConSpecCoder.commit();
        sqlConSpecCoder.setAutoCommit(true);
        sqlConSpecCoder.close();
        if (lockBone.isLocked()) {
            lockBone.unlock();
        }
        logger.info(resultCountBone + " MS2 finished in Main Search.");

        logger.info("Supplementary searching...");
        p_thread_num = 3 + Runtime.getRuntime().availableProcessors();
        ExecutorService threadPoolSupp = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<SuppSearch.Entry>> taskListSupp = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockSupp = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
        int submitTimeSupp = 0;

        for (String thisScanName : validScanSet) {
            int thisFileId = fileNameIdMap.get(thisScanName.split("\\.")[0]);
            double mass = precursorMassMap.get(thisScanName);
            int thisScanNum = Integer.valueOf(thisScanName.split("\\.")[2]);

            String thisPtmSeq = "XFAKEPEP";
            if (scanName_TopPeptide_Map.containsKey(thisScanName)) {
                thisPtmSeq = scanName_TopPeptide_Map.get(thisScanName).getVarPtmContainingSeqNow();
            }
            Map<String, TreeMap<Integer, VarPtm>> ptmSeq_posVarPtmMap_Map = new HashMap<>();
            Map<String, PeptideInfo> local_PepSeq_PeptideInfo_Map = new HashMap<>();

            TreeMap<Double, Set<String>> local_other_pcMassScanNameMap = fileId_pcMassScanNameMap.get(thisFileId);
            for (int i = 0; i <= 0; i++) {
                for (Set<String> otherScanNameSet : local_other_pcMassScanNameMap.subMap(mass + i * MassTool.PROTON - mass * ms1Tol * 1e-6, true, mass + i * MassTool.PROTON + mass * ms1Tol * 1e-6, true).values()) {
                    for (String otherScanName : otherScanNameSet) {
                        int otherScanNum = Integer.valueOf(otherScanName.split("\\.")[2]);
                        if (otherScanNum < thisScanNum + 2000 && otherScanNum > thisScanNum - 2000 && otherScanNum != thisScanNum) {
                            String otherPtmSeq = scanName_TopPeptide_Map.get(otherScanName).getVarPtmContainingSeqNow();
                            String otherFreeSeq = scanName_TopPeptide_Map.get(otherScanName).getFreeSeq();
                            if (otherPtmSeq.contentEquals(thisPtmSeq)) continue;

                            if (!ptmSeq_posVarPtmMap_Map.containsKey(otherPtmSeq)) {
                                ptmSeq_posVarPtmMap_Map.put(otherPtmSeq, scanName_TopPeptide_Map.get(otherScanName).posVarPtmResMap);
                                local_PepSeq_PeptideInfo_Map.put(otherFreeSeq, allPeptideInfoMap.get(otherFreeSeq).clone());
                            }
                        }
                    }
                }
            }

            if (ptmSeq_posVarPtmMap_Map.isEmpty()) {
                continue;
            }
            int precursorCharge = precursorChargeMap.get(thisScanName);
            double precursorMass = precursorMassMap.get(thisScanName);
            submitTimeSupp++;
            taskListSupp.add(threadPoolSupp.submit(new SuppSearch(thisScanNum, buildIndex, massTool
                    , ms2Tol, spectraParserArray[thisFileId], lockSupp, thisScanName, precursorCharge
                    , precursorMass, specProcessor, binomial, ptmSeq_posVarPtmMap_Map, local_PepSeq_PeptideInfo_Map)));
        }
        logger.info(submitTimeSupp + " MS2 submitted to Supplementary Search.");
        int lastProgressSupp = 0;
        int resultCountSupp = 0;
        int totalCountSupp = taskListSupp.size();
        int countSupp = 0;
        Connection sqlConSupp = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatementPTM = sqlConSupp.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient,  precursorScanNo, labelling, peptide, theoMass, score, peptideSet) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSupp.setAutoCommit(false);
        while (countSupp < totalCountSupp) {
            List<Future<SuppSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountSupp - countSupp);
            for (Future<SuppSearch.Entry> task : taskListSupp) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        SuppSearch.Entry entry = task.get();
                        if (scanName_TopPeptide_Map.containsKey(entry.scanName)
                                && entry.score < scanName_TopPeptide_Map.get(entry.scanName).getScore()) {
                            toBeDeleteTaskList.add(task);
                            continue;
                        }

                        sqlPreparedStatementPTM.setInt(1, entry.scanNum);
                        sqlPreparedStatementPTM.setString(2, entry.scanName);
                        sqlPreparedStatementPTM.setInt(3, entry.precursorCharge);
                        sqlPreparedStatementPTM.setDouble(4, entry.precursorMass);
                        sqlPreparedStatementPTM.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setInt(8, precursorScanNoMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setString(9, entry.labelling);
                        sqlPreparedStatementPTM.setString(10, entry.peptide);
                        sqlPreparedStatementPTM.setDouble(11, entry.theoMass);
                        sqlPreparedStatementPTM.setDouble(12, entry.score);
                        String peptideSetStr = null;
                        if (scanName_PepString_Map.containsKey(entry.scanName)) {
                            peptideSetStr = entry.peptideSet + "," + scanName_PepString_Map.get(entry.scanName);
                        } else {
                            peptideSetStr = entry.peptideSet;
                        }
                        sqlPreparedStatementPTM.setString(13, peptideSetStr);
                        sqlPreparedStatementPTM.executeUpdate();
                        ++resultCountSupp;
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countSupp += toBeDeleteTaskList.size();
            taskListSupp.removeAll(toBeDeleteTaskList);
            taskListSupp.trimToSize();

            sqlConSupp.commit();

            int progress = countSupp * 20 / totalCountSupp;
            if (progress != lastProgressSupp) {
                logger.info("Supplementary Searching {}%...", progress * 5);
                lastProgressSupp = progress;
            }

            if (countSupp == totalCountSupp) {
                break;
            }
            Thread.sleep(6000);
        }
        threadPoolSupp.shutdown();
        if (!threadPoolSupp.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolSupp.shutdownNow();
            if (!threadPoolSupp.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        sqlConSupp.commit();
        sqlConSupp.setAutoCommit(true);
        sqlConSupp.close();
        if (lockSupp.isLocked()) {
            lockSupp.unlock();
        }

        postProcessing(outputDir, allPeptideInfoMap, sqlPath, buildIndex.protSeqMap, massTool, hostName, varPtmCountMap, ptmWithPos_countMap);
        logger.info("Saving results...");
    }

    private void buildBiDirectionFMIndex(BuildIndex buildIndex) throws IOException {
        BufferedWriter writerProtNormal = new BufferedWriter(new FileWriter("catProtNormal.txt"));
        BufferedWriter writerProtReverse = new BufferedWriter(new FileWriter("catProtReverse.txt"));
        StringBuilder normalSb = new StringBuilder();
        int dotPos = 0;
        int dotNum = 0;
        buildIndex.dotPosArrNormal = new int[buildIndex.protSeqMap.keySet().size()];
        for (String protId : buildIndex.protSeqMap.keySet()) {
            buildIndex.dotPosArrNormal[dotNum] = dotPos;
            buildIndex.posProtMapNormal.put(dotNum, protId);
            String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
            buildIndex.protSeqMap.put(protId, protSeq);
            normalSb.append(".");
            normalSb.append(protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length() + 1;
        }
        writerProtNormal.write(normalSb.toString());
        writerProtReverse.write(normalSb.reverse().toString());
        writerProtNormal.close();
        writerProtReverse.close();
        char[] textNormal = buildIndex.loadFile("catProtNormal.txt", true);
        char[] textReverse = buildIndex.loadFile("catProtReverse.txt", true);
        buildIndex.fmIndexNormal = new FMIndex(textNormal);
        buildIndex.fmIndexReverse = new FMIndex(textReverse);
        buildIndex.textNormalLen = textNormal.length - 1;
    }

    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "Identification of Peptide with Multiple PTMs using Combinatorial Optimization.\r\n"
                + "Author: Shengzhi Lai\r\n"
                + "Email: slaiad@connect.ust.hk\r\n"
                + "java -Xmx8g -jar PIPI.jar parameter.def spectra_file output_directory\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }

    class ScanRes {
        public int scanNum;
        public double qValue = -0.1;
        public List<CandiScore> peptideInfoScoreList;
        public double expMass;
        public int charge;
        public String scanName;

        ScanRes(String scanName, int scanNum, double expMass, List<CandiScore> peptideInfoScoreList, int charge) {
            this.scanName = scanName;
            this.scanNum = scanNum;
            this.expMass = expMass;
            this.peptideInfoScoreList = peptideInfoScoreList;
            this.charge = charge;
        }
    }

    class CandiScore implements Comparable<CandiScore> {
        public String ptmContainingSeq;
        public PeptideInfo peptideInfo;
        public double pepScore;
        public double protScore = 0;
        public double varPtmTotalScore = 0;
        public String ptmPosStr;

        CandiScore(PeptideInfo peptideInfo, double pepScore, String ptmContainingSeq, String ptmPosStr) {
            this.peptideInfo = peptideInfo;
            this.pepScore = pepScore;
            this.ptmContainingSeq = ptmContainingSeq;
            this.ptmPosStr = ptmPosStr;
        }

        public void setVarPtmTotalScore(Map<String, Double> ptmWithPos_ScoreMap) {
            int startI = -1;
            double varPtmTotalScore = 0;
            Set<String> countedPtmStr = new HashSet<>();
            int n_PtmOnPep = 0;
            String[] ptmPositions = ptmPosStr.split("_");
            for (int anyI = 0; anyI < ptmContainingSeq.length(); anyI++) {
                char thisChar = ptmContainingSeq.charAt(anyI);
                if (thisChar == '(') {
                    n_PtmOnPep++;
                    startI = anyI - 1;
                } else if (thisChar == ')') {
                    String thisPtmWithPosStr;
                    String thisPtmPos = ptmPositions[n_PtmOnPep - 1];
                    if (thisPtmPos.contentEquals("4")) {
                        thisPtmWithPosStr = ptmContainingSeq.substring(startI, anyI + 1) + thisPtmPos;
                    } else if (thisPtmPos.contentEquals("3")) {
                        thisPtmWithPosStr = "pepC" + ptmContainingSeq.substring(startI + 1, anyI + 1) + thisPtmPos;
                    } else if (thisPtmPos.contentEquals("2")) {
                        thisPtmWithPosStr = "pepN" + ptmContainingSeq.substring(startI + 1, anyI + 1) + thisPtmPos;
                    } else if (thisPtmPos.contentEquals("1")) {
                        thisPtmWithPosStr = "protC" + ptmContainingSeq.substring(startI + 1, anyI + 1) + thisPtmPos;
                    } else {
                        thisPtmWithPosStr = "protN" + ptmContainingSeq.substring(startI + 1, anyI + 1) + thisPtmPos;
                    }
                    if (ptmWithPos_ScoreMap.containsKey(thisPtmWithPosStr)) {
                        varPtmTotalScore += ptmWithPos_ScoreMap.get(thisPtmWithPosStr);
                    }
                }
            }
            this.varPtmTotalScore = varPtmTotalScore / n_PtmOnPep;
        }

        public int compareTo(CandiScore o2) {
            double coe = 1;
            if (this.protScore < o2.protScore * coe) {
                return -1;
            } else if (this.protScore * coe > o2.protScore) {
                return 1;
            } else {
                if (this.varPtmTotalScore < o2.varPtmTotalScore) {
                    return -1;
                } else if (this.varPtmTotalScore > o2.varPtmTotalScore) {
                    return 1;
                } else {
                    if (this.pepScore < o2.pepScore) {
                        return -1;
                    } else if (this.pepScore > o2.pepScore) {
                        return 1;
                    }
                }
            }
            return 0;
        }
    }

    private void postProcessing(String outputDir, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath, Map<String, String> protSeqMap, MassTool massTool, String hostName, Map<VarPtm, Integer> varPtmCountMap, Map<String, Integer> ptmWithPos_countMap) throws IOException, SQLException {
        Map<String, Double> varPtmRefScoreMap = new HashMap<>();
        Set<String> byUserPtmStr = new HashSet<>();
        for (VarPtm varPtm : varPtmCountMap.keySet()) {
            if (varPtmCountMap.get(varPtm) == 1) {
                continue;
            }
            if (varPtm.classification.contentEquals("ByUser")) {
                byUserPtmStr.add(varPtm.site + "(" + InferPTM.df3.format(varPtm.mass) + ")");
            }
        }
        for (String ptmWithPos : ptmWithPos_countMap.keySet()) {
            if (ptmWithPos_countMap.get(ptmWithPos) == 1) {
                continue;
            }
            varPtmRefScoreMap.put(ptmWithPos, Math.sqrt(ptmWithPos_countMap.get(ptmWithPos)));
        }
        List<Map.Entry<String, Double>> testList = new ArrayList<>(varPtmRefScoreMap.entrySet());
        Collections.sort(testList, Map.Entry.comparingByValue(Comparator.reverseOrder()));
        int j = 0;
        for (Map.Entry<String, Double> entry : testList) {
            String varPtmStr = entry.getKey();
            if ((j > NUM_SIG_PTM || entry.getValue() < 10.0) && !byUserPtmStr.contains(varPtmStr)) {
                varPtmRefScoreMap.remove(varPtmStr);
            }
            j++;
        }

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, scanName, precursorCharge, precursorMass, peptide, theoMass, score, peptideSet FROM spectraTable");
        Map<String, Map<String, Double>> protPepScoreMap = new HashMap<>();

        List<ScanRes> scanResList = new ArrayList<>();
        List<Double> topScoreList = new LinkedList<>();
        while (sqlResultSet.next()) {
            String topPeptide = sqlResultSet.getString("peptide");
            if (!sqlResultSet.wasNull()) {
                int charge = sqlResultSet.getInt("precursorCharge");
                double expMass = sqlResultSet.getDouble("precursorMass");
                double score = sqlResultSet.getDouble("score");
                String peptideSet = sqlResultSet.getString("peptideSet");
                int scanNum = sqlResultSet.getInt("scanNum");
                String scanName = sqlResultSet.getString("scanName");

                String[] candiSetStr = peptideSet.split(",");
                int numPep = candiSetStr.length / 4;

                PeptideInfo pepInfo = allPeptideInfoMap.get(topPeptide.replaceAll("[^A-Z]+", ""));
                List<CandiScore> candiScoreList = new ArrayList<>();
                for (int i = 0; i < numPep; i++) {
                    String ptmContainingSeq = candiSetStr[4 * i + 0];
                    PeptideInfo candiPeptideInfo = allPeptideInfoMap.get(ptmContainingSeq.replaceAll("[^A-Z]+", ""));
                    double thisScore = Double.valueOf(candiSetStr[4 * i + 1]);
                    String ptmPosStr = candiSetStr[4 * i + 3];
                    CandiScore candiScore = new CandiScore(candiPeptideInfo, thisScore, ptmContainingSeq, ptmPosStr);
                    if (!ptmPosStr.contentEquals("emtry")) {
                        candiScore.setVarPtmTotalScore(varPtmRefScoreMap);
                    }
                    candiScoreList.add(candiScore);
                }
                Collections.sort(candiScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
                double topPepScore = candiScoreList.get(0).pepScore;
                topScoreList.add(topPepScore);
                Iterator<CandiScore> iter = candiScoreList.iterator();
                while (iter.hasNext()) {
                    CandiScore candiScore = iter.next();
                    if (candiScore.pepScore < 0.85 * topPepScore) {
                        iter.remove();
                    }
                }
                scanResList.add(new ScanRes(scanName, scanNum, expMass, candiScoreList, charge));

                for (String protId : pepInfo.protIdSet) {
                    if (protPepScoreMap.containsKey(protId)) {
                        Map<String, Double> pepScoreMap = protPepScoreMap.get(protId);
                        if (pepScoreMap.containsKey(topPeptide)) {
                            pepScoreMap.put(topPeptide, Math.max(pepScoreMap.get(topPeptide), score));
                        } else {
                            pepScoreMap.put(topPeptide, score);
                        }
                    } else {
                        Map<String, Double> pepScoreMap = new HashMap<>();
                        pepScoreMap.put(topPeptide, score);
                        protPepScoreMap.put(protId, pepScoreMap);
                    }
                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Collections.sort(topScoreList);
        int numScores = topScoreList.size();
        double medianScore = numScores % 2 == 0 ? (topScoreList.get(numScores / 2) + topScoreList.get(numScores / 2 - 1)) / 2 : topScoreList.get((numScores - 1) / 2);
        medianScore *= 1.5;
        Map<String, Double> protScoreMap = new HashMap<>();
        for (String protId : protPepScoreMap.keySet()) {
            Map<String, Double> pepScoreMap = protPepScoreMap.get(protId);
            for (String pep : pepScoreMap.keySet()) {
                if (pepScoreMap.get(pep) > medianScore) {
                    if (protScoreMap.containsKey(protId)) {
                        protScoreMap.put(protId, protScoreMap.get(protId) + pepScoreMap.get(pep));
                    } else {
                        protScoreMap.put(protId, pepScoreMap.get(pep));
                    }
                }
            }
        }

        Map<String, String> shortProtSeqMap = new HashMap<>();
        for (String protId : protSeqMap.keySet()) {
            String shortProtId = protId.split(" ")[0];
            shortProtSeqMap.put(shortProtId, protSeqMap.get(protId));
        }
        for (String protId : protScoreMap.keySet()) {
            protScoreMap.put(protId, protScoreMap.get(protId) / Math.log(shortProtSeqMap.get(protId).length()));
        }

        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            for (CandiScore candiScore : candiScoreList) {
                double protScoreForCand = -1;
                for (String protId : candiScore.peptideInfo.protIdSet) {
                    if (protScoreMap.getOrDefault(protId, 0.0) > protScoreForCand) {
                        protScoreForCand = protScoreMap.getOrDefault(protId, 0.0);
                    }
                }
                candiScore.protScore = protScoreForCand;
            }
        }

        DecimalFormat df = new DecimalFormat(".00000");

        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.reverseOrder());
        }
        Collections.sort(scanResList, Comparator.comparing(o -> o.peptideInfoScoreList.get(0).pepScore * Math.log10(o.peptideInfoScoreList.get(0).protScore + 1), Comparator.reverseOrder()));
        List<Double> fdrList = new ArrayList<>(scanResList.size());
        int numTPlusD = 0;
        int numD = 0;
        for (ScanRes scanRes : scanResList) {
            CandiScore candiScore = scanRes.peptideInfoScoreList.get(0);
            numTPlusD++;
            if (candiScore.peptideInfo.isDecoy) numD++;
            fdrList.add(2 * (double) numD / numTPlusD);
        }
        double minQ = 1.0;
        boolean found = false;
        for (int i = fdrList.size() - 1; i >= 0; i--) {
            minQ = Math.min(fdrList.get(i), minQ);
            scanResList.get(i).qValue = minQ;
            if (!found && minQ < 0.01) {
                found = true;
            }
        }

        if (true) {
            List<ScanRes> copy_scanResList = new ArrayList<>(scanResList);
            Set<String> topPeps = new HashSet<>();
            List<Integer> indexesToDel = new ArrayList<>();

            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                if (topPeps.contains(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq)) {
                    indexesToDel.add(i);
                } else {
                    topPeps.add(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq);
                }
            }
            indexesToDel.sort(Comparator.reverseOrder());
            for (int i : indexesToDel) {
                copy_scanResList.remove(i);
            }

            int newNumT = 0;
            int newNumD = 0;
            List<Double> newFdrList = new ArrayList<>();
            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                if (thisRes.peptideInfoScoreList.get(0).peptideInfo.isDecoy) {
                    newNumD++;
                } else {
                    newNumT++;
                }
                newFdrList.add(Math.min(1.0, (1.0 + newNumD) / newNumT));
            }
            for (int i = newFdrList.size() - 2; i >= 0; i--) {
                newFdrList.set(i, Math.min(newFdrList.get(i), newFdrList.get(i + 1)));
            }

            Map<String, Double> pepQMap = new HashMap<>();
            for (int i = 0; i < copy_scanResList.size(); i++) {
                ScanRes thisRes = copy_scanResList.get(i);
                pepQMap.put(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq, newFdrList.get(i));
            }
            for (int i = 0; i < scanResList.size(); i++) {
                ScanRes thisRes = scanResList.get(i);
                thisRes.qValue = pepQMap.get(thisRes.peptideInfoScoreList.get(0).peptideInfo.freeSeq);
            }
        }
        scanResList.sort(Comparator.comparing(o -> o.qValue, Comparator.naturalOrder()));

        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanResList.size());
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            CandiScore topCandi = candiScoreList.get(0);

            double theoMass = massTool.calResidueMass(topCandi.ptmContainingSeq) + massTool.H2O;
            double massDiff = getMassDiff(scanRes.expMass, theoMass, MassTool.C13_DIFF);
            double ppm = Math.abs(massDiff * 1e6 / theoMass);
            double finalScore = topCandi.pepScore * Math.log10(topCandi.protScore + 1);
            String finalStr = String.format(Locale.US, "%s,%d,%f,%d,%s,%s,%s,%s,%s,%s,%f,%f,%f,%d\n"
                    , scanRes.scanName, scanRes.scanNum, scanRes.qValue, topCandi.peptideInfo.isDecoy ? 0 : 1, df.format(finalScore), topCandi.ptmContainingSeq, topCandi.peptideInfo.freeSeq, df.format(topCandi.pepScore)
                    , String.join(";", topCandi.peptideInfo.protIdSet), df.format(topCandi.protScore), ppm, theoMass, scanRes.expMass, scanRes.charge
            );

            finalExcelList.add(new Pair(finalScore, finalStr));
        }

        Collections.sort(finalExcelList, Comparator.comparing(o -> o.getFirst(), Comparator.reverseOrder()));

        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(outputDir + "PIPI3." + hostName + ".csv");
        } catch (Exception e) {
            fileWriter = new FileWriter(outputDir + "PIPI3." + hostName + df.format(Math.random()) + ".csv");
        }
        BufferedWriter writer = new BufferedWriter(fileWriter);
        writer.write("scanName,scanNum,qValue,TorD,finalScore,peptide,freeSeq,pepScore,proteins,protScore,ppm,theoMass,expMass,charge\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
    }

    public static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }
}
