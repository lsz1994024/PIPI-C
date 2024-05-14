package proteomics;

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Score;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import com.google.common.collect.Sets;
import gurobi.GRB;
import gurobi.GRBEnv;
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMRes;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.MIN_TAG_LEN;
import static proteomics.PIPI.lszDebugScanNum;
import static proteomics.PTM.InferPTM.C_PART;
import static proteomics.PTM.InferPTM.N_PART;
import static proteomics.Segment.InferSegment.*;

public final class MainSearch implements Callable<MainSearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final JMzReader spectraParser;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor specProcessor;
    private final int scanNum;
    private final double ms2Tolerance;
    private final InferPTM inferPTM;

    public MainSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double minPtmMass, double maxPtmMass
            , JMzReader spectraParser, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.spectraParser = spectraParser;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.specProcessor = specProcessor;
        this.scanNum = scanNum;
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
        } finally {
            lock.unlock();
        }

        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass * ms1Tolerance / 1000000));

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, DatasetReader.topN);

        if (plMap.isEmpty()) return null;
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = specProcessor.prepareXCorr(finalPlMap, false);
        } else {
            expProcessedPL = specProcessor.digitizePL(finalPlMap);
        }
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, MIN_TAG_LEN);


        if (allLongTagList.isEmpty()) return null;

        double totalMass = precursorMass + 2 * MassTool.PROTON;


        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        int minTagLen = 3;

        GRBEnv env = new GRBEnv(true);
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") < 0) {
            env.set(GRB.IntParam.OutputFlag, 0);
            env.set(GRB.IntParam.LogToConsole, 0);
        }
        env.start();


        int tagNumsToTest = Math.min(20, allLongTagList.size());
        List<ExpTag> tagsToTest = new ArrayList<>();
        TreeSet<Peptide> resPeptTreeSet = new TreeSet<>(Collections.reverseOrder());
        ExpTag tagInfo;
        for (int i = 0; i < tagNumsToTest; i++) {
            tagInfo = allLongTagList.get(i);
            if (tagInfo.isNorC == N_TAG) {
                tagsToTest.add(tagInfo);
            } else if (tagInfo.isNorC == C_TAG) {
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
            } else {
                tagsToTest.add(tagInfo);
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
            }
        }
        List<OccGroup> occGroupList = new LinkedList<>();
        searchWithMultiTags(scanNum, tagsToTest, occGroupList, minTagLen, totalMass);
        occGroupList.sort(Comparator.comparing(o -> o.totalScore, Comparator.reverseOrder()));


        int count = 0;
        if (occGroupList.isEmpty()) {
            return null;
        }
        double topOccScore = occGroupList.get(0).totalScore;
        for (OccGroup occG : occGroupList) {
            if (occG.totalScore < 0.7 * topOccScore || count > 20) break;
            occG.tagPosList.sort(Comparator.comparing(o -> o.getSecond()));
            for (Pair<ExpTag, Integer> tagPosPair : occG.tagPosList) {
                ExpTag tmpTag = tagPosPair.getFirst();
                if (Math.abs(tmpTag.getHeadLocation() - MassTool.PROTON) < ms2Tolerance) {
                    tmpTag.isNorC = N_TAG;
                } else if (Math.abs(tmpTag.getTailLocation() - precursorMass + massTool.H2O - MassTool.PROTON) < ms1TolAbs + ms2Tolerance) {
                    tmpTag.isNorC = C_TAG;
                } else {
                    tmpTag.isNorC = NON_NC_TAG;
                }
            }
            String protId = occG.protId;
            count++;
            addCandisWithMultiMaxPeakMILP(scanNum, protId, occG.tagPosList, ms1TolAbs, resPeptTreeSet, peptideInfoMap, expProcessedPL, plMap, env);
        }
        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (resPeptTreeSet.isEmpty()) {
            env.dispose();
            return null;
        }
        List<Peptide> pepList = new ArrayList<>(resPeptTreeSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size() - 1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
                try {
                    if (isHomo(pepList.get(i), pepList.get(j), peptideInfoMap) || pepList.get(j).getScore() > 0.6 * pepList.get(i).getScore()) {
                        int iPriority = pepList.get(i).getPriority();
                        int jPriority = pepList.get(j).getPriority();
                        if (iPriority < jPriority) {
                            pepIdsToRemove.add(i);
                        } else if (iPriority > jPriority) {
                            pepIdsToRemove.add(j);
                        } else {
                            if (onlyDifferUnsettledPtm(pepList.get(i), pepList.get(j))) {
                                pepIdsToRemove.add(j);
                            }
                        }
                    }
                } catch (Exception e) {
                    System.out.println(scanNum + ", isHomo");
                }
            }
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++) {
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            env.dispose();
            return null;
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);
        Peptide topPep = peptideArray[0];

        Entry entry;
        String pepSetString = "";
        for (Peptide peptide : peptideArray) {
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            String ptmPosString;
            if (peptide.posVarPtmResMap.isEmpty()) {
                ptmPosString = "empty";
            } else {
                ptmPosString = "";
            }
            for (int pos : peptide.posVarPtmResMap.keySet()) {
                ptmPosString += peptide.posVarPtmResMap.get(pos).position + "_";
            }
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) + "," + ptmPosString.substring(0, ptmPosString.length() - 1) + ",";
        }

        double deltaLCn = 1;
        if (peptideArray.length > candisNum - 1) {
            deltaLCn = (peptideArray[0].getScore() - peptideArray[candisNum - 1].getScore()) / peptideArray[0].getScore();
        }
        double deltaCn = 1;
        if (peptideArray.length > 1) {
            for (int i = 0; i < peptideArray.length; i++) {
                if (peptideArray[i].getScore() != peptideArray[0].getScore()) {
                    deltaCn = (peptideArray[0].getScore() - peptideArray[i].getScore()) / peptideArray[0].getScore();
                    break;
                }
            }
        }
        String otherPtmPatterns = "-";
        entry = new MainSearch.Entry(
                scanNum, scanName, precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore()
                , pepSetString.substring(0, pepSetString.length() - 1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        entry.varPtmList.addAll(topPep.posVarPtmResMap.values());
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }
        int c = 1;
        env.dispose();
        return entry;
    }

    private boolean mergeTagPosList(List<Pair<ExpTag, Integer>> tagRelPosList, String protSeq) {
        boolean successAdded = true;
        tagRelPosList.sort(Comparator.comparing(o -> o.getSecond()));
        int lastTagId = 0;
        int contTagId = 1;
        List<Integer> tagIdAppended = new ArrayList<>();
        while (contTagId < tagRelPosList.size()) {
            int lastTagPos = tagRelPosList.get(lastTagId).getSecond();
            ExpTag lastTag = tagRelPosList.get(lastTagId).getFirst();
            int contTagPos = tagRelPosList.get(contTagId).getSecond();
            ExpTag contTag = tagRelPosList.get(contTagId).getFirst();

            if (contTagPos <= lastTagPos + lastTag.size()) {
                double lastFlagMass;
                double contFlagMass = contTag.getHeadLocation();
                if (lastTagPos == contTagPos) {
                    lastFlagMass = lastTag.getHeadLocation();
                } else {
                    lastFlagMass = lastTag.expAaList.get(contTagPos - lastTagPos - 1).getTailLocation();
                }
                if (Math.abs(contFlagMass - lastFlagMass) < ms2Tolerance) {
                    lastTag.appendTag(contTagPos - lastTagPos, contTag);
                    tagIdAppended.add(contTagId);
                    contTagId++;
                } else {
                    int badTadId = -1;
                    if (lastTag.isNorC == N_TAG) {
                        badTadId = contTagId;
                    } else if (contTag.isNorC == C_TAG) {
                        badTadId = lastTagId;
                    } else {
                        badTadId = lastTag.getTotalIntensity() < contTag.getTotalIntensity() ? lastTagId : contTagId;
                    }
                    tagIdAppended.add(badTadId);
                    if (badTadId == lastTagId) {
                        lastTagId = contTagId;
                        contTagId = lastTagId + 1;
                    } else {
                        contTagId++;
                    }
                    successAdded = false;
                }
            } else {
                int badTadId = -1;
                if (lastTag.isNorC == N_TAG) {
                    badTadId = contTagId;
                } else if (contTag.isNorC == C_TAG) {
                    badTadId = lastTagId;
                } else {
                    badTadId = lastTag.getTotalIntensity() < contTag.getTotalIntensity() ? lastTagId : contTagId;
                }

                if (isValidGap(lastTag, lastTagPos, contTag, contTagPos, protSeq)) {
                    lastTagId = contTagId;
                    contTagId = lastTagId + 1;
                } else {
                    tagIdAppended.add(badTadId);
                    if (badTadId == lastTagId) {
                        lastTagId = contTagId;
                        contTagId = lastTagId + 1;
                    } else {
                        contTagId++;
                    }
                    successAdded = false;
                }
            }
        }
        tagIdAppended.sort(Comparator.reverseOrder());
        for (int i : tagIdAppended) {
            tagRelPosList.remove(i);
        }
        Collections.sort(tagRelPosList, Comparator.comparing(o -> o.getSecond()));
        return successAdded;
    }

    private boolean isValidGap(ExpTag lastTag, int lastTagPos, ExpTag contTag, int contTagPos, String protSeq) {
        double deltaMass = contTag.getHeadLocation() - lastTag.getTailLocation();
        double sumAaMass = massTool.calResidueMass(protSeq.substring(lastTagPos + lastTag.size(), contTagPos));
        return deltaMass < sumAaMass + 250 && deltaMass > Math.max(57 * (contTagPos - lastTagPos - lastTag.size()), sumAaMass - 250);
    }

    final class OccGroup implements Cloneable {
        public List<Pair<ExpTag, Integer>> tagPosList = new ArrayList<>();
        public int lPos;
        public int rPos;
        public int spanLen;
        public String protId;
        public double totalScore;
        private int hashcode;

        public OccGroup(ExpTag initialTag, String protId, int tagPosInProt) {
            this.tagPosList.add(new Pair<>(initialTag, tagPosInProt));
            this.spanLen = initialTag.size();
            this.lPos = tagPosInProt;
            this.rPos = tagPosInProt + this.spanLen;
            this.protId = protId;
            this.totalScore = initialTag.getTotalIntensity();

            StringBuilder sb = new StringBuilder();
            sb.append(protId);
            sb.append(lPos);
            sb.append(rPos);
            sb.append(totalScore);
            for (int i = 0; i < this.tagPosList.size(); ++i) {
                sb.append(tagPosList.get(i).getFirst().getFreeAaString());
                sb.append(tagPosList.get(i).getFirst().getHeadLocation());
            }
            hashcode = sb.toString().hashCode();
        }

        public OccGroup() {
        }

        public boolean addTag(ExpTag newTag, int newTagPos) {
            boolean shouldAdd = true;
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList) {
                if (tagPosPair.getFirst().getFreeAaString().contentEquals(newTag.getFreeAaString())
                        && tagPosPair.getSecond() == newTagPos
                        && tagPosPair.getFirst().getHeadLocation() == newTag.getHeadLocation()) {
                    shouldAdd = false;
                    break;
                }
            }
            if (!shouldAdd) return true;
            this.tagPosList.add(new Pair<>(newTag, newTagPos));
            boolean successAdded = mergeTagPosList(this.tagPosList, buildIndex.protSeqMap.get(protId));

            this.totalScore = 0;
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList)
                this.totalScore += tagPosPair.getFirst().getTotalIntensity();
            this.lPos = this.tagPosList.get(0).getSecond();
            this.rPos = this.tagPosList.get(this.tagPosList.size() - 1).getSecond() + this.tagPosList.get(this.tagPosList.size() - 1).getFirst().size();
            this.spanLen = this.rPos - this.lPos;

            StringBuilder sb = new StringBuilder();
            sb.append(protId);
            sb.append(lPos);
            sb.append(rPos);
            sb.append(totalScore);
            for (int i = 0; i < this.tagPosList.size(); ++i) {
                sb.append(tagPosList.get(i).getFirst().getFreeAaString());
                sb.append(tagPosList.get(i).getFirst().getHeadLocation());
            }
            hashcode = sb.toString().hashCode();
            return successAdded;
        }

        public OccGroup clone() {
            OccGroup newOccG = new OccGroup();
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList) {
                newOccG.tagPosList.add(new Pair<>(tagPosPair.getFirst().clone(), tagPosPair.getSecond()));
            }
            newOccG.spanLen = this.spanLen;
            newOccG.lPos = this.lPos;
            newOccG.rPos = this.rPos;
            newOccG.protId = this.protId;
            newOccG.totalScore = this.totalScore;
            newOccG.hashcode = this.hashcode;
            return newOccG;
        }

        @Override
        public int hashCode() {
            return hashcode;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other)
                return true;
            if (other == null || getClass() != other.getClass())
                return false;
            OccGroup tmp = (OccGroup) other;
            return hashCode() == tmp.hashCode();
        }

    }

    private void searchWithMultiTags(int scanNum, List<ExpTag> tagsToTest, List<OccGroup> occGroupList, int minTagLen, double totalMass) {

        char[] tagChar;
        FMRes fmRes;
        int matchedPos, tagLen;
        String protId;
        ExpTag tagInfo;
        int tagsNum = tagsToTest.size();
        int occId = -1;
        for (int i = 0; i < tagsNum; i++) {
            boolean shouldTryForwardSearch = false;
            tagInfo = tagsToTest.get(i);
            tagLen = tagInfo.size();
            tagChar = tagInfo.getFreeAaString().toCharArray();
            fmRes = buildIndex.fmIndexNormal.fmSearchFuzzy(tagChar);
            if (fmRes != null) {
                matchedPos = fmRes.matchedPos;
                if (matchedPos != 0) {
                    shouldTryForwardSearch = true;
                }
                tagLen = tagInfo.size();
                if (matchedPos == 0 || tagLen - matchedPos >= minTagLen) {
                    ExpTag resTag = matchedPos == 0 ? tagInfo : tagInfo.subTag(matchedPos, tagLen);
                    int absTagPos, dotIndex, lPos, rPos;

                    Set<OccGroup> oldOccGToDel = new HashSet<>();
                    Set<OccGroup> newOccGToAdd = new HashSet<>();
                    for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                        absTagPos = buildIndex.fmIndexNormal.SA[ii];
                        dotIndex = -2 - Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                        protId = buildIndex.posProtMapNormal.get(dotIndex);
                        lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                        rPos = lPos + tagLen - matchedPos;
                        boolean belongToExistingOcc = false;

                        for (OccGroup occG : occGroupList) {
                            if (occG.protId.contentEquals(protId) && (lPos <= occG.rPos + occG.spanLen) && (rPos >= occG.lPos - occG.spanLen)) {
                                OccGroup newOccG = occG.clone();
                                boolean successAdded = newOccG.addTag(resTag.clone(), lPos);
                                if (successAdded) {
                                    oldOccGToDel.add(occG);
                                    newOccGToAdd.add(newOccG);
                                    belongToExistingOcc = true;
                                }
                            }
                        }
                        if (!belongToExistingOcc) {
                            OccGroup occG = new OccGroup(resTag.clone(), protId, lPos);
                            newOccGToAdd.add(occG);
                        }
                    }
                    for (OccGroup oldOccG : oldOccGToDel) occGroupList.remove(oldOccG);
                    for (OccGroup newOccG : newOccGToAdd) {
                        occGroupList.add(newOccG);
                    }
                }
            }

            if (shouldTryForwardSearch) {
                int absTagPos;
                int dotIndex;
                int lPos;
                int rPos;

                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagChar = revTagInfo.getFreeAaString().toCharArray();
                fmRes = buildIndex.fmIndexReverse.fmSearchFuzzy(tagChar);
                if (fmRes != null) {
                    matchedPos = fmRes.matchedPos;
                    if (matchedPos > 0 && tagInfo.size() - matchedPos < minTagLen) continue;

                    ExpTag resRevTag = matchedPos == 0 ? tagInfo : tagInfo.subTag(0, tagLen - matchedPos);
                    Set<OccGroup> oldOccGToDel = new HashSet<>();
                    Set<OccGroup> newOccGToAdd = new HashSet<>();
                    for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                        absTagPos = buildIndex.fmIndexReverse.SA[ii];
                        absTagPos = buildIndex.textNormalLen - absTagPos - tagLen + matchedPos;
                        dotIndex = -2 - Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                        protId = buildIndex.posProtMapNormal.get(dotIndex);
                        lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                        rPos = lPos + tagLen - matchedPos;

                        boolean foundRange = false;
                        for (OccGroup occG : occGroupList) {
                            if (occG.protId.contentEquals(protId) && (lPos <= occG.rPos + occG.spanLen) && (rPos >= occG.lPos - occG.spanLen)) {
                                OccGroup newOccG = occG.clone();
                                boolean successAdded = newOccG.addTag(resRevTag.clone(), lPos);
                                if (successAdded) {
                                    oldOccGToDel.add(occG);
                                    newOccGToAdd.add(newOccG);
                                    foundRange = true;
                                }
                            }
                        }
                        if (!foundRange) {
                            OccGroup occG = new OccGroup(resRevTag.clone(), protId, lPos);
                            newOccGToAdd.add(occG);
                        }
                    }
                    for (OccGroup oldOccG : oldOccGToDel) occGroupList.remove(oldOccG);
                    for (OccGroup newOccG : newOccGToAdd) {
                        occGroupList.add(newOccG);
                    }
                }
            }
        }

    }

    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {
        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3 * Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }

    private boolean onlyDifferUnsettledPtm(Peptide p1, Peptide p2) {
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        if (sbv1.dot(sbv2) < 0.3 * Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length())) {
            return false;
        }

        if (p1.posVarPtmResMap.size() != p2.posVarPtmResMap.size()) {
            return false;
        }
        if (!p1.posVarPtmResMap.keySet().containsAll(p2.posVarPtmResMap.keySet())) {
            return false;
        }
        byte n_SameMass = 0;
        for (int pos : p2.posVarPtmResMap.keySet()) {
            if (p1.posVarPtmResMap.get(pos).mass == p2.posVarPtmResMap.get(pos).mass) {
                n_SameMass++;
            } else if (Math.abs(p1.posVarPtmResMap.get(pos).mass - p2.posVarPtmResMap.get(pos).mass) < 0.02) {
                if (p1.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled") || p2.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled")) {
                    n_SameMass++;
                }
            }
        }
        return n_SameMass == p1.posVarPtmResMap.size();
    }

    private double addCandisWithMultiMaxPeakMILP(int scanNum, String protId, List<Pair<ExpTag, Integer>> tagRelPosList, double ms1TolAbs, TreeSet<Peptide> resPepTreeSet
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) {

        int tagNum;
        String protSeq = buildIndex.protSeqMap.get(protId);
        int protLen = protSeq.length();
        if (tagRelPosList.get(0).getSecond() < 0
                || (tagRelPosList.get(0).getSecond() == 0 && tagRelPosList.get(0).getFirst().isNorC != N_TAG)) {
            tagRelPosList.remove(0);
        }
        tagNum = tagRelPosList.size();
        if (tagNum == 0) return 0;
        if (tagRelPosList.get(tagNum - 1).getSecond() >= protLen
                || (tagRelPosList.get(tagNum - 1).getSecond() + tagRelPosList.get(tagNum - 1).getFirst().size() == protLen && tagRelPosList.get(tagNum - 1).getFirst().isNorC != C_TAG)) {
            tagRelPosList.remove(tagNum - 1);
        }
        tagNum = tagRelPosList.size();
        if (tagNum == 0) return 0;

        double maxScore = 0;

        Pair<ExpTag, Integer> lastTagPosPair = tagRelPosList.get(tagNum - 1);
        int remainMC = massTool.missedCleavage -
                getNumOfMCFromStr(protSeq.substring(tagRelPosList.get(0).getSecond(), lastTagPosPair.getSecond() + lastTagPosPair.getFirst().size()));
        if (lastTagPosPair.getFirst().isNorC == C_TAG && isKR(lastTagPosPair.getFirst().getFreeAaString().charAt(lastTagPosPair.getFirst().size() - 1))) {
            remainMC++;
        }

        List<TreeSet<Peptide>> segResList = new ArrayList<>();
        for (int gId = 0; gId < tagNum + 1; ++gId) {
            ExpTag lTag = null;
            int lTagPos = 0;
            if (gId != 0) {
                lTag = tagRelPosList.get(gId - 1).getFirst();
                lTagPos = tagRelPosList.get(gId - 1).getSecond();
            }
            ExpTag rTag = null;
            int rTagPos = 0;
            if (gId != tagNum) {
                rTag = tagRelPosList.get(gId).getFirst();
                rTagPos = tagRelPosList.get(gId).getSecond();
            }
            boolean solved;
            if (gId == 0) {
                if (Math.abs(rTag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs + ms2Tolerance) {
                    solved = true;
                } else {
                    TreeSet<Peptide> nModPepsSet = new TreeSet<>(Comparator.reverseOrder());
                    solved = solveGapN(scanNum, rTag, rTagPos, protSeq, remainMC, ms1TolAbs, expProcessedPL, plMap, env, nModPepsSet);
                    if (solved) segResList.add(nModPepsSet);
                }
                if (solved) segResList.add(getPepFromTag(rTag));
            } else if (gId == tagNum) {
                if (Math.abs(lTag.getTailLocation() - precursorMass + massTool.H2O - MassTool.PROTON) < ms1TolAbs + ms2Tolerance) {
                    solved = true;
                } else {
                    TreeSet<Peptide> cModPepsSet = new TreeSet<>(Comparator.reverseOrder());
                    solved = solveGapC(scanNum, lTag, lTagPos, protId, protSeq, remainMC,
                            ms1TolAbs, expProcessedPL, plMap, env, cModPepsSet);
                    if (solved) segResList.add(cModPepsSet);
                }
            } else {
                TreeSet<Peptide> midModPepsSet = new TreeSet<>(Comparator.reverseOrder());
                String gapSeq = protSeq.substring(tagRelPosList.get(gId - 1).getSecond() + tagRelPosList.get(gId - 1).getFirst().size(), tagRelPosList.get(gId).getSecond());
                double deltaMass = tagRelPosList.get(gId).getFirst().getHeadLocation() - tagRelPosList.get(gId - 1).getFirst().getTailLocation() - massTool.calResidueMass(gapSeq);
                solved = solveGapMid(scanNum, lTag, lTagPos, rTag, rTagPos, gapSeq, deltaMass, ms1TolAbs, expProcessedPL, plMap, env, midModPepsSet);
                if (solved) {
                    segResList.add(midModPepsSet);
                    segResList.add(getPepFromTag(rTag));
                }
            }
            if (!solved) {
                return 0;
            }
        }
        if (!segResList.isEmpty()) {
            collectResult(segResList, resPepTreeSet, expProcessedPL, plMap, peptideInfoMap, protId, protSeq, tagRelPosList.get(0).getSecond(), tagRelPosList.get(0).getFirst().isNorC == N_TAG, tagRelPosList.get(tagNum - 1).getSecond() + tagRelPosList.get(tagNum - 1).getFirst().size(), tagRelPosList.get(tagNum - 1).getFirst().isNorC == C_TAG);
        }
        return maxScore;
    }

    private void collectResult(List<TreeSet<Peptide>> segResList, TreeSet<Peptide> resPepTreeSet, SparseVector expProcessedPL, TreeMap<Double, Double> plMap,
                               Map<String, PeptideInfo> peptideInfoMap, String protId, String protSeq, int firstTagStartPos, boolean nIsTag, int lastTagEndPos, boolean cIsTag) {
        List<Set<Integer>> idSetList = new ArrayList<>(segResList.size());
        for (int segNum = 0; segNum < segResList.size(); segNum++) {
            Set<Integer> idSet = new HashSet<>();
            for (int id = 0; id < segResList.get(segNum).size(); ++id) {
                idSet.add(id);
            }
            idSetList.add(idSet);
        }

        for (List<Integer> resIdList : Sets.cartesianProduct(idSetList)) {
            StringBuilder pepSeqSB = new StringBuilder();
            PosMassMap posMassMap = new PosMassMap();
            TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();

            int curRelPos = 0;
            for (int i = 0; i < resIdList.size(); ++i) {
                int id = resIdList.get(i);
                Peptide[] segArr = new Peptide[segResList.get(i).size()];
                segResList.get(i).toArray(segArr);
                Peptide curPepSeg = segArr[id];
                pepSeqSB.append(curPepSeg.getFreeSeq());

                if (!curPepSeg.posVarPtmResMap.isEmpty()) {
                    for (int j : curPepSeg.posVarPtmResMap.keySet()) {
                        posMassMap.put(j + curRelPos, curPepSeg.posVarPtmResMap.get(j).mass);
                        posVarPtmResMap.put(j + curRelPos, curPepSeg.posVarPtmResMap.get(j));
                    }
                }
                curRelPos += segArr[id].getFreeSeq().length();
            }
            Peptide resPep = new Peptide(pepSeqSB.toString(), false, massTool);
            if (!posMassMap.isEmpty()) {
                resPep.setVarPTM(posMassMap);
                resPep.posVarPtmResMap.putAll(posVarPtmResMap);
            }
            double calScore = massTool.buildVectorAndCalXCorr(resPep.getIonMatrixNow(), 1, expProcessedPL, resPep.matchedBions, resPep.matchedYions);
            resPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, resPep.getIonMatrixNow(), ms2Tolerance));
            resPep.setScore(calScore * (1 - resPep.posVarPtmResMap.size() * 0.02));

            updatePeptideTreeSet(resPep, resPepTreeSet, peptideInfoMap, protId, protSeq, protSeq.indexOf(pepSeqSB.toString()), protSeq.indexOf(pepSeqSB.toString()) + pepSeqSB.length() - 1);
        }
    }

    private TreeSet<Peptide> getPepFromTag(ExpTag tag) {
        String freeSeq = tag.getFreeAaString();
        String modSeq = tag.getPtmAaString();
        Peptide peptide = new Peptide(freeSeq, false, massTool);
        TreeSet<Peptide> modPepPool = new TreeSet<>();
        modPepPool.add(peptide);
        if (!freeSeq.contentEquals(modSeq)) {
            PosMassMap posMassMap = new PosMassMap();
            int idOfAa = -1;
            for (char aaChar : tag.getPtmAaString().toCharArray()) {
                if (Character.isUpperCase(aaChar)) {
                    idOfAa += 1;
                } else {
                    posMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                    peptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                }
            }
        }
        return modPepPool;
    }

    private boolean solveGapC(int scanNum, ExpTag tag, int tagPosInProt, String protId, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> cModPepsSet) {

        double tagCMass = tag.getTailLocation() + massTool.H2O - MassTool.PROTON;
        if (tag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs + ms2Tolerance) return true;
        double cCutMass = tag.getTailLocation() - MassTool.PROTON;
        int tagLen = tag.size();
        int protLen = protSeq.length();
        List<Integer> cPoses = new ArrayList<>();
        List<Integer> krPoses = new ArrayList<>();
        int mcInC = 0;
        boolean startRecord = false;
        char aaChar;
        boolean isKR;

        for (int cPos = tagPosInProt + tagLen; cPos < protLen; cPos++) {
            aaChar = protSeq.charAt(cPos);
            if (badAA(aaChar)) break;
            isKR = isKR(aaChar);
            tagCMass += massTool.getMassTable().get(aaChar);
            if (mcInC > remainMC || tagCMass > precursorMass + maxPtmMass) break;
            if (tagCMass >= precursorMass - maxPtmMass) {
                if (startRecord) {
                    cPoses.add(cPos);
                } else {
                    startRecord = true;
                }
                if (isKR) krPoses.add(cPos);

                if (Math.abs(tagCMass - precursorMass) <= ms1TolAbs + ms2Tolerance) {
                    String cPartSeq = protSeq.substring(tagPosInProt + tagLen, cPos + 1);
                    storeCleanPartPeptides(cPartSeq, cCutMass, C_PART, unUsedExpProcessedPL, unUsedPlMap, cModPepsSet);
                    return true;
                }
            }
            if (isKR) mcInC++;
        }
        if (!startRecord) return false;
        Map<Integer, Integer> posYIdMap = new HashMap<>();
        Map<Integer, Integer> yIdMaxAbsPosMap = new HashMap<>();
        int optStartPos = 1;
        int optEndPosP1 = 0;
        if (buildIndex.cTermSpecific) {
            if (!krPoses.isEmpty()) {
                optEndPosP1 = krPoses.get(krPoses.size() - 1) + 1;
                optStartPos = krPoses.get(0) + 1;
                if (!cPoses.isEmpty() && cPoses.get(cPoses.size() - 1) == protLen - 1) {
                    optEndPosP1 = protLen;
                }
                int yId = 0;
                for (int pos = optStartPos; pos < optEndPosP1; pos++) {
                    posYIdMap.put(pos, yId);
                    int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
                    yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
                    if (krPoses.contains(pos)) {
                        yId++;
                    }
                }
            } else {
                if (cPoses.isEmpty()) {
                    return false;
                } else if (cPoses.get(cPoses.size() - 1) == protLen - 1) {
                    optEndPosP1 = protLen;
                    optStartPos = protLen;
                }
            }
        } else {
            optEndPosP1 = cPoses.get(cPoses.size() - 1) + 1;
            optStartPos = cPoses.get(0);
            int yId = 0;
            for (int pos = optStartPos; pos < optEndPosP1; pos++) {
                posYIdMap.put(pos, yId);
                int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
                yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
                yId++;
            }
        }

        if (optStartPos <= optEndPosP1) {
            int cPartStartPos = tagPosInProt + tagLen;
            String cPartSeq = protSeq.substring(cPartStartPos, optEndPosP1);
            int cPartSeqLen = cPartSeq.length();
            double flexiableMass = precursorMass - (tag.getTailLocation() + massTool.H2O - MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(cPartStartPos, optStartPos));
            Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
            Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(cPartSeqLen, 1);
            inferPTM.prepareInfoCTerm(scanNum, cPartSeq, absPos_MassVarPtm_Map, yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos, absPos_ptmPositions_Map);

            if (cPartSeqLen == 1) {
                findPtmOnOneAa(cModPepsSet, flexiableMass, absPos_MassVarPtm_Map, cPartSeq, cPartStartPos, ms1TolAbs);
            } else {
                inferPTM.findBestPtmMIPExtCPWL(scanNum, env, flexiableMass, cPartStartPos,
                        cPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map, optStartPos, optEndPosP1 == protLen
                        , yIdMaxAbsPosMap, unUsedPlMap, cModPepsSet, absPos_ptmPositions_Map);
            }
            if (yIdMaxAbsPosMap.isEmpty()) {
                inferPTM.findPossible1Ptm(scanNum, cPartSeq, absPos_MassVarPtm_Map, cPartStartPos, cModPepsSet, flexiableMass, ms1TolAbs);
            } else {
                inferPTM.findPossible1Ptm(scanNum, cPartSeq.substring(0, optStartPos - cPartStartPos), absPos_MassVarPtm_Map, cPartStartPos, cModPepsSet, flexiableMass, ms1TolAbs);
                for (int tmpYId : yIdMaxAbsPosMap.keySet()) {
                    int maxAbsPos = yIdMaxAbsPosMap.get(tmpYId);
                    double massToDeduct = massTool.calResidueMass(cPartSeq.substring(optStartPos - cPartStartPos, maxAbsPos + 1 - cPartStartPos));
                    inferPTM.findPossible1Ptm(scanNum, cPartSeq.substring(0, maxAbsPos + 1 - cPartStartPos), absPos_MassVarPtm_Map, cPartStartPos, cModPepsSet, flexiableMass - massToDeduct, ms1TolAbs);
                }
            }
        }

        return !cModPepsSet.isEmpty();
    }

    private boolean solveGapN(int scanNum, ExpTag tag, int tagPosInProt, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> nModPepsSet) {
        double nDeltaMass = tag.getHeadLocation() - MassTool.PROTON;
        double nCutMass = precursorMass + MassTool.PROTON - tag.getHeadLocation() - massTool.H2O;
        if (tag.isNorC == N_TAG || Math.abs(tag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs) return true;

        int mcInN = 0;
        char aaChar;

        int fixedStartPos = -1;
        int optStartPos = -1;
        for (int nPos = tagPosInProt - 1; nPos >= 0; nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (badAA(aaChar)) break;
            nDeltaMass -= massTool.getMassTable().get(aaChar);
            if (Math.abs(nDeltaMass) <= ms1TolAbs + ms2Tolerance) {
                String nPartSeq = protSeq.substring(nPos, tagPosInProt);
                storeCleanPartPeptides(nPartSeq, nCutMass, N_PART, unUsedExpProcessedPL, unUsedPlMap, nModPepsSet);
                return true;
            }
            if (nDeltaMass < maxPtmMass) {
                if (fixedStartPos == -1) {
                    fixedStartPos = nPos;
                    optStartPos = nPos;
                }
                if (nDeltaMass >= minPtmMass) {
                    optStartPos = nPos;
                } else {
                    break;
                }
            }
        }

        if (fixedStartPos == -1) return false;
        for (int nPos = tagPosInProt - 1; nPos >= fixedStartPos; nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isKR(aaChar)) {
                mcInN++;
                if (mcInN > remainMC) {
                    return false;
                }
            }
        }
        List<Integer> krPoses = new ArrayList<>();
        for (int nPos = fixedStartPos - 1; nPos >= Math.max(0, optStartPos - 1); nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isKR(aaChar)) {
                krPoses.add(nPos);
                mcInN++;
                if (mcInN > remainMC) {
                    optStartPos = nPos + 1;
                    break;
                }
            }
        }

        if (buildIndex.nTermSpecific) {
            if (fixedStartPos >= 2) {
                if (krPoses.isEmpty()) {
                    return false;
                } else {
                    fixedStartPos = 1 + krPoses.get(0);
                    optStartPos = 1 + krPoses.get(krPoses.size() - 1);
                }
            }
        }

        Map<Integer, Integer> posYIdMap = new HashMap<>();
        Map<Integer, Integer> yIdMinAbsPosMap = new HashMap<>();
        int yId = -1;

        if (fixedStartPos >= 2) {
            for (int pos = fixedStartPos - 1; pos >= optStartPos; pos--) {
                if (buildIndex.nTermSpecific) {
                    if (krPoses.contains(pos)) yId++;
                } else {
                    yId++;
                }
                int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                posYIdMap.put(pos, yId);
            }
        } else {
            for (int pos = fixedStartPos - 1; pos >= optStartPos; pos--) {
                yId++;
                int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                posYIdMap.put(pos, yId);
            }
        }

        if (optStartPos <= fixedStartPos) {
            String nPartSeq = protSeq.substring(optStartPos, tagPosInProt);
            int nPartSeqLen = tagPosInProt - optStartPos;
            double flexiableMass = tag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(fixedStartPos, tagPosInProt));
            Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);
            Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(nPartSeqLen, 1);

            inferPTM.prepareInfoNTerm(scanNum, nPartSeq,
                    absPos_MassVarPtm_Map, yIdMinAbsPosMap, optStartPos == 0, fixedStartPos, tagPosInProt, protSeq, absPos_ptmPositions_Map);

            if (nPartSeqLen == 1) {
                findPtmOnOneAa(nModPepsSet, flexiableMass, absPos_MassVarPtm_Map, nPartSeq, tagPosInProt - 1, ms1TolAbs);
            } else {
                inferPTM.findBestPtmMIPExtNPWL(scanNum, env, flexiableMass, tagPosInProt, nPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map
                        , fixedStartPos, optStartPos == 0, yIdMinAbsPosMap, protSeq, unUsedPlMap, nModPepsSet, absPos_ptmPositions_Map);
            }
            if (yIdMinAbsPosMap.isEmpty()) {
                inferPTM.findPossible1Ptm(scanNum, nPartSeq, absPos_MassVarPtm_Map, optStartPos, nModPepsSet, flexiableMass, ms1TolAbs);
            } else {
                inferPTM.findPossible1Ptm(scanNum, nPartSeq.substring(fixedStartPos - optStartPos), absPos_MassVarPtm_Map, fixedStartPos, nModPepsSet, flexiableMass, ms1TolAbs);
                for (int tmpYId : yIdMinAbsPosMap.keySet()) {
                    int minAbsPos = yIdMinAbsPosMap.get(tmpYId);
                    double massToDeduct = massTool.calResidueMass(nPartSeq.substring(minAbsPos - optStartPos, fixedStartPos - optStartPos));
                    inferPTM.findPossible1Ptm(scanNum, nPartSeq.substring(minAbsPos - optStartPos), absPos_MassVarPtm_Map, minAbsPos, nModPepsSet, flexiableMass - massToDeduct, ms1TolAbs);
                }
            }
        }
        return !nModPepsSet.isEmpty();
    }

    private boolean solveGapMid(int scanNum, ExpTag lTag, int lTagPos, ExpTag rTag, int rTagPos, String gapSeq, double deltaMass, double ms1TolAbs,
                                SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> midModPepsSet) {

        if (Math.abs(deltaMass) < 0.02) {
            Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
            tmpPeptide.setScore(0);
            midModPepsSet.add(tmpPeptide);
            return true;
        }
        int gapLen = gapSeq.length();
        Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(gapLen, 1);
        inferPTM.prepareInfoMid(scanNum, gapSeq, absPos_MassVarPtm_Map, lTagPos + lTag.size());
        if (gapLen == 1) {
            findPtmOnOneAa(midModPepsSet, deltaMass, absPos_MassVarPtm_Map, gapSeq, lTagPos + lTag.size(), ms1TolAbs);
        } else {
            double bIonRefMass = lTag.getTailLocation();
            double yIonRefMass = precursorMass + 2 * MassTool.PROTON - rTag.getHeadLocation();
            inferPTM.findBestPtmInGapPWL(scanNum, env, deltaMass, lTagPos + lTag.size(), rTagPos,
                    gapSeq, ms1TolAbs, absPos_MassVarPtm_Map, unUsedPlMap, midModPepsSet, bIonRefMass, yIonRefMass, precursorMass);
        }

        inferPTM.findPossible1Ptm(scanNum, gapSeq, absPos_MassVarPtm_Map, lTagPos + lTag.size(), midModPepsSet, deltaMass, ms1TolAbs);

        return !midModPepsSet.isEmpty();
    }

    private void findPtmOnOneAa(TreeSet<Peptide> midModPepsSet, double deltaMass, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int absPos, double ms1TolAbs) {
        TreeMap<Double, VarPtm> massPtmMap = absPos_MassVarPtm_Map.get(absPos);
        if (massPtmMap == null) return;
        for (double ptmMass : massPtmMap.keySet()) {
            if (Math.abs(ptmMass - deltaMass) < ms1TolAbs + ms2Tolerance) {
                Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(0, ptmMass);
                tmpPeptide.posVarPtmResMap.put(0, massPtmMap.get(ptmMass));
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(0);
                midModPepsSet.add(tmpPeptide);
            }
        }
    }

    private void updatePeptideTreeSet(Peptide newPeptide, TreeSet<Peptide> peptideTreeSet, Map<String, PeptideInfo> peptideInfoMap,
                                      String protId, String protSeq, int pepStartPos, int pepEndPos) {

        String shortProtId = protId.split(" ")[0];
        boolean added = false;
        if (peptideTreeSet.size() < candisNum) {
            peptideTreeSet.add(newPeptide);
            added = true;
        } else if (newPeptide.compareTo(peptideTreeSet.last()) == 1) {
            added = true;
            peptideTreeSet.pollLast();
            peptideTreeSet.add(newPeptide);
        }
        if (added) {
            char leftFlank = pepStartPos == 0 ? '-' : protSeq.charAt(pepStartPos - 1);
            char rightFlank = pepEndPos == protSeq.length() - 1 ? '-' : protSeq.charAt(pepEndPos + 1);
            String freePepSeq = newPeptide.getFreeSeq();
            PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
            if (peptideInfo != null) {
                if (!peptideInfo.protIdSet.contains(shortProtId)) {
                    peptideInfo.leftFlank = leftFlank;
                    peptideInfo.rightFlank = rightFlank;
                    peptideInfo.protIdSet.add(shortProtId);
                    if (!shortProtId.startsWith("DECOY_")) {
                        peptideInfo.isDecoy = false;
                    }
                }
            } else {
                peptideInfo = new PeptideInfo(freePepSeq, shortProtId.startsWith("DECOY_"), leftFlank, rightFlank);
                peptideInfo.protIdSet.add(shortProtId);
                peptideInfoMap.put(freePepSeq, peptideInfo);
            }
        }
    }

    private void storeCleanPartPeptides(String partSeq, double cutMass, byte isNCPart, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, TreeSet<Peptide> modPepsSet) {

        Peptide partPeptide = new Peptide(partSeq, false, massTool);
        double[][] ionMatrix = partPeptide.getIonMatrixNow();
        inferPTM.updateIonMatrix(ionMatrix, cutMass, isNCPart);
        Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length() - 1).boxed().collect(Collectors.toSet());
        double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, partPeptide.matchedBions, partPeptide.matchedYions, jRange);
        partPeptide.setScore(score * (1 - partPeptide.posVarPtmResMap.size() * 0.01));
        partPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ionMatrix, ms2Tolerance));
        if (modPepsSet.size() < 2) {
            modPepsSet.add(partPeptide);
        } else if (modPepsSet.last().compareTo(partPeptide) < 0) {
            modPepsSet.pollLast();
            modPepsSet.add(partPeptide);
        }
    }

    private int getNumOfMCFromStr(String tagStr) {
        String str2 = tagStr.replaceAll("[KR]", "");
        int numMC = tagStr.length() - str2.length();
        return numMC;
    }

    private boolean isKR(char aa) {
        return aa == 'K' || aa == 'R';
    }

    private boolean badAA(char aa) {
        return aa == 'X' || aa == 'B' || aa == 'J' || aa == 'Z' || aa == 'O' || aa == 'U';
    }

    public class Entry {
        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public Peptide topPeptide = null;
        final int scanNum;
        final String scanName;
        final int precursorCharge;
        final double precursorMass;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns;
        final String aScore;
        final String peptideSet;
        List<VarPtm> varPtmList = new ArrayList<>();

        Entry(int scanNum, String scanName, int precursorCharge, double precursorMass
                , String labelling, String peptide, double theoMass, int isDecoy
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.peptideSet = peptideSet;
        }
    }
}
