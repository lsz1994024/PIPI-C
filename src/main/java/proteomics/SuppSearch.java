package proteomics;

import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import proteomics.Index.BuildIndex;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.Peptide;
import proteomics.Types.PeptideInfo;
import proteomics.Types.PosMassMap;
import proteomics.Types.VarPtm;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.lszDebugScanNum;

public class SuppSearch implements Callable<SuppSearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms2Tolerance;
    private final JMzReader spectraParser;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor preSpectrum;
    private final Binomial binomial;
    private final int scanNum;
    private final Map<String, TreeMap<Integer, VarPtm>> ptmSeq_posVarPtmMap_Map;
    private Map<String, PeptideInfo> peptideInfoMap;

    private boolean isHomo(Peptide p1, Peptide p2) {

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

    public SuppSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, JMzReader spectraParser, ReentrantLock lock, String scanName, int precursorCharge
            , double precursorMass, SpecProcessor preSpectrum, Binomial binomial, Map<String, TreeMap<Integer, VarPtm>> ptmSeq_posVarPtmMap_Map, Map<String, PeptideInfo> peptideInfoMap) {
        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms2Tolerance = ms2Tolerance;
        this.spectraParser = spectraParser;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.preSpectrum = preSpectrum;
        this.binomial = binomial;
        this.scanNum = scanNum;
        this.ptmSeq_posVarPtmMap_Map = ptmSeq_posVarPtmMap_Map;
        this.peptideInfoMap = peptideInfoMap;
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

        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
        } else {
            expProcessedPL = preSpectrum.digitizePL(plMap);
        }


        TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
        Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();

        for (String pepPtmSeq : ptmSeq_posVarPtmMap_Map.keySet()) {
            String pepFreeSeq = pepPtmSeq.replaceAll("[(\\[][^A-Znc()\\[\\]]+[)\\]]", "");
            TreeMap<Integer, VarPtm> refVarPtmMap = ptmSeq_posVarPtmMap_Map.get(pepPtmSeq);

            PosMassMap fullPosMassMap = new PosMassMap();
            for (int pos : refVarPtmMap.keySet()) {
                fullPosMassMap.put(pos, refVarPtmMap.get(pos).mass);
            }
            Peptide peptide = new Peptide(pepFreeSeq, true, massTool);
            if (!fullPosMassMap.isEmpty()) {
                peptide.setVarPTM(fullPosMassMap);
            }

            double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions);
            peptide.setScore(score * (1 - fullPosMassMap.size() * 0.02));
            peptide.posVarPtmResMap = refVarPtmMap;

            if (score > 0) {
                if (peptideSet.size() < candisNum) {
                    peptideSet.add(peptide);
                } else if (peptide.compareTo(peptideSet.last()) == 1) {
                    peptideSet.pollLast();
                    peptideSet.add(peptide);
                }
                TreeSet<Peptide> peptideTreeSet = new TreeSet<>();
                peptideTreeSet.add(peptide);
                modSequences.put(pepFreeSeq, peptideTreeSet);
            }
        }

        if (peptideSet.isEmpty()) {
            return null;
        }

        List<Peptide> pepList = new ArrayList<>(peptideSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size() - 1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;

            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;

                if (isHomo(pepList.get(i), pepList.get(j))
                        || pepList.get(j).getScore() > 0.6 * pepList.get(i).getScore()
                ) {
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
            }
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++) {
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            return null;
        }


        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);

        Peptide topPep = peptideArray[0];
        Entry entry;
        TreeSet<Peptide> ptmPatterns = null;
        if (topPep.hasVarPTM()) {
            ptmPatterns = modSequences.get(topPep.getFreeSeq());
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
        if (ptmPatterns != null) {
            List<String> tempList = new LinkedList<>();
            Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
            ptmPatternsIterator.next();
            while (ptmPatternsIterator.hasNext()) {
                Peptide temp = ptmPatternsIterator.next();
                tempList.add(String.format(Locale.US, "%s-%.4f", temp.getPtmContainingSeq(buildIndex.returnFixModMap()), temp.getScore()));
            }
            otherPtmPatterns = String.join(";", tempList);
        }
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
        entry = new Entry(
                scanNum, scanName
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length() - 1)
        );

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        return entry;
    }

    public class Entry {

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
        final String candidates;
        final String peptideSet;

        Entry(int scanNum, String scanName, int precursorCharge, double precursorMass
                , String labelling, String peptide, double theoMass, int isDecoy, double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet) {
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
            this.candidates = candidates;
            this.peptideSet = peptideSet;
        }
    }
}
