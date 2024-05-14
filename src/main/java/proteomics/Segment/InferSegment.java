package proteomics.Segment;


import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.SparseBooleanVector;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Types.ExpAa;
import proteomics.Types.ExpTag;
import proteomics.Types.Segment;
import proteomics.Types.VarPtm;

import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.MIN_AA_INTENSITY;

public class InferSegment {
    private static final Logger logger = LoggerFactory.getLogger(InferSegment.class);
    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;
    public static final byte N_TAG = -1;
    public static final byte NON_NC_TAG = 0;
    public static final byte C_TAG = 1;
    private static final Pattern pattern = Pattern.compile("([nc][0-9a-i])?([A-Z#$].?)");
    private final double ms2Tolerance;
    private final TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private final Map<Double, String> augedMassAaMap = new HashMap<>(35, 1);
    private double maxAugedMass = 0;
    private final Double[] augedMassArray;
    private final Map<String, Double> extraAaMassMap = new HashMap<>(35, 1);
    private final double[] nTermPossibleMod = null;
    private final double[] cTermPossibleMod = null;
    private final MassTool massTool;
    public String aaModLabel = "~!@#$%^&*";
    public String nModLabel = "qwertyuio";
    public String cModLabel = "asdfghjkl";

    public InferSegment(MassTool massTool, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        this.massTool = massTool;
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        Map<Character, Double> massTable = massTool.getMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};

        Map<Double, Character> oriMassAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            oriMassAaMap.put(massTable.get(aa), aa);
        }

        Character[] aaArray = oriMassAaMap.values().toArray(new Character[0]);

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    aaVectorTemplate.put(new Segment(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3)), 0);
                }
            }
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }

        for (double k : oriMassAaMap.keySet()) {
            augedMassAaMap.put(k, oriMassAaMap.get(k).toString());
        }
        int iLabelAa = 0;
        int iLabelN = 0;
        int iLabelC = 0;
        parameterMap.put("modCommonForTag", "114.042927,K,4,GG");
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;

            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;

            String[] modStr = v.split(",");
            double modMass = Double.valueOf(modStr[0]);
            char modSite = modStr[1].charAt(0);
            char label;
            byte modPosition = Byte.valueOf(modStr[2]);
            VarPtm varPtmByUser = new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", 0);
            if (modPosition == 4) {
                boolean added = false;
                label = aaModLabel.charAt(iLabelAa);
                double augedMass = massTable.get(modSite) + modMass;
                boolean shouldAdd = true;
                for (double tempMass : augedMassAaMap.keySet()) {
                    if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                        logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, modSite, augedMassAaMap.get(tempMass), augedMass, tempMass));
                        shouldAdd = false;
                    }
                }
                if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                    String augedAa = "" + modSite + label;
                    augedMassAaMap.put(augedMass, augedAa);
                    extraAaMassMap.put(augedAa, modMass);
                    added = true;
                }

                if (added) {
                    iLabelAa++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }

            } else if (modPosition == 0 || modPosition == 2) {
                label = nModLabel.charAt(iLabelN);
                boolean added = false;
                if (modSite == 'X') {
                    for (char oriAa : aaArray) {
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            String augedAa = "" + oriAa + label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else {
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        String augedAa = "" + modSite + label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }

                if (added) {
                    iLabelN++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            } else {
                label = cModLabel.charAt(iLabelC);
                boolean added = false;
                if (modSite == 'X') {
                    for (char oriAa : aaArray) {
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            String augedAa = "" + oriAa + label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else {
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        String augedAa = "" + modSite + label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }
                if (added) {
                    iLabelC++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            }
        }
        augedMassArray = augedMassAaMap.keySet().toArray(new Double[0]);
        maxAugedMass = Collections.max(augedMassAaMap.keySet()) + ms2Tolerance;
    }


    private final int minTagLen = 3;
    private final int maxTagLen = 6;


    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {

        String normalizedPeptide = normalizeSequence(DbTool.getSequenceOnly(peptide));

        Set<Integer> tempSet = new HashSet<>(DbTool.getSequenceOnly(peptide).length() + 1, 1);
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
        }
        return new SparseBooleanVector(tempSet);
    }

    public static String normalizeSequence(String seq) {
        return seq.replaceAll("I", "L");
    }

    public List<ExpTag> getLongTag(TreeMap<Double, Double> plMap, int minTagLenToExtractLocal) {
        int maxTagLenToExtractLocal = 99;
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        List<ExpTag> outputList = new LinkedList<>();

        Set<Pair<Integer, Integer>> edgeSet = new HashSet<>();
        Set<Integer> nodeSet = new HashSet<>();
        Set<Integer> startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Set<Integer> endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Map<Pair<Integer, Integer>, ExpAa> edgeInfoMap = new HashMap<>();
        for (int i = 0; i < mzArray.length - 1; ++i) {
            double mz1 = mzArray[i];
            int isNorC = NON_NC_TAG;
            if (Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance) {
                isNorC = N_TAG;
            } else if (Math.abs(mz1 - MassTool.PROTON - massTool.H2O) <= ms2Tolerance) {
                isNorC = C_TAG;
            }
            double intensity1 = intensityArray[i];
            if (intensity1 < 0.1) continue;
            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                if (mz2 > mz1 + maxAugedMass + 5) break;
                double intensity2 = intensityArray[j];
                if (intensity2 < 0.1) continue;
                if ((intensity1 + intensity2) < MIN_AA_INTENSITY) continue;
                String aa = inferAA(mz1, mz2, isNorC);

                if (aa != null) {
                    Pair<Integer, Integer> edge = new Pair<>(i, j);
                    edgeSet.add(edge);
                    nodeSet.add(i);
                    nodeSet.add(j);
                    startNodeSet.remove(j);
                    endNodeSet.remove(i);
                    edgeInfoMap.put(edge, new ExpAa(aa, aa.charAt(0), mz1, mz2, intensity1, intensity2, isNorC));
                }
            }
        }
        if (edgeInfoMap.size() > 120) {
            List<Map.Entry<Pair<Integer, Integer>, ExpAa>> edgeInfoList = new ArrayList<>(edgeInfoMap.entrySet());
            Collections.sort(edgeInfoList, Comparator.comparing(o -> o.getValue().getTotalIntensity(), Comparator.reverseOrder()));
            nodeSet.clear();
            edgeSet.clear();
            startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            int i = 0;
            for (Map.Entry<Pair<Integer, Integer>, ExpAa> entry : edgeInfoList) {
                if (i > 120) {
                    edgeInfoMap.remove(entry.getKey());
                } else {
                    edgeSet.add(entry.getKey());
                    nodeSet.add(entry.getKey().getFirst());
                    nodeSet.add(entry.getKey().getSecond());
                    startNodeSet.remove(entry.getKey().getSecond());
                    endNodeSet.remove(entry.getKey().getFirst());
                }
                i++;
            }

        }

        startNodeSet.retainAll(nodeSet);
        endNodeSet.retainAll(nodeSet);
        Graph g = new Graph(edgeSet, nodeSet);
        ArrayList<ArrayList<Integer>> allPath = g.getAllPaths(startNodeSet, endNodeSet);
        for (ArrayList<Integer> path : allPath) {
            if (path.size() < minTagLenToExtractLocal + 1) continue;
            if (path.size() > maxTagLenToExtractLocal + 1) {
                for (int i = 0; i < path.size() - maxTagLenToExtractLocal; i++) {
                    List<ExpAa> expAaList = new ArrayList<>();
                    for (int ii = i; ii < i + maxTagLenToExtractLocal; ii++) {
                        int j = ii + 1;
                        expAaList.add(edgeInfoMap.get(new Pair<>(path.get(ii), path.get(j))));
                    }
                    outputList.add(new ExpTag(expAaList));
                }
            } else {
                List<ExpAa> expAaList = new ArrayList<>();
                for (int i = 0; i < path.size() - 1; i++) {
                    int j = i + 1;
                    expAaList.add(edgeInfoMap.get(new Pair<>(path.get(i), path.get(j))));
                }
                outputList.add(new ExpTag(expAaList));
            }
        }
        outputList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
        boolean[] shouldKeep = new boolean[outputList.size()];
        Set<String> addedTags = new HashSet<>();
        int i = 0;
        for (ExpTag tag : outputList) {
            if (addedTags.contains(tag.getFreeAaString())) {
                shouldKeep[i] = false;
            } else {
                shouldKeep[i] = true;
                addedTags.add(tag.getFreeAaString());
            }
            i++;
        }
        List<ExpTag> finalList = new LinkedList<>();
        for (int j = 0; j < outputList.size(); j++) {
            if (shouldKeep[j]) {
                finalList.add(outputList.get(j));
            }
        }
        return finalList;
    }

    public int getLongTagNumForProt(String prot) {
        String normalizedProt = normalizeSequence(DbTool.getSequenceOnly(prot));
        Set<String> tagSet = new HashSet<>();
        for (int i = 0; i <= normalizedProt.length() - 7; ++i) {
            Segment seg = new Segment(normalizedProt.substring(i, i + 7));
            String tag = seg.toString();
            tagSet.add(tag);
        }

        return tagSet.size();
    }

    private String inferAA(double mz1, double mz2, int isNorC) {
        double mzDiff = mz2 - mz1;
        String aa = null;
        for (Map.Entry<Double, String> massAa : augedMassAaMap.entrySet()) {
            String augedAa = massAa.getValue();
            if ((isNorC == N_TAG && cModLabel.contains(augedAa.substring(augedAa.length() - 1))) || (isNorC == C_TAG && nModLabel.contains(augedAa.substring(augedAa.length() - 1))) || (isNorC == NON_NC_TAG && cModLabel.contains(augedAa.substring(augedAa.length() - 1))) || (isNorC == NON_NC_TAG && nModLabel.contains(augedAa.substring(augedAa.length() - 1)))) {
                continue;
            }

            if (Math.abs(mzDiff - massAa.getKey()) <= 1.5 * ms2Tolerance) {
                aa = massAa.getValue();
            }
        }

        if (aa != null && isNorC == C_TAG && (!"KR".contains(aa))) {
            aa = null;
        }

        if (aa != null && aa.length() > 1) {
            int a = 1;
        }
        return aa;
    }

    public TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {
            }

        }

        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {
        }
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }
        finalPlMap.put(MassTool.PROTON + massTool.H2O, 1d);


        return finalPlMap;
    }
}
