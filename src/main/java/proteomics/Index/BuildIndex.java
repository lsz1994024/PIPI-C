package proteomics.Index;

import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import proteomics.FM.FMIndex;
import proteomics.FM.UnicodeReader;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


public class BuildIndex {
    public final MassTool massTool;
    private final Map<Character, Double> fixModMap = new HashMap<>(25, 1);
    public double minPeptideMass = 9999;
    public double maxPeptideMass = 0;
    public InferSegment inferSegment;
    public TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();
    private final String labelling;
    public final DbTool dbTool;
    private final InferPTM inferPTM;
    public Map<String, Integer> protLengthMap = new HashMap<>();
    public FMIndex fmIndexFull;
    public int[] dotPosArrFull;
    public Map<Integer, String> posProtMapFull = new HashMap<>();
    public FMIndex fmIndexNormal;
    public int[] dotPosArrNormal;
    public Map<Integer, String> posProtMapNormal = new HashMap<>();
    public FMIndex fmIndexReverse;
    public int textNormalLen;
    public final boolean nTermSpecific;
    public final boolean cTermSpecific;
    public Map<String, String> protSeqMap;

    public BuildIndex(Map<String, String> parameterMap) throws Exception {
        boolean addContaminant = parameterMap.get("add_contaminant").contentEquals("1");
        String dbPath = parameterMap.get("database");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        nTermSpecific = Integer.valueOf(parameterMap.get("n_specific")) == 1;
        cTermSpecific = Integer.valueOf(parameterMap.get("c_specific")) == 1;
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = parameterMap.get("15N").trim().contentEquals("1") ? "N15" : "N14";
        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('X', 0.0);
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Double.valueOf(parameterMap.get("c")));

        dbTool = new DbTool(dbPath, parameterMap.get("database_type"));
        DbTool contaminantsDb = null;
        if (addContaminant) {
            contaminantsDb = new DbTool(null, "contaminants");
            protSeqMap = contaminantsDb.getProtSeqMap();
            protSeqMap.putAll(dbTool.getProtSeqMap());
        } else {
            protSeqMap = dbTool.getProtSeqMap();
        }
        massTool = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site_1").trim(), parameterMap.get("protection_site_1").trim(), parameterMap.get("is_from_C_term_1").trim().contentEquals("1"), parameterMap.getOrDefault("cleavage_site_2", null), parameterMap.getOrDefault("protection_site_2", null), parameterMap.containsKey("is_from_C_term_2") ? parameterMap.get("is_from_C_term_2").trim().contentEquals("1") : null, ms2Tolerance, oneMinusBinOffset, labelling);

        inferPTM = new InferPTM(massTool, fixModMap, parameterMap);

        inferSegment = new InferSegment(massTool, parameterMap, fixModMap);

        BufferedWriter writerProt = new BufferedWriter(new FileWriter("catProt.txt"));
        int dotPos = 0;
        int dotNum = 0;
        dotPosArrFull = new int[protSeqMap.keySet().size()];

        for (String protId : protSeqMap.keySet()) {
            dotPosArrFull[dotNum] = dotPos;
            posProtMapFull.put(dotNum, protId);
            String protSeq = protSeqMap.get(protId).replace('I', 'L');
            protSeqMap.put(protId, protSeq);
            writerProt.write("." + protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length() + 1;
            int numOfTags = inferSegment.getLongTagNumForProt(protSeq);
            protLengthMap.put(protId, numOfTags);
        }
        writerProt.close();
        char[] text = loadFile("catProt.txt", true);
        fmIndexFull = new FMIndex(text);
    }

    public static char[] loadFile(String file, boolean appendTerminalCharacter) throws IOException {
        BufferedReader reader = null;
        CharArrayWriter writer = null;
        UnicodeReader r = new UnicodeReader(new FileInputStream(file), null);

        char[] buffer = new char[20 * 1024];
        int read;
        try {
            reader = new BufferedReader(r);
            writer = new CharArrayWriter();
            while ((read = reader.read(buffer)) != -1) {
                writer.write(buffer, 0, read);
            }
            if (appendTerminalCharacter) {
                writer.append('\0');
            }
            writer.flush();
            return writer.toCharArray();
        } catch (IOException ex) {
            throw ex;
        } finally {
            try {
                writer.close();
                reader.close();
                r.close();
            } catch (Exception ex) {
            }
        }
    }

    public MassTool returnMassTool() {
        return massTool;
    }

    public Map<Character, Double> returnFixModMap() {
        return fixModMap;
    }

    public InferSegment getInferSegment() {
        return inferSegment;
    }

    public InferPTM getInferPTM() {
        return inferPTM;
    }

    public String getLabelling() {
        return labelling;
    }

}
