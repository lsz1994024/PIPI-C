package ProteomicsLibrary;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DbTool {
    private static final Map<Character, Character> nextAaMap = new HashMap<>();
    private static final Random random = new Random();

    private final Map<String, String> protSeqMap = new HashMap<>();
    private final Map<String, String> proteinAnnotateMap = new HashMap<>();

    public DbTool(String dbName, String databaseType) throws IOException {
        nextAaMap.put('A', 'C');
        nextAaMap.put('C', 'D');
        nextAaMap.put('D', 'E');
        nextAaMap.put('E', 'F');
        nextAaMap.put('F', 'G');
        nextAaMap.put('G', 'H');
        nextAaMap.put('H', 'L');
        nextAaMap.put('L', 'M');
        nextAaMap.put('M', 'N');
        nextAaMap.put('N', 'P');
        nextAaMap.put('P', 'Q');
        nextAaMap.put('Q', 'R');
        nextAaMap.put('R', 'S');
        nextAaMap.put('S', 'T');
        nextAaMap.put('T', 'V');
        nextAaMap.put('V', 'W');
        nextAaMap.put('W', 'Y');
        nextAaMap.put('Y', 'A');


        String id = "";
        String annotate;
        StringBuilder sequence = new StringBuilder(99999);
        databaseType = databaseType.trim().toLowerCase();

        boolean newPro = true;

        Pattern headerPattern;
        if (databaseType.contentEquals("tair")) {
            headerPattern = Pattern.compile("^>([^\\s]+)[\\s|]*(.*)$");
        } else if (databaseType.contentEquals("uniprot") || databaseType.contentEquals("swissprot")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)$");
        } else if (databaseType.contentEquals("nextprot")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)");
        } else if (databaseType.contentEquals("contaminants") || databaseType.contentEquals("itag") || databaseType.contentEquals("refseq")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)$");
        } else if (databaseType.contentEquals("others")) {
            headerPattern = Pattern.compile("^>(.+)$");
        } else {
            throw new NullPointerException(String.format(Locale.US, "Incorrect database type (%s) in the parameter file.", databaseType));
        }

        BufferedReader dbReader;
        if (databaseType.contentEquals("contaminants")) {
            InputStream inputStream = getClass().getClassLoader().getResourceAsStream("contaminants.fasta");
            dbReader = new BufferedReader(new InputStreamReader(inputStream));
        } else {
            dbReader = new BufferedReader(new FileReader(dbName));
        }
        String line;
        while ((line = dbReader.readLine()) != null) {
            line = line.trim();
            Matcher headMatcher = headerPattern.matcher(line);
            if (headMatcher.matches()) {
                if (!newPro) {
                    protSeqMap.put(id, sequence.toString());
                }
                id = headMatcher.group(1).trim();
                if (databaseType.contentEquals("others")) {
                    annotate = id;
                } else {
                    annotate = headMatcher.group(2).trim();
                }
                proteinAnnotateMap.put(id, annotate);
                newPro = true;
            } else if (!line.isEmpty()) {
                if (newPro) {
                    sequence = new StringBuilder(99999);
                    sequence.append(line);
                    newPro = false;
                } else {
                    sequence.append(line);
                }
            }
        }
        dbReader.close();

        protSeqMap.put(id, sequence.toString());
    }

    public Map<String, String> getProtSeqMap() {
        return protSeqMap;
    }

    public Map<String, String> getProteinAnnotateMap() {
        return proteinAnnotateMap;
    }


    public static String shuffleProtKeepKR(String sequence, String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
        StringBuilder seqToBeShuffled;
        if (sequence.startsWith("M")) {
            seqToBeShuffled = new StringBuilder(sequence.substring(1));
        } else {
            seqToBeShuffled = new StringBuilder(sequence);
        }

        Pattern digestSitePattern = MassTool.getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
        List<Integer> cutSiteList = new ArrayList<>();
        Matcher matcher = digestSitePattern.matcher(seqToBeShuffled);
        while (matcher.find()) {
            cutSiteList.add(matcher.start());
        }
        Collections.sort(cutSiteList);
        LinkedList<Character> aaToShuffleList = new LinkedList<>();
        LinkedList<Character> krToShuffleList = new LinkedList<>();
        for (int i = 0; i < seqToBeShuffled.length(); i++) {
            if (cutSiteList.contains(i)) {
                krToShuffleList.add(seqToBeShuffled.charAt(i));
            } else {
                aaToShuffleList.add(seqToBeShuffled.charAt(i));
            }
        }
        Collections.shuffle(krToShuffleList);
        Collections.shuffle(aaToShuffleList);
        StringBuilder shuffledSeq = new StringBuilder();
        for (int i = 0; i < seqToBeShuffled.length(); i++) {
            if (cutSiteList.contains(i)) {
                shuffledSeq.append(krToShuffleList.poll());
            } else {
                shuffledSeq.append(aaToShuffleList.poll());
            }
        }

        if (sequence.startsWith("M")) {
            return "M" + shuffledSeq;
        } else {
            return shuffledSeq.toString();
        }
    }


    public static String getSequenceOnly(String peptide) {
        return peptide.replaceAll("[^A-Z]+", "");
    }

}
