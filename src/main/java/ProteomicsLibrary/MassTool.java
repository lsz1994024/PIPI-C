package ProteomicsLibrary;

import ProteomicsLibrary.Types.AA;
import ProteomicsLibrary.Types.SparseVector;
import proteomics.Types.Peptide;
import proteomics.Types.PosMassMap;
import proteomics.Types.VarPtm;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MassTool {

    private static final Pattern modAAPattern = Pattern.compile("([A-Znc])(([(\\[])?([0-9.\\-+]+)([)\\]])?)?");
    private static final Pattern leftFlankPattern = Pattern.compile("^[A-Z-]\\.");
    private static final Pattern rightFlankPattern = Pattern.compile("\\.[A-Z-]$");

    public static final double PROTON = 1.00727646688;
    public static final double C13_DIFF = 1.00335483;
    public static final double N1514_DIFF = 15.0001088 - 14.0030732;
    public Map<Character, VarPtm> labelVarPtmMap = new HashMap<>();
    public static final Map<Character, Double> n1514DiffMap = new HashMap<>();

    static {
        n1514DiffMap.put('G', N1514_DIFF);
        n1514DiffMap.put('A', N1514_DIFF);
        n1514DiffMap.put('S', N1514_DIFF);
        n1514DiffMap.put('P', N1514_DIFF);
        n1514DiffMap.put('V', N1514_DIFF);
        n1514DiffMap.put('T', N1514_DIFF);
        n1514DiffMap.put('C', N1514_DIFF);
        n1514DiffMap.put('I', N1514_DIFF);
        n1514DiffMap.put('L', N1514_DIFF);
        n1514DiffMap.put('N', N1514_DIFF * 2);
        n1514DiffMap.put('D', N1514_DIFF);
        n1514DiffMap.put('Q', N1514_DIFF * 2);
        n1514DiffMap.put('K', N1514_DIFF * 2);
        n1514DiffMap.put('E', N1514_DIFF);
        n1514DiffMap.put('M', N1514_DIFF);
        n1514DiffMap.put('H', N1514_DIFF * 3);
        n1514DiffMap.put('F', N1514_DIFF);
        n1514DiffMap.put('R', N1514_DIFF * 4);
        n1514DiffMap.put('Y', N1514_DIFF);
        n1514DiffMap.put('W', N1514_DIFF * 2);
        n1514DiffMap.put('U', N1514_DIFF);
        n1514DiffMap.put('O', N1514_DIFF * 3);
        n1514DiffMap.put('B', 0d);
        n1514DiffMap.put('J', 0d);
        n1514DiffMap.put('X', 0d);
        n1514DiffMap.put('Z', 0d);
        n1514DiffMap.put('n', 0d);
        n1514DiffMap.put('c', 0d);
        n1514DiffMap.put('*', 0d);
    }

    public final double H2O;

    private final Map<String, Double> elementTable = new HashMap<>();
    private final Map<Character, Double> massTable = new HashMap<>(30, 1);
    public final int missedCleavage;
    private final double ms2Tolerance;
    private final double inverse2Ms2Tolerance;
    private final double oneMinusBinOffset;
    private final Pattern digestSitePattern1;
    private final Pattern digestSitePattern2;
    private final Boolean cleavageFromCTerm1;
    private final Boolean cleavageFromCTerm2;
    private final Pattern digestSitePatternForLinkSiteChecking1;
    private final Pattern digestSitePatternForLinkSiteChecking2;
    private final String labelling;
    private final Map<Character, Double> fixModMap;

    public MassTool(int missedCleavage, Map<Character, Double> fixModMap, String cleavageSite1, String protectionSite1, boolean cleavageFromCTerm1, String cleavageSite2, String protectionSite2, Boolean cleavageFromCTerm2, double ms2Tolerance, double oneMinusBinOffset, String labelling) {
        this.labelling = labelling;
        this.fixModMap = fixModMap;

        elementTable.put("-", 0d);
        elementTable.put("H", 1.0078246);
        elementTable.put("He", 3.01603);
        elementTable.put("Li", 6.015121);
        elementTable.put("Be", 9.012182);
        elementTable.put("B", 10.012937);
        elementTable.put("C", 12.0000000);
        elementTable.put("N", 14.0030732);
        if (labelling.contentEquals("N15")) {
            elementTable.put("N", 15.0001088);
        }
        elementTable.put("O", 15.9949141);
        elementTable.put("F", 18.9984032);
        elementTable.put("Ne", 19.992435);
        elementTable.put("Na", 22.989767);
        elementTable.put("Mg", 23.985042);
        elementTable.put("Al", 26.981539);
        elementTable.put("Si", 27.976927);
        elementTable.put("P", 30.973762);
        elementTable.put("S", 31.972070);
        elementTable.put("Cl", 34.9688531);
        elementTable.put("Ar", 35.967545);
        elementTable.put("K", 38.963707);
        elementTable.put("Ca", 39.962591);
        elementTable.put("Sc", 44.955910);
        elementTable.put("Ti", 45.952629);
        elementTable.put("V", 49.947161);
        elementTable.put("Cr", 49.946046);
        elementTable.put("Mn", 54.938047);
        elementTable.put("Fe", 53.939612);
        elementTable.put("Co", 58.933198);
        elementTable.put("Ni", 57.935346);
        elementTable.put("Cu", 62.939598);
        elementTable.put("Zn", 63.929145);
        elementTable.put("Ga", 68.925580);
        elementTable.put("Ge", 69.924250);
        elementTable.put("As", 74.921594);
        elementTable.put("Se", 73.922475);
        elementTable.put("Br", 78.918336);
        elementTable.put("Kr", 77.914);
        elementTable.put("Rb", 84.911794);
        elementTable.put("Sr", 83.913430);
        elementTable.put("Y", 88.905849);
        elementTable.put("Zr", 89.904703);
        elementTable.put("Nb", 92.906377);
        elementTable.put("Mo", 91.906808);
        elementTable.put("Tc", 98.0);
        elementTable.put("Ru", 95.907599);
        elementTable.put("Rh", 102.905500);
        elementTable.put("Pd", 101.905634);
        elementTable.put("Ag", 106.905092);
        elementTable.put("Cd", 105.906461);
        elementTable.put("In", 112.904061);
        elementTable.put("Sn", 111.904826);
        elementTable.put("Sb", 120.903821);
        elementTable.put("Te", 119.904048);
        elementTable.put("I", 126.904473);
        elementTable.put("Xe", 123.905894);
        elementTable.put("Cs", 132.905429);
        elementTable.put("Ba", 129.906282);
        elementTable.put("La", 137.90711);
        elementTable.put("Ce", 135.907140);
        elementTable.put("Pr", 140.907647);
        elementTable.put("Nd", 141.907719);
        elementTable.put("Pm", 145.0);
        elementTable.put("Sm", 143.911998);
        elementTable.put("Eu", 150.919847);
        elementTable.put("Gd", 151.919786);
        elementTable.put("Tb", 158.925342);
        elementTable.put("Dy", 155.925277);
        elementTable.put("Ho", 164.930319);
        elementTable.put("Er", 161.928775);
        elementTable.put("Tm", 168.934212);
        elementTable.put("Yb", 167.933894);
        elementTable.put("Lu", 174.940770);
        elementTable.put("Hf", 173.940044);
        elementTable.put("Ta", 179.947462);
        elementTable.put("W", 179.946701);
        elementTable.put("Re", 184.952951);
        elementTable.put("Os", 183.952488);
        elementTable.put("Ir", 190.960584);
        elementTable.put("Pt", 189.959917);
        elementTable.put("Au", 196.966543);
        elementTable.put("Hg", 195.965807);
        elementTable.put("Tl", 202.972320);
        elementTable.put("Pb", 203.973020);
        elementTable.put("Bi", 208.980374);
        elementTable.put("Po", 209.0);
        elementTable.put("At", 210.0);
        elementTable.put("Rn", 222.0);
        elementTable.put("Fr", 223.0);
        elementTable.put("Ra", 226.025);
        elementTable.put("Th", 232.038054);
        elementTable.put("Pa", 231.0359);
        elementTable.put("U", 234.040946);
        elementTable.put("Np", 237.048);
        elementTable.put("Pu", 244.0);
        elementTable.put("Am", 243.0);
        elementTable.put("Cm", 247.0);
        elementTable.put("Bk", 247.0);
        elementTable.put("Cf", 251.0);
        elementTable.put("Es", 252.0);
        elementTable.put("Fm", 257.0);
        elementTable.put("Md", 258.0);
        elementTable.put("No", 259.0);
        elementTable.put("Lr", 260.0);
        elementTable.put("13C", 13.0033554);
        elementTable.put("15N", 15.0001088);
        elementTable.put("18O", 17.9991616);
        elementTable.put("2H", 2.0141021);
        elementTable.put("dHex", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 10);
        elementTable.put("Hep", elementTable.get("C") * 7 + elementTable.get("O") * 6 + elementTable.get("H") * 12);
        elementTable.put("Hex", elementTable.get("C") * 6 + elementTable.get("O") * 5 + elementTable.get("H") * 10);
        elementTable.put("HexA", elementTable.get("C") * 6 + elementTable.get("O") * 6 + elementTable.get("H") * 8);
        elementTable.put("HexN", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 11 + elementTable.get("N"));
        elementTable.put("HexNAc", elementTable.get("C") * 8 + elementTable.get("O") * 5 + +elementTable.get("N") + elementTable.get("H") * 13);
        elementTable.put("Kdn", elementTable.get("C") * 9 + elementTable.get("H") * 14 + elementTable.get("O") * 8);
        elementTable.put("Kdo", elementTable.get("C") * 8 + elementTable.get("H") * 12 + elementTable.get("O") * 7);
        elementTable.put("NeuAc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 8 + elementTable.get("N"));
        elementTable.put("NeuGc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 9 + elementTable.get("N"));
        elementTable.put("Pent", elementTable.get("C") * 5 + elementTable.get("O") * 4 + elementTable.get("H") * 8);
        elementTable.put("Phos", elementTable.get("O") * 3 + elementTable.get("H") + elementTable.get("P"));
        elementTable.put("Sulf", elementTable.get("S") + elementTable.get("O") * 3);
        elementTable.put("Water", elementTable.get("H") * 2 + elementTable.get("O"));
        elementTable.put("Me", elementTable.get("C") + elementTable.get("H") * 2);
        elementTable.put("Ac", elementTable.get("C") * 2 + elementTable.get("H") * 2 + elementTable.get("O"));
        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        inverse2Ms2Tolerance = 1 / (2 * ms2Tolerance);
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm1 = cleavageFromCTerm1;
        this.cleavageFromCTerm2 = cleavageFromCTerm2;
        massTable.put('G', (elementTable.get("C") * 2 + elementTable.get("H") * 3 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('G')));
        massTable.put('A', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('A')));
        massTable.put('S', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('S')));
        massTable.put('P', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('P')));
        massTable.put('V', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('V')));
        massTable.put('T', (elementTable.get("C") * 4 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('I')));
        massTable.put('C', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('C')));
        massTable.put('I', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('I')));
        massTable.put('L', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('L')));
        massTable.put('N', (elementTable.get("C") * 4 + elementTable.get("H") * 6 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('N')));
        massTable.put('D', (elementTable.get("C") * 4 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('D')));
        massTable.put('Q', (elementTable.get("C") * 5 + elementTable.get("H") * 8 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('Q')));
        massTable.put('K', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('K')));
        massTable.put('E', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('E')));
        massTable.put('M', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('M')));
        massTable.put('H', (elementTable.get("C") * 6 + elementTable.get("H") * 7 + elementTable.get("N") * 3 + elementTable.get("O") + fixModMap.get('H')));
        massTable.put('F', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('F')));
        massTable.put('R', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 4 + elementTable.get("O") + fixModMap.get('R')));
        massTable.put('Y', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('Y')));
        massTable.put('W', (elementTable.get("C") * 11 + elementTable.get("H") * 10 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('W')));
        massTable.put('U', (elementTable.get("C") * 3 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + elementTable.get("Se") + fixModMap.get('U')));
        massTable.put('O', (elementTable.get("C") * 12 + elementTable.get("H") * 21 + elementTable.get("N") * 3 + elementTable.get("O") * 3 + fixModMap.get('O')));
        massTable.put('n', fixModMap.get('n'));
        massTable.put('c', fixModMap.get('c'));
        massTable.put('$', (massTable.get('Q') + massTable.get('K')) * 0.5);
        H2O = elementTable.get("H") * 2 + elementTable.get("O");
        massTable.put('B', 0d);
        massTable.put('J', 0d);
        massTable.put('X', 0d);
        massTable.put('Z', 0d);
        massTable.put('*', 0d);

        digestSitePattern1 = getDigestSitePattern(cleavageSite1, protectionSite1, cleavageFromCTerm1);
        digestSitePatternForLinkSiteChecking1 = getDigestSitePatternForLinkSiteChecking(cleavageSite1, protectionSite1, cleavageFromCTerm1);
        if (cleavageSite2 != null) {
            digestSitePattern2 = getDigestSitePattern(cleavageSite2, protectionSite2, cleavageFromCTerm2);
            digestSitePatternForLinkSiteChecking2 = getDigestSitePatternForLinkSiteChecking(cleavageSite2, protectionSite2, cleavageFromCTerm2);
        } else {
            digestSitePattern2 = null;
            digestSitePatternForLinkSiteChecking2 = null;
        }
    }

    public double calResidueMass(String sequence) {
        double totalMass = 0;
        Matcher matcher = modAAPattern.matcher(sequence);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double deltaMass = 0;
            if (matcher.group(4) != null) {
                deltaMass = Double.valueOf(matcher.group(4));
            }
            totalMass += massTable.get(aa) + deltaMass;
        }

        return totalMass;
    }


    public static AA[] seqToAAList(String sequence) {
        Matcher matcher = modAAPattern.matcher(sequence);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double deltaMass = 0;
            if (matcher.group(4) != null) {
                deltaMass = Double.valueOf(matcher.group(4));
            }
            temp.add(new AA(aa, deltaMass));
        }
        return temp.toArray(new AA[0]);
    }


    public double[][] buildIonArray(String sequence) {
        AA[] aaArray = seqToAAList(sequence);

        double[][] ionArray = new double[2][aaArray.length];
        double bIonMass = PROTON;
        for (int i = 0; i < aaArray.length; ++i) {
            bIonMass += massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            ionArray[0][i] = bIonMass;
        }

        double yIonMass = H2O + PROTON;
        for (int i = aaArray.length - 1; i >= 0; --i) {
            yIonMass += massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            ionArray[1][i] = yIonMass;
        }
        return ionArray;
    }

    public double buildVectorAndCalXCorr(double[][] ionMatrix, int precursorCharge, SparseVector xcorrPL, Map<Integer, Double> matchedBions, Map<Integer, Double> matchedYions) {
        int colNum = ionMatrix[0].length;
        int rowNum = Math.min(ionMatrix.length / 2, precursorCharge - 1) * 2;
        if (precursorCharge == 1) {
            rowNum = 2;
        }
        double xcorr = 0;
        for (int i = 0; i < rowNum; ++i) {
            for (int j = 0; j < colNum; ++j) {
                double peakIntensity = xcorrPL.get(mzToBin(ionMatrix[i][j]));
                if (peakIntensity > 0.0) {
                    if (0 == i) {
                        matchedBions.put(j, peakIntensity);
                    } else if (1 == i) {
                        matchedYions.put(j, peakIntensity);
                    }
                }
                xcorr += peakIntensity;
            }
        }
        return xcorr * 0.25;
    }

    public double buildVectorAndCalXCorr(double[][] ionMatrix, int precursorCharge, SparseVector xcorrPL, Map<Integer, Double> matchedBions, Map<Integer, Double> matchedYions, Set<Integer> jRange) {
        int colNum = ionMatrix[0].length;
        int rowNum = Math.min(ionMatrix.length / 2, precursorCharge - 1) * 2;
        if (precursorCharge == 1) {
            rowNum = 2;
        }
        double xcorr = 0;
        for (int i = 0; i < rowNum; ++i) {
            for (int j : jRange) {
                double peakIntensity = xcorrPL.get(mzToBin(ionMatrix[i][j]));
                if (peakIntensity > 0.0) {
                    if (0 == i) {
                        matchedBions.put(j, peakIntensity);
                    } else if (1 == i) {
                        matchedYions.put(j, peakIntensity);
                    }
                }
                xcorr += peakIntensity;
            }
        }
        return xcorr * 0.25;
    }


    public Map<Character, Double> getMassTable() {
        return massTable;
    }

    public Map<String, Double> getElementTable() {
        return elementTable;
    }

    public int mzToBin(double mz) {
        return (int) Math.floor(mz * inverse2Ms2Tolerance + oneMinusBinOffset);
    }

    public String getLabelling() {
        return labelling;
    }

    static Pattern getDigestSitePattern(String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
        Pattern digestSitePattern;
        if (cleavageFromCTerm) {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
            } else {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "](?![" + protectionSite + "])");
            }
        } else {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
            } else {
                digestSitePattern = Pattern.compile("(?<![" + protectionSite + "])" + "[" + cleavageSite + "]");
            }
        }
        return digestSitePattern;
    }

    private static Pattern getDigestSitePatternForLinkSiteChecking(String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
        Pattern digestSitePattern2;
        if (cleavageFromCTerm) {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern2 = Pattern.compile("[" + cleavageSite + "]$");
            } else {
                digestSitePattern2 = Pattern.compile("[" + cleavageSite + "](?![" + protectionSite + "])$");
            }
        } else {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern2 = Pattern.compile("^[" + cleavageSite + "]");
            } else {
                digestSitePattern2 = Pattern.compile("^(?<![" + protectionSite + "])" + "[" + cleavageSite + "]");
            }
        }
        return digestSitePattern2;
    }

}
