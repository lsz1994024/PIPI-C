package proteomics.Types;

import ProteomicsLibrary.MassTool;
import proteomics.Segment.InferSegment;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

public class Peptide implements Comparable<Peptide>, Cloneable {
    public double nDeltaMass = -0.99;
    public double cDeltaMass = -0.99;
    public ExpTag finderTag = null;
    public int tagPosInPep = -1;

    public int scanNum = 0;
    private final String freeSeq;
    public String ptmSeq;
    public boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massTool;
    private int hashCode;

    private PosMassMap varPtmMap = null;
    public TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
    private double theoMass = -1;
    private double[][] ionMatrix = null;
    private String varPtmContainingSeq = null;

    private String ptmContainingSeq = null;

    private double score = -1;
    private int matchedPeakNum = -1;
    private double ionFrac = -1;
    private double matchedHighestIntensityFrac = -1;
    private double explainedAaFrac = -1;
    private double qValue = -1;
    private String aScore = "-";
    public Map<Integer, Double> matchedBions = new HashMap<>();
    public Map<Integer, Double> matchedYions = new HashMap<>();
    public double precursorMass;
    public double calculatedMass;

    public Peptide(String freeSeq, boolean isDecoy, MassTool massTool) {
        this.freeSeq = freeSeq;
        this.ptmSeq = freeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferSegment.normalizeSequence(freeSeq);
        this.massTool = massTool;

        hashCode = freeSeq.hashCode();
    }


    public int getPriority() {
        int res = 0;
        for (VarPtm varPtm : posVarPtmResMap.values()) {
            res += varPtm.priority;
        }
        return res;
    }

    public double[][] getIonMatrixNow() {
        varPtmContainingSeq = getVarPtmContainingSeqNow();
        ionMatrix = massTool.buildIonArray(varPtmContainingSeq);
        theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        return ionMatrix;
    }

    public String getVarPtmContainingSeqNow() {
        if (varPtmMap != null) {
            StringBuilder sb = new StringBuilder(freeSeq.length() * 5);
            int tempIdx = varPtmMap.firstKey() + 1;
            if (tempIdx > 1) {
                sb.append(freeSeq, 0, tempIdx - 1);
            }
            int i = tempIdx - 1;
            tempIdx = varPtmMap.lastKey() + 1;
            while (i < freeSeq.length()) {
                boolean hasMod = false;
                if (tempIdx > i) {
                    for (Integer co : varPtmMap.keySet()) {
                        if (co == i) {
                            sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPtmMap.get(co)));
                            hasMod = true;
                            ++i;
                            break;
                        }
                    }
                    if (!hasMod) {
                        sb.append(freeSeq.charAt(i));
                        ++i;
                    }
                } else {
                    break;
                }
            }
            if (tempIdx < freeSeq.length()) {
                sb.append(freeSeq.substring(tempIdx));
            }
            varPtmContainingSeq = sb.toString();
        } else {
            varPtmContainingSeq = freeSeq;
        }
        return varPtmContainingSeq;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public double getTheoMass() {
        if (theoMass < 0) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        }
        return theoMass;
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode == otherPeptide.hashCode;
    }

    public Peptide clone() throws CloneNotSupportedException {
        super.clone();
        Peptide other = new Peptide(freeSeq, isDecoy, massTool);
        if (varPtmMap != null) {
            other.setVarPTM(varPtmMap.clone());
            other.setScore(score);
            other.setMatchedHighestIntensityFrac(matchedHighestIntensityFrac);
            other.setExplainedAaFrac(explainedAaFrac);
            other.setIonFrac(ionFrac);
            other.setaScore(aScore);
            other.setQValue(qValue);
            if (this.posVarPtmResMap != null) {
                other.posVarPtmResMap = new TreeMap<>(this.posVarPtmResMap);
            }
        }
        other.nDeltaMass = this.nDeltaMass;
        other.cDeltaMass = this.cDeltaMass;
        other.finderTag = null;
        if (this.finderTag != null) {
            other.finderTag = new ExpTag(this.finderTag.expAaList);
        }
        other.tagPosInPep = this.tagPosInPep;
        return other;
    }

    public int hashCode() {
        return hashCode;
    }

    public void setVarPTM(PosMassMap ptmMap) {
        this.varPtmMap = ptmMap;
        if (ptmMap != null) {
            theoMass = -1;
            ionMatrix = null;
            varPtmContainingSeq = null;
            ptmContainingSeq = null;

            String toString = freeSeq + "." + ptmMap;
            hashCode = toString.hashCode();
        }
    }

    public String toString() {
        if (varPtmMap == null) {
            return freeSeq;
        }
        return freeSeq + "." + varPtmMap;
    }

    public boolean hasVarPTM() {
        return varPtmMap != null;
    }

    private int getVarPTMNum() {
        if (hasVarPTM()) {
            return varPtmMap.size();
        } else {
            return 0;
        }
    }

    public String getFreeSeq() {
        return freeSeq;
    }

    private String getVarPtmContainingSeq() {
        if (varPtmContainingSeq == null) {
            if (varPtmMap != null) {
                StringBuilder sb = new StringBuilder(freeSeq.length() * 5);

                int tempIdx = 0;
                try {
                    tempIdx = varPtmMap.firstKey() + 1;

                } catch (Exception ex) {
                    System.out.println("error");
                }
                if (tempIdx > 1) {
                    sb.append(freeSeq, 0, tempIdx - 1);
                }
                int i = tempIdx - 1;
                tempIdx = varPtmMap.lastKey() + 1;
                while (i < freeSeq.length()) {
                    boolean hasMod = false;
                    if (tempIdx > i) {
                        for (Integer co : varPtmMap.keySet()) {
                            if (co == i) {
                                sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPtmMap.get(co)));
                                hasMod = true;
                                ++i;
                                break;
                            }
                        }
                        if (!hasMod) {
                            sb.append(freeSeq.charAt(i));
                            ++i;
                        }
                    } else {
                        break;
                    }
                }
                if (tempIdx < freeSeq.length()) {
                    sb.append(freeSeq.substring(tempIdx));
                }
                varPtmContainingSeq = sb.toString();
            } else {
                varPtmContainingSeq = freeSeq;
            }
        }

        return varPtmContainingSeq;
    }

    public String getPtmContainingSeq(Map<Character, Double> fixModMap) {
        if (ptmContainingSeq == null) {
            ptmContainingSeq = getVarPtmContainingSeq();
            for (char aa : fixModMap.keySet()) {
                if (Math.abs(fixModMap.get(aa)) > 0.01) {
                    ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format(Locale.US, "%c(%.3f)", aa, fixModMap.get(aa)));
                }
            }
        }

        return ptmContainingSeq;
    }


    public void setScore(double score) {
        this.score = score;
    }

    public void setMatchedPeakNum(int matchedPeakNum) {
        this.matchedPeakNum = matchedPeakNum;
    }

    public void setIonFrac(double ionFrac) {
        this.ionFrac = ionFrac;
    }

    public void setMatchedHighestIntensityFrac(double matchedHighestIntensityFrac) {
        this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
    }

    public void setExplainedAaFrac(double explainedAaFrac) {
        this.explainedAaFrac = explainedAaFrac;
    }

    public void setaScore(String aScore) {
        this.aScore = aScore;
    }

    private void setQValue(double qValue) {
        this.qValue = qValue;
    }

    public double getScore() {
        return score;
    }

    public int getMatchedPeakNum() {
        return matchedPeakNum;
    }

    public double getIonFrac() {
        return ionFrac;
    }

    public double getMatchedHighestIntensityFrac() {
        return matchedHighestIntensityFrac;
    }

    public double getExplainedAaFrac() {
        return explainedAaFrac;
    }

    public String getaScore() {
        return aScore;
    }


    public int compareTo(Peptide peptide) {
        if (this.getVarPtmContainingSeqNow().contentEquals(peptide.getVarPtmContainingSeqNow())) return 0;
        if (hasVarPTM() && peptide.hasVarPTM()) {
            if (score * this.getPriority() > peptide.score * peptide.getPriority()) {
                return 1;
            } else if (score * this.getPriority() < peptide.score * peptide.getPriority()) {
                return -1;
            } else {
                return this.subCompareTo(peptide);
            }
        } else {
            return this.subCompareTo(peptide);
        }
    }

    public int subCompareTo(Peptide peptide) {
        if (score > peptide.getScore()) {
            return 1;
        } else if (score < peptide.getScore()) {
            return -1;
        } else {
            if (matchedPeakNum > peptide.getMatchedPeakNum()) {
                return 1;
            } else if (matchedPeakNum < peptide.getMatchedPeakNum()) {
                return -1;
            } else {
                if (getPriority() > peptide.getPriority()) {
                    return 1;
                } else if (getPriority() < peptide.getPriority()) {
                    return -1;
                } else {
                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
                        return 1;
                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
                        return -1;
                    } else {
                        if (this.hashCode < peptide.hashCode) {
                            return -1;
                        } else if (this.hashCode > peptide.hashCode) {
                            return 1;
                        } else {
                            return 0;
                        }
                    }
                }
            }
        }
    }
}
