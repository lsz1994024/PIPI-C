package proteomics.Types;


import java.util.ArrayList;
import java.util.List;

public class ExpTag implements Comparable<ExpTag> {
    public int isNorC = 0;
    private int hashCode;
    private double totalIntensity;
    private String freeAaString;
    private String ptmAaString;
    public List<ExpAa> expAaList;
    public List<Double> intensityList;

    public ExpTag(List<ExpAa> expAaList) {
        isNorC = expAaList.get(0).isNorC;
        intensityList = new ArrayList<>(expAaList.size() + 1);
        this.expAaList = expAaList;
        intensityList.add(expAaList.get(0).getHeadIntensity());
        for (int i = 0; i < expAaList.size(); i++) {
            intensityList.add(expAaList.get(i).getTailIntensity());
        }
        String toString = expAaList.get(0).toString();
        for (int i = 1; i < expAaList.size(); i++) {
            toString += "-" + expAaList.get(i).toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaList.size());
        StringBuilder sbPtmSeq = new StringBuilder(2 * expAaList.size() - 1);
        for (ExpAa aa : expAaList) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaList.get(0).getHeadIntensity();
        for (ExpAa aa : expAaList) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public ExpTag() {
    }

    public boolean equals(Object other) {
        return (other instanceof ExpTag) && (this.hashCode() == other.hashCode());
    }


    public int compareTo(ExpTag other) {
        return Double.compare(this.expAaList.get(0).getHeadLocation(), other.expAaList.get(0).getHeadLocation());
    }

    public String getFreeAaString() {
        return freeAaString;
    }

    public String getPtmAaString() {
        return ptmAaString;
    }

    public String toString() {
        return isNorC + "," + ptmAaString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return expAaList.get(0).getHeadLocation();
    }

    public double getTailLocation() {
        return expAaList.get(expAaList.size() - 1).getTailLocation();
    }

    public ExpTag clone() {
        ExpTag newTag = new ExpTag();
        newTag.isNorC = this.isNorC;
        newTag.hashCode = this.hashCode;
        newTag.totalIntensity = this.totalIntensity;
        newTag.freeAaString = this.freeAaString;
        newTag.ptmAaString = this.ptmAaString;
        newTag.expAaList = new ArrayList<>(this.expAaList);
        newTag.intensityList = new ArrayList<>(this.intensityList);
        return newTag;
    }

    public int size() {
        return expAaList.size();
    }


    public ExpTag revTag(double totalMass) {
        List<ExpAa> revExpAaList = new ArrayList<>(expAaList.size());
        for (int i = expAaList.size() - 1; i >= 0; i--) {
            revExpAaList.add(expAaList.get(i).revAA(totalMass));
        }
        ExpTag revedTag = new ExpTag(revExpAaList);
        revedTag.isNorC = isNorC;
        return revedTag;
    }

    public void appendTag(int pos, ExpTag newTag) {
        int diffLen = this.size() - pos;
        if (newTag.size() <= diffLen) return;

        int theLastIntensId = Math.min(this.size(), pos + newTag.size());
        for (int i = pos; i < theLastIntensId; ++i) {
            this.intensityList.set(pos, Math.max(this.intensityList.get(i), newTag.intensityList.get(i - pos)));
        }
        this.intensityList.set(theLastIntensId, Math.max(this.intensityList.get(theLastIntensId), newTag.intensityList.get(theLastIntensId - pos)));
        for (int i = diffLen; i < newTag.size(); i++) {
            this.intensityList.add(newTag.intensityList.get(i + 1));
            this.expAaList.add(newTag.expAaList.get(i));
        }

        String toString = expAaList.get(0).toString();
        for (int i = 1; i < expAaList.size(); i++) {
            toString += "-" + expAaList.get(i).toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaList.size());
        StringBuilder sbPtmSeq = new StringBuilder(2 * expAaList.size() - 1);
        for (ExpAa aa : expAaList) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        totalIntensity = 0;
        for (double intens : intensityList) {
            totalIntensity += intens;
        }
    }

    public ExpTag subTag(int start, int end) {
        List<ExpAa> subExpAaList = new ArrayList<>(end - start);
        for (int i = start; i < end; i++) {
            subExpAaList.add(expAaList.get(i));
        }
        ExpTag subTag = new ExpTag(subExpAaList);
        subTag.isNorC = isNorC;
        return subTag;
    }

}
