package proteomics.Types;

public class ExpAa implements Comparable<ExpAa> {
    private final String aa;
    private final char ptmFreeAA;
    private final double headLocation;
    private final double tailLocation;
    private final double headIntensity;
    private final double tailIntensity;
    public int isNorC;
    private final int hashCode;

    public ExpAa(String aa, char ptmFreeAA, double headLocation, double tailLocation, double headIntensity, double tailIntensity, int isNorC) {
        this.aa = aa;
        this.ptmFreeAA = ptmFreeAA;
        this.headLocation = headLocation;
        this.tailLocation = tailLocation;
        this.headIntensity = headIntensity;
        this.tailIntensity = tailIntensity;
        this.isNorC = isNorC;
        String toString = headLocation + "." + aa + "." + isNorC + "." + tailLocation;
        hashCode = toString.hashCode();
    }

    public ExpAa revAA(double totalMass) {
        return new ExpAa(aa, ptmFreeAA, totalMass - tailLocation, totalMass - headLocation, tailIntensity, headIntensity, isNorC);
    }


    public String getAA() {
        return aa;
    }

    public char getPtmFreeAA() {
        return ptmFreeAA;
    }


    public int compareTo(ExpAa other) {
        if (headLocation > other.headLocation) {
            return 1;
        } else if (headLocation < other.headLocation) {
            return -1;
        } else {
            return Double.compare(tailIntensity, other.tailIntensity);
        }
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ExpAa) {
            ExpAa temp = (ExpAa) other;
            return this.hashCode() == temp.hashCode();
        } else {
            return false;
        }
    }

    public ExpAa clone() throws CloneNotSupportedException {
        super.clone();
        return new ExpAa(aa, ptmFreeAA, headLocation, tailLocation, headIntensity, tailIntensity, isNorC);
    }

    public double getHeadLocation() {
        return headLocation;
    }

    public double getTailLocation() {
        return tailLocation;
    }

    public double getHeadIntensity() {
        return headIntensity;
    }

    public double getTailIntensity() {
        return tailIntensity;
    }

    public double getTotalIntensity() {
        return tailIntensity + headIntensity;
    }
}