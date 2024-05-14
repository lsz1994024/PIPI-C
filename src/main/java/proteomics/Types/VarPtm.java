package proteomics.Types;


import proteomics.PTM.InferPTM;

public class VarPtm implements Comparable<VarPtm> {

    public final double mass;
    public final char site;
    public final int priority;
    public final boolean onlyProteinTerminalIfnc;
    public String name = null;
    public byte position;
    public String classification = null;
    private final String toString;

    private final int hashCode;

    public VarPtm(double mass, char site, byte position, String name, String classification, int priority) {
        this.mass = mass;
        this.site = site;
        this.priority = priority;
        this.position = position;
        this.classification = classification;
        this.onlyProteinTerminalIfnc = false;
        this.name = name;
        toString = InferPTM.df3.format(mass) + "@" + site + "@" + position;
        hashCode = (InferPTM.df3.format(mass) + "@" + site).hashCode();
    }

    public String getStr() {
        return InferPTM.df3.format(mass) + "@" + site;
    }

    public int hashCode() {
        return hashCode;
    }

    @Override
    public String toString() {
        return toString;
    }

    public boolean equals(Object other) {
        if (other instanceof VarPtm) {
            VarPtm temp = (VarPtm) other;
            return temp.hashCode == hashCode;
        } else {
            return false;
        }
    }

    @Override
    public int compareTo(VarPtm other) {
        return 0;
    }
}
