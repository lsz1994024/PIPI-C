package proteomics.Types;

import java.util.HashSet;
import java.util.Set;

public class PeptideInfo implements Cloneable {
    public String freeSeq;
    public boolean isDecoy;
    public Set<String> protIdSet = new HashSet<>();
    public char leftFlank;
    public char rightFlank;


    public PeptideInfo(String freeSeq, boolean isDecoy, char leftFlank, char rightFlank) {
        this.isDecoy = isDecoy;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.freeSeq = freeSeq;
    }

    public PeptideInfo clone() throws CloneNotSupportedException {
        super.clone();
        PeptideInfo other = new PeptideInfo(freeSeq, isDecoy, leftFlank, rightFlank);
        other.protIdSet.addAll(this.protIdSet);
        return other;
    }

}
