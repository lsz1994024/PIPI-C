package ProteomicsLibrary.Types;

import java.util.HashSet;
import java.util.Set;

public class SparseBooleanVector {

    private final Set<Integer> sparseVector;

    public SparseBooleanVector(Set<Integer> sparseVector) {
        this.sparseVector = new HashSet<>(sparseVector);
    }

    public int dot(SparseBooleanVector other) {
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector);
        intersectedKeys.retainAll(other.sparseVector);
        return intersectedKeys.size();
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(sparseVector.size() * 6);
        for (int idx : sparseVector) {
            sb.append(idx);
            sb.append(";");
        }
        return sb.toString();
    }
}
