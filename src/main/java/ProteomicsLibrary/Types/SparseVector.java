package ProteomicsLibrary.Types;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SparseVector {

    private final Map<Integer, Double> sparseVector = new HashMap<>();

    public SparseVector() {
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public double get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }
}
