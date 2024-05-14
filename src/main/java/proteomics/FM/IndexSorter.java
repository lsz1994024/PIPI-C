package proteomics.FM;

import java.util.Arrays;
import java.util.Comparator;

public class IndexSorter implements Comparator<Integer> {

    private final char[] values;
    private final Integer[] indexes;

    public IndexSorter(char[] text, Integer[] indexes) {
        this.values = text;
        this.indexes = indexes;
    }

    public void sort() {
        Arrays.sort(indexes, this);
    }

    @Override
    public int compare(Integer arg0, Integer arg1) {
        return Character.compare(values[arg0], values[arg1]);
    }

}

