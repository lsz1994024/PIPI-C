package proteomics.Types;


import java.util.Locale;
import java.util.TreeMap;

public class PosMassMap extends TreeMap<Integer, Double> {

    public PosMassMap() {
        super();
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(1000);
        for (Integer co : this.keySet()) {
            sb.append(String.format(Locale.US, "%.3f", this.get(co)));
            sb.append("@");
            sb.append(co.toString());
            sb.append(";");
        }
        return sb.toString();
    }

    public int hashCode() {
        return this.toString().hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof PosMassMap) {
            PosMassMap temp = (PosMassMap) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }

    public PosMassMap clone() {
        super.clone();
        PosMassMap other = new PosMassMap();
        other.clear();
        for (Integer co : this.keySet()) {
            other.put(co, this.get(co));
        }
        return other;
    }
}