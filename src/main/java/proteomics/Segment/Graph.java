package proteomics.Segment;

import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Graph {

    private Set<Pair<Integer, Integer>> g = null;
    private final int length;
    private final Map<Integer, Boolean> visit = new HashMap<>();
    private final ArrayList<Integer> path = new ArrayList<>();
    private final ArrayList<ArrayList<Integer>> allPaths = new ArrayList<>();
    private int start;
    private int end;

    private Set<Integer> nodeSet = null;

    Graph(Set<Pair<Integer, Integer>> g, Set<Integer> nodeSet) {
        this.g = g;
        this.length = 4;
        for (int node : nodeSet) {
            visit.put(node, false);
        }
        this.nodeSet = nodeSet;
    }


    void dfs(int s, int e) {
        visit.put(s, true);
        path.add(s);
        if (s == e) {
            ArrayList<Integer> newPath = new ArrayList<>(path);
            allPaths.add(newPath);
        } else {
            for (int i : nodeSet) {
                if (!visit.get(i) && i != s && g.contains(new Pair<>(s, i))) {
                    dfs(i, e);
                }
            }
        }
        path.remove(path.size() - 1);
        visit.put(s, false);
    }

    public ArrayList<ArrayList<Integer>> getAllPaths(Set<Integer> startNodeSet, Set<Integer> endNodeSet) {
        for (int start : startNodeSet) {
            for (int end : endNodeSet) {
                dfs(start, end);
            }
        }
        return allPaths;
    }
}
