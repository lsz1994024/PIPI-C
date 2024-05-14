package proteomics.FM;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;

public class WaveletTreeNode {
    private PreprocessRank preprocessRank;
    private WaveletTreeNode leftNode;
    private WaveletTreeNode rightNode;
    private HashMap<Character, Boolean> characterBitMap;

    public WaveletTreeNode(ArrayList<Character> S, ArrayList<Character> sortedAlphabet){
        int n = sortedAlphabet.size();
        if (n < 2) {
            return;
        }
        characterBitMap = new HashMap<Character, Boolean>();
        ArrayList<Character> sortedLeftAlphabet = new ArrayList<Character>();
        ArrayList<Character> sortedRightAlphabet = new ArrayList<Character>();
        ArrayList<Character> leftS = new ArrayList<Character>();
        ArrayList<Character> rightS = new ArrayList<Character>();

        for (int i = 0; i < (n + 1) / 2; i++) {
            sortedLeftAlphabet.add(sortedAlphabet.get(i));
            characterBitMap.put(sortedAlphabet.get(i), false);
        }
        for (int i = (n + 1) / 2; i < n; i++) {
            sortedRightAlphabet.add(sortedAlphabet.get(i));
            characterBitMap.put(sortedAlphabet.get(i), true);
        }


        try {
            Field tableField = HashMap.class.getDeclaredField("table");
            tableField.setAccessible(true);
            Object[] table;
            table = (Object[]) tableField.get(characterBitMap);

            FMIndex.space += 32 * characterBitMap.size() +  4 * (table == null ? 0 : table.length);
        } catch (Exception e) {
            e.printStackTrace();
        }

        boolean[] B = new boolean[S.size()];
        for (int i = 0; i < S.size(); i++) {
            B[i] = characterBitMap.get(S.get(i));
            if (B[i]) {
                rightS.add(S.get(i));
            } else {
                leftS.add(S.get(i));
            }
        }

        preprocessRank = new PreprocessRank(B);

        leftNode = new WaveletTreeNode(leftS, sortedLeftAlphabet);
        rightNode = new WaveletTreeNode(rightS, sortedRightAlphabet);
    }

    public int Occ(int i, char c){
        if (leftNode == null) {
            return i;
        }
        if (characterBitMap.get(c)) {
            return rightNode.Occ(preprocessRank.rank1(i - 1), c);
        } else {
            return leftNode.Occ(i - preprocessRank.rank1(i - 1), c);
        }
    }
}
