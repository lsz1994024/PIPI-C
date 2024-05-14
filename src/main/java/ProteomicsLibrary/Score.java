package ProteomicsLibrary;


import java.util.*;

public class Score {

    public static int getMatchedIonNum(TreeMap<Double, Double> plMap, int localMaxMs2Charge, double[][] ionMatrix, double ms2Tolerance) {
        int K1 = 0;
        for (int i = 0; i < localMaxMs2Charge * 2; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : plMap.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++K1;
                        break;
                    }
                }
            }
        }
        return K1;
    }
}
