package ProteomicsLibrary;

import ProteomicsLibrary.Types.SparseVector;

import java.util.*;

public class SpecProcessor {
    private static final double defaultIntensity = 1;
    private static final int xcorrOffset = 75;
    private static final double removePrecursorPeakTolerance = 1.5;
    private final MassTool massTool;

    public SpecProcessor(MassTool massTool) {
        this.massTool = massTool;
    }

    public TreeMap<Double, Double> preSpectrumTopNStyleWithChargeLimit(Map<Double, Double> inputPL, double precursorMass, int precursorCharge, int topN) {

        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge);
        List<Map.Entry<Double, Double>> plList = new ArrayList<>(outputPL.entrySet());
        Collections.sort(plList, Map.Entry.comparingByValue(Comparator.reverseOrder()));
        plList = plList.subList(0, Math.min(plList.size(), 600));
        Collections.sort(plList, Map.Entry.comparingByKey());

        TreeMap<Double, Double> newPlMap = new TreeMap<>();

        for (Map.Entry<Double, Double> peak : plList) {
            newPlMap.put(peak.getKey(), peak.getValue());
        }
        if (newPlMap.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return topNStyleNormalization(sqrtPL(new TreeMap<>(newPlMap.subMap(0d, precursorMass))), topN);
        }
    }

    public SparseVector prepareXCorr(TreeMap<Double, Double> plMap, boolean flankingPeaks) {
        if (plMap.isEmpty()) {
            return new SparseVector();
        } else {
            double[] plArray = new double[massTool.mzToBin(plMap.lastKey()) + 1];
            for (double mz : plMap.keySet()) {
                if (Math.abs(plMap.get(mz)) > 1e-6) {
                    int idx = massTool.mzToBin(mz);
                    plArray[idx] = Math.max(plMap.get(mz), plArray[idx]);
                }
            }
            return prepareXCorr(plArray, flankingPeaks);
        }
    }

    public SparseVector prepareXCorr(double[] plArray, boolean flankingPeaks) {
        SparseVector xcorrPL = new SparseVector();
        int offsetRange = 2 * xcorrOffset + 1;
        double factor = 1 / (double) (offsetRange - 1);
        double mySum = 0;
        for (int i = 0; i < xcorrOffset; ++i) {
            mySum += plArray[i];
        }

        double[] tempArray = new double[plArray.length];
        for (int i = xcorrOffset; i < plArray.length + xcorrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            tempArray[i - xcorrOffset] = (mySum - plArray[i - xcorrOffset]) * factor;
        }

        for (int i = 1; i < plArray.length; ++i) {
            double temp = plArray[i] - tempArray[i];
            if (flankingPeaks) {
                temp += (plArray[i - 1] - tempArray[i - 1]) * 0.5;
                if (i + 1 < plArray.length) {
                    temp += (plArray[i + 1] - tempArray[i + 1]) * 0.5;
                }
            }
            if (Math.abs(temp) > 1e-6) {
                xcorrPL.put(i, temp);
            }
        }

        return xcorrPL;
    }

    public SparseVector digitizePL(TreeMap<Double, Double> plMap) {
        SparseVector digitizedPL = new SparseVector();
        for (double mz : plMap.keySet()) {
            int idx = massTool.mzToBin(mz);
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                digitizedPL.put(idx, Math.max(plMap.get(mz), digitizedPL.get(idx)));
            }
        }
        return digitizedPL;
    }

    public static TreeMap<Double, Double> topNStyleNormalization(TreeMap<Double, Double> inputPL, int localTopN) {
        if (inputPL.isEmpty()) {
            return new TreeMap<>();
        } else {
            TreeMap<Double, Double> outputPL = new TreeMap<>();
            double minMz = inputPL.firstKey();
            double maxMz = inputPL.lastKey();
            double leftMz = minMz;
            while (leftMz < maxMz) {
                double rightMz = Math.min(leftMz + 100, maxMz);
                NavigableMap<Double, Double> subMap;
                if (rightMz < maxMz) {
                    subMap = inputPL.subMap(leftMz, true, rightMz, false);
                } else {
                    subMap = inputPL.subMap(leftMz, true, rightMz, true);
                }
                if (!subMap.isEmpty()) {
                    Double[] intensityArray = subMap.values().toArray(new Double[0]);
                    Arrays.sort(intensityArray, Comparator.reverseOrder());
                    double temp1 = defaultIntensity / intensityArray[0];
                    double temp2 = subMap.size() > localTopN ? intensityArray[localTopN] : 0;
                    for (double mz : subMap.keySet()) {
                        if (subMap.get(mz) > temp2) {
                            outputPL.put(mz, subMap.get(mz) * temp1);
                        }
                    }
                }
                leftMz = rightMz;
            }

            return outputPL;
        }
    }

    private static TreeMap<Double, Double> removeCertainPeaks(Map<Double, Double> peakMap, double precursorMass, int precursorCharge) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();
        double precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
        for (double mz : peakMap.keySet()) {
            if (((mz < 112.5) || (mz > 121.5)) && (mz > 50)) {
                if ((peakMap.get(mz) > 1e-6) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put(mz, peakMap.get(mz));
                }
            }
        }

        return mzIntensityMap;
    }

    private static TreeMap<Double, Double> sqrtPL(TreeMap<Double, Double> plMap) {
        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            if (plMap.get(mz) > 1e-6) {
                sqrtPlMap.put(mz, Math.sqrt(plMap.get(mz)));
            }
        }
        return sqrtPlMap;
    }
}
