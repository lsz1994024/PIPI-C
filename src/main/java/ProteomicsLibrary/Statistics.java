package ProteomicsLibrary;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

public class Statistics {

    private static final NormalDistribution normalDistribution = new NormalDistribution(0, 1);

    public static double calMean(double[] input) throws Exception {
        if (input.length > 0) {
            double mean = 0;
            for (double v : input) {
                mean += v;
            }
            return mean / input.length;
        } else {
            throw new Exception("There is no element in the input array.");
        }
    }

    public static double calPearsonCorrelationCoefficient(double[] input1, double[] input2) throws Exception {
        if (input1.length != input2.length) {
            throw new Exception("Two vectors' lengths are different.");
        }

        double mean1 = calMean(input1);
        double mean2 = calMean(input2);
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        for (int i = 0; i < input1.length; ++i) {
            double c1 = input1[i] - mean1;
            double c2 = input2[i] - mean2;
            temp1 += c1 * c2;
            temp2 += Math.pow(c1, 2);
            temp3 += Math.pow(c2, 2);
        }
        return (temp1 == 0 || temp2 == 0) ? 0 : temp1 / (Math.sqrt(temp2 * temp3));
    }


}
