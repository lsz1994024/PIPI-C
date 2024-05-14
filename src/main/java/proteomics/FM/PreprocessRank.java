package proteomics.FM;

public class PreprocessRank {
    public boolean[] B;

    private int n;
    private int[] boundaryRank;
    private int[][] smallRank;
    private int t;

    public PreprocessRank(boolean[] B) {
        this.B = B;

        n = B.length;
        FMIndex.space += n;

        t = (int) Math.ceil((Math.log(n) / Math.log(2)) / 2);

        int numberOfBlocks = (int) Math.ceil((double) n / t);

        boundaryRank = new int[numberOfBlocks];

        FMIndex.space += numberOfBlocks * 4;

        boundaryRank[0] = 0;

        int numberOfOnes = 0;
        int boundaryRankIndex = 1;
        for (int i = 0; i < (numberOfBlocks - 1) * t; i++) {
            if (B[i]) {
                numberOfOnes++;
            }
            if ((i + 1) % t == 0) {
                boundaryRank[boundaryRankIndex] = numberOfOnes;
                boundaryRankIndex++;
            }
        }

        int smallRankRows = (int) Math.pow(2, t);
        smallRank = new int[smallRankRows][t];

        FMIndex.space += smallRankRows * t * 4;


        for (int i = 0; i < smallRankRows; i++) {
            numberOfOnes = 0;
            for (int j = 0; j < t; j++) {
                numberOfOnes += (i >> (t - j - 1)) & 1;
                smallRank[i][j] = numberOfOnes;
            }
        }
    }

    public int rank1(int i) {
        if (i < 0) {
            return 0;
        }
        int j = Math.floorDiv(i, t);
        int y = 0;
        for (int p = 0; p < t && (j * t + p) < n; p++) {
            if (B[j * t + p]) {
                y |= 1 << (t - p - 1);
            }
        }
        int k = Math.floorMod(i, t);

        return boundaryRank[j] + smallRank[y][k];
    }
}
