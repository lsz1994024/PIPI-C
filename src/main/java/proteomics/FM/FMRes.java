package proteomics.FM;

public class FMRes {
    public int sp;
    public int ep;
    public boolean settled = true;
    public int matchedPos = 0;
    public FMRes(int sp, int ep){
        this.sp = sp;
        this.ep = ep;
    }

}