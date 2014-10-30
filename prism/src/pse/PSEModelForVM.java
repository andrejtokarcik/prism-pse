package pse;

public class PSEModelForVM
{
    public PSEModelForVM
      ( int stCnt, int trCnt
      , double[] trRateLower
      , double[] trRateUpper
      , double[] trRatePopul
      , int[] trStSrc
      , int[] trStTrg
      , int[] trsI
      , int[] trsO
      , int[] trsIO
      , double[] trsNPVal
      , int[] trsNPTrg
      , int[] trsNPSrcBeg
      )
    {
        this.stCnt = stCnt;
        this.trCnt = trCnt;

        this.trRateLower = trRateLower;
        this.trRateUpper = trRateUpper;
        this.trRatePopul = trRatePopul;

        this.trStSrc = trStSrc;
        this.trStTrg = trStTrg;
        this.trsI = trsI;
        this.trsO = trsO;
        this.trsIO = trsIO;
        this.trsNPVal = trsNPVal;
        this.trsNPTrg = trsNPTrg;
        this.trsNPSrcBeg = trsNPSrcBeg;
    }

    public void vmMult(double min[], double resMin[], double max[], double resMax[])
    {
        System.arraycopy(min, 0, resMin, 0, min.length);
        System.arraycopy(max, 0, resMax, 0, max.length);

        for (int t : trsI) {
            final int v0 = trStSrc[t];
            final int v1 = trStTrg[t];

            resMin[v1] += trRateLower[t] * trRatePopul[t] * min[v0];
            resMax[v1] += trRateUpper[t] * trRatePopul[t] * max[v0];
        }

        for (int t : trsO) {
            final int v0 = trStSrc[t];

            resMin[v0] -= trRateUpper[t] * trRatePopul[t] * min[v0];
            resMax[v0] -= trRateLower[t] * trRatePopul[t] * max[v0];
        }

        for (int ii = 0; ii < trsIO.length; ) {
            final int t0 = trsIO[ii++]; // t0 goes from v0 to v1
            final int t1 = trsIO[ii++]; // t1 goes from v1 to v2

            final int v0 = trStSrc[t0];
            final int v1 = trStTrg[t0];
            // int v1 = st; // == trStTrg[t0] == trStSrc[t1]
            // int v2 = trStTrg[t1];

            // The rate params of t1 and t2 must be identical
            // assert trRateLower[t0] == trRateLower[t1];
            // assert trRateUpper[t0] == trRateUpper[t1];

            final double midSumNumeratorMin = trRatePopul[t0] * min[v0] - trRatePopul[t1] * min[v1];
            if (midSumNumeratorMin > 0.0) {
                resMin[v1] += trRateLower[t1] * midSumNumeratorMin;
            } else {
                resMin[v1] += trRateUpper[t1] * midSumNumeratorMin;
            }

            final double midSumNumeratorMax = trRatePopul[t0] * max[v0] - trRatePopul[t1] * max[v1];
            if (midSumNumeratorMax > 0.0) {
                resMax[v1] += trRateUpper[t1] * midSumNumeratorMax;
            } else {
                resMax[v1] += trRateLower[t1] * midSumNumeratorMax;
            }
        }

        for (int v0 = 0; v0 < stCnt; ++v0) {
            for (int ii = trsNPSrcBeg[v0]; ii < trsNPSrcBeg[v0 + 1]; ++ii) {
                final int v1 = trsNPTrg[ii];
                final double rate = trsNPVal[ii];

                resMin[v0] -= rate * min[v0];
                resMax[v0] -= rate * max[v0];

                resMin[v1] += rate * min[v0];
                resMax[v1] += rate * max[v0];
            }
        }
    }

    final private int stCnt;
    final private int trCnt;

    private double[] trRateLower;
    private double[] trRateUpper;
    final private double[] trRatePopul;

    final private int[] trStSrc;
    final private int[] trStTrg;

    final private int[] trsI;
    final private int[] trsO;
    final private int[] trsIO;
    final private double[] trsNPVal;
    final private int[] trsNPTrg;
    final private int[] trsNPSrcBeg;
}

