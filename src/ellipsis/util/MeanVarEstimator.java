package ellipsis.util;

public class MeanVarEstimator
{
    private long n = 0; // number of samples
    private double m; // mean
    private double ssd; // sum of squares of differences
    private double sv; // sample variance
    private double pv; // population variance
    
    public MeanVarEstimator()
    {
        // nothing to do
    }
    
    public void addSample(double x)
    {
        ++n;
        
        double m_new = m + (x-m)/n;
        if(n > 1)
        {
            ssd += (x-m)*(x-m_new);
            sv = ssd/(n-1);
            pv = ssd/n;
        }
        m = m_new;
    }
    
    public long getSampleCount()
    {
        return n;
    }
    
    public double getMean()
    {
        return m;
    }
    public double getSumOfSquaresOfDifferences()
    {
        return ssd;
    }
    public double getSampleVariance()
    {
        return sv;
    }
    public double getPopulationVariance()
    {
        return pv;
    }
}
