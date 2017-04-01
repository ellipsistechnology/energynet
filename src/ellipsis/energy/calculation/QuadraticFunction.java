package ellipsis.energy.calculation;

public class QuadraticFunction
{
    private double alpha, beta, gamma;
    
    public QuadraticFunction(double alpha, double beta, double gamma)
    {
        setAlpha(alpha);
        setBeta(beta);
        setGamma(gamma);
    }

    public double f(double x)
    {
        return alpha*x*x + beta*x + gamma;
    }
    
    public double df(double x)
    {
        return 2*alpha*x + beta;
    }

    public double getAlpha()
    {
        return alpha;
    }

    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
    }
    
    public double getBeta()
    {
        return beta;
    }

    public void setBeta(double beta)
    {
        this.beta = beta;
    }

    public double getGamma()
    {
        return gamma;
    }

    public void setGamma(double gamma)
    {
        this.gamma = gamma;
    }
}
