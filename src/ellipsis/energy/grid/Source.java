package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.QuadraticFunction;

public abstract class Source extends Unit implements Child
{   
    private static final double M = 1000000; // Mega
    
    private double Pmin, Pmax, alpha, beta, gamma;
    private Complex powerOutput = Complex.ZERO;
    private Bus parent;

    @Override
    public Bus getParent()
    {
        return parent;
    }

    @Override
    public void setParent(Bus parent)
    {
        this.parent = parent;
    }

    /**
     * Convenience setter for all cost data.
     * @param costData
     */
    public void setData(double Pmin, double Pmax, double alpha, double beta, double gamma)
    {
        setPmin(Pmin);
        setPmax(Pmax);
        setAlpha(alpha);
        setBeta(beta);
        setGamma(gamma);
    }
    
    public double costPerMW()
    {
        // Convert alpha and beta from Mega units
        return new QuadraticFunction(alpha/M/M, beta/M, gamma).f(getPowerOutput().getReal());
    }
    
    public double costPerMWh()
    {
        return new QuadraticFunction(alpha/M/M, beta/M, gamma).df(getPowerOutput().getReal()); 
    }
    
    /**
     * The current power available.
     * Returns {@link #getPmax()} by default.
     * @return
     */
    public double getAvailablePower()
    {
        return getPmax();
    }
    
    /**
     * Current power output.
     * @return
     */
    public Complex getPowerOutput()
    {
        return powerOutput;
    }

    public void setPowerOutput(Complex powerOutput)
    {
        if(powerOutput == null)
            throw new NullPointerException("Source power output cannot be null.");
        
        // Scale back to max power output if required:
        double s = powerOutput.abs();
        double max = getAvailablePower();
        if(s > max)
        {
            double angle = powerOutput.getArgument();
            powerOutput = new Complex(
                    max*Math.cos(angle),
                    max*Math.sin(angle));
        }
        
        this.powerOutput = powerOutput;
    }

    public void setReactivePowerOutput(double Q)
    {
        double P = getPowerOutput().getReal();
        double pmax = getAvailablePower();
        
        if(Q >= pmax)
        {
            Q = pmax;
            P = 0;
        } 
        else
        {
            double Psquared = P*P;
            double Qsquared = Q*Q;
            if(Math.sqrt(Psquared+Qsquared) > pmax)
            {
                P = Math.sqrt(pmax*pmax-Qsquared);
            }
        }
        
        setPowerOutput(P, Q);
    }

    public void setPowerOutput(double watts, double vars)
    {
        setPowerOutput(new Complex(watts, vars));
    }

    /**
     * The device minimum power.
     * @return
     */
    public double getPmin()
    {
        return Pmin;
    }

    public void setPmin(double pmin)
    {
        Pmin = pmin;
    }

    /**
     * The device maximum power.
     * @return
     */
    public double getPmax()
    {
        return Pmax;
    }

    public void setPmax(double pmax)
    {
        Pmax = pmax;
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

	public void copyInto(Source copy) 
	{
		super.copyInto(copy);
		copy.setAlpha(getAlpha());
		copy.setBeta(getBeta());
		copy.setGamma(getGamma());
		copy.setPmax(getPmax());
		copy.setPmin(getPmin());
		copy.setPowerOutput(getPowerOutput());
	}
}
