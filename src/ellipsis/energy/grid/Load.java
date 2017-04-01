package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;

import ellipsis.energy.smartgrid.global_local.Forecast;

public class Load extends Unit implements Child
{
    private Complex load = Complex.ZERO;
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

    public Complex getLoad()
    {
        return load;
    }

    public void setLoad(Complex load)
    {
        if(load == null)
            throw new NullPointerException("Load cannot be null");
        this.load = load;
    }
    
    public Complex setLoad(double watts, double var)
    {
        this.load = new Complex(watts, var);
        return this.load;
    }

    public Complex getExpectedMinimumPower(Forecast forecast)
    {
        double angle = load.getArgument();
        double loadMean = forecast.getLoadMean();
        double loadSD = forecast.getLoadStandardDeviation();
        double power = loadMean - 3*loadSD;
        if(power < 0)
            power = 0;
        return ComplexUtils.polar2Complex(power, angle);
    }

    public Complex getExpectedMaximumPower(Forecast forecast)
    {
        double angle = load.getArgument();
        double loadMean = forecast.getLoadMean();
        double loadSD = forecast.getLoadStandardDeviation();
        double power = loadMean + 3*loadSD;
        return ComplexUtils.polar2Complex(power, angle);
    }

	@Override
	public Unit copy(Grid grid) 
	{
		Load copy = new Load();
		copyInto(copy);
		return copy;
	}

	protected void copyInto(Load copy) 
	{
		super.copyInto(copy);
		copy.setLoad(getLoad());
	}

	@Override
	public Complex netPower()
	{
		return load.negate();
	}
}
