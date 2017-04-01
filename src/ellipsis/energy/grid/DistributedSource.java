package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

public class DistributedSource extends Source
{
    private double maxAngle = Math.PI/2;
    private double minAngle = -Math.PI/2;
    
    public double getMaxAngle()
    {
        return maxAngle;
    }
    public void setMaxAngle(double maxAngle)
    {
        this.maxAngle = maxAngle;
    }
    public double getMinAngle()
    {
        return minAngle;
    }
    public void setMinAngle(double minAngle)
    {
        this.minAngle = minAngle;
    }
    
    public Complex getAvailableComplexPower()
    {
        double availablePower = getAvailablePower();
        double angle = getPowerOutput().getArgument();
        return new Complex(
                availablePower*Math.cos(angle),
                availablePower*Math.sin(angle));
    }
    
	@Override
	public Unit copy(Grid grid) 
	{
		DistributedSource copy = new DistributedSource();
		return copyInto(copy);
	}
	protected Unit copyInto(DistributedSource copy) 
	{
		super.copyInto(copy);
		copy.setMaxAngle(getMaxAngle());
		copy.setMinAngle(getMinAngle());
		return copy;
	}
	
	@Override
	public Complex netPower() 
	{
		return getPowerOutput();
	}
}
