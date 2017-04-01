package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

public class Capacitor extends Load
{
    // TO DO add constraints: min, max, step size
    private double targetPowerFactor = 0.9;
    private double min = 0;
    private double max;
    
    public double getVAR()
    {
        return -getLoad().getImaginary();
    }
    public void setVAR(double var)
    {
        setLoad(new Complex(0, -var));
    }
    
    @Override
    public String toString()
    {
        return super.toString()+"("+getVAR()+")";
    }
    
    public double getTargetPowerFactor()
    {
        return targetPowerFactor;
    }
    
    public void setTargetPowerFactor(double targetPowerFactor)
    {
        this.targetPowerFactor = targetPowerFactor;
    }
    
	public double getMin() {
		return min;
	}
	public void setMin(double min) {
		this.min = min;
	}
	public double getMax() {
		return max;
	}
	public void setMax(double max) {
		this.max = max;
	}
	
	@Override
	public Unit copy(Grid grid) 
	{
		Capacitor copy = new Capacitor();
		copyInto(copy);
		copy.setMax(getMax());
		copy.setMin(getMin());
		copy.setTargetPowerFactor(targetPowerFactor);
		return copy;
	}
}
