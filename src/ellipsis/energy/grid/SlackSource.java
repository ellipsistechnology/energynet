package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

public class SlackSource extends Source
{
    private double resistance, inductance;
    private Complex voltage;

    public double getResistance()
    {
        return resistance;
    }

    public void setResistance(double resistance)
    {
        this.resistance = resistance;
    }

    public double getInductance()
    {
        return inductance;
    }

    public void setInductance(double inductance)
    {
        this.inductance = inductance;
    }

    public Complex getVoltage()
    {
        return voltage;
    }

    public void setVoltage(Complex voltage)
    {
        this.voltage = voltage;
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	SlackSource copy = new SlackSource();
    	copyInto(copy);
    	copy.setResistance(getResistance());
    	copy.setInductance(getInductance());
    	copy.setVoltage(getVoltage());
    	return copy;
    }

	@Override
	public Complex netPower() 
	{
		return Complex.ZERO;
	}
}
