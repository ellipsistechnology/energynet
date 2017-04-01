package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

public class ShuntAdmittance extends Unit implements Child
{
    private Complex admittance;
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

    public Complex getAdmittance()
    {
        return admittance;
    }

    public void setAdmittance(Complex admittance)
    {
        this.admittance = admittance;
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	ShuntAdmittance copy = new ShuntAdmittance();
    	copyInto(copy);
    	copy.setAdmittance(getAdmittance());
    	return copy;
    }

	@Override
	public Complex netPower() 
	{
		return Complex.ZERO;
	}
}
