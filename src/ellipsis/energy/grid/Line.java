package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

/**
 * A line should be used to connect buses. From and To busses are parent and child units
 * respectively; only a single child can be added, subsequent additions of children will
 * overwrite the current child. As such, parents and children must only be Busses. A {@link ClassCastException}
 * will be thrown for other types of {@link Unit}s.
 * @author bmillar
 *
 */
public class Line extends Unit
{
    private double length;
    private double resistancePerMetre, inductancePerMetre;
    
    private Bus fromBus;
    private Bus toBus;

    public double getLength()
    {
        return length;
    }

    public void setLength(double length)
    {
        this.length = length;
    }

    public Bus getToBus()
    {
        return toBus;
    }

    public void setToBus(Bus toBus)
    {
        if(this.fromBus != null && this.fromBus.equals(toBus))
            throw new RuntimeException(getName()+": from and to busses must not be the same");
        
        if(this.toBus != null)
            this.toBus.removeChild(this);
        
        this.toBus = toBus;
        
        if(toBus != null)
            toBus.addChild(this);
    }

    public Bus getFromBus()
    {
        return fromBus;
    }

    public void setFromBus(Bus fromBus)
    {
        if(this.toBus != null && this.toBus.equals(fromBus))
            throw new RuntimeException(getName()+": from and to busses must not be the same");
     
        if(this.fromBus != null)
            this.fromBus.removeChild(this);
        
        this.fromBus = fromBus;
        
        if(fromBus != null)
            fromBus.addChild(this);
    }

    public double getResistancePerMetre()
    {
        return resistancePerMetre;
    }

    public void setResistancePerMetre(double resistancePerMetre)
    {
        this.resistancePerMetre = resistancePerMetre;
    }

    public double getInductancePerMetre()
    {
        return inductancePerMetre;
    }

	public void setImpedencePerMetre(Complex impedence601)
	{
		this.resistancePerMetre = impedence601.getReal();
		this.inductancePerMetre = impedence601.getImaginary();
	}

    public void setInductancePerMetre(double inductancePerMetre)
    {
        this.inductancePerMetre = inductancePerMetre;
    }
    
    public Complex impedance()
    {
        return impedancePerMetre().multiply(length);
    }
    
    public Complex impedancePerMetre()
    {
        return new Complex(resistancePerMetre, inductancePerMetre);
    }
    
    public Complex admittance()
    {
        return impedance().reciprocal();
    }

	@Override
	public Unit copy(Grid grid) 
	{
		Line copy = new Line();
		copyInto(copy);
		return copy;
	}

	/**
	 * Copies data into the given Line.
	 * This will set the to and from busses to be place-holders with the same name.
	 * @param copy
	 */
	protected void copyInto(Line copy) 
	{
		super.copyInto(copy);
		copy.setResistancePerMetre(getResistancePerMetre());
		copy.setInductancePerMetre(getInductancePerMetre());
		copy.setLength(getLength());
		
		Bus from = new Bus.BusPlaceHolder(null);
		from.setName(getFromBus().getName());
		copy.setFromBus(from);
		
		Bus to = new Bus.BusPlaceHolder(null);
		to.setName(getToBus().getName());
		copy.setToBus(to);
	}

	@Override
	public Complex netPower() 
	{
		return Complex.ZERO;
	}
}