package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;


/**
 * Represents various units within the power grid. All units in a {@link Grid} must
 * implement this interface.
 * @author bmillar
 *
 */
public abstract class Unit implements Comparable<Unit>
{
    private String name;
    
    public abstract Complex netPower();

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }
    
    public void wasRemovedFromBus(Bus parent)
    {
        // Do nothing by default.
    }
    
    public String toString() 
    { 
        return getName(); 
    }

	public abstract Unit copy(Grid grid);
	
	protected void copyInto(Unit copy)
	{
		copy.setName(getName());
	}
	
	@Override
	public int compareTo(Unit o) 
	{
        if(o == null)
            return Integer.MIN_VALUE;
		return this.name.compareTo(o.name);
	}
	
	@Override
	public boolean equals(Object u2) 
	{
		if(u2 == this)
			return true;
		if(u2 == null)
			return false;
		if(!(u2 instanceof Unit))
			return false;
		
		String name2 = ((Unit)u2).getName();
		if(name == null && name2 == null)
			return super.equals(u2);
		else if(name == null)
			return false;
		
		return name.equals(name2);
	}
	
	@Override
	public int hashCode() 
	{
		return name.hashCode();
	}
}
