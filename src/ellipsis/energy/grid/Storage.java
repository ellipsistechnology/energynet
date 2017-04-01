package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

@Deprecated
public class Storage extends Unit
{
	@Override
	public Unit copy(Grid grid) 
	{
		return new Storage();
	}

	@Override
	public Complex netPower() 
	{
		return Complex.ZERO;
	}
}
