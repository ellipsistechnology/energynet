package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;


/**
 * 
 * @author bmillar
 * @deprecated Use Bus
 */
@Deprecated
public class Junction extends Unit
{
	@Override
	public Unit copy(Grid grid) {
		throw new NotImplementedException();
	}

	@Override
	public Complex netPower() {
		return null;
	}
}
