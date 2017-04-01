package ellipsis.energy.util;

import org.apache.commons.math3.linear.RealVector;

public interface StateDependent
{
	void setState(RealVector x);
}