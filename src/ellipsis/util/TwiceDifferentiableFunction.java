package ellipsis.util;

import org.apache.commons.math3.linear.RealVector;

public interface TwiceDifferentiableFunction extends HessianFunction
{
	public double value(RealVector x);
}