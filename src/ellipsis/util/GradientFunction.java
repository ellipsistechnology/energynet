package ellipsis.util;

import org.apache.commons.math3.linear.RealVector;

public interface GradientFunction
{
	RealVector gradient(RealVector x);
}