package ellipsis.util;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public interface HessianFunction extends GradientFunction
{
	RealMatrix hessian(RealVector x);
}