package ellipsis.util;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class QuadraticFunction implements TwiceDifferentiableFunction
{
	private RealMatrix A;
	
	public RealMatrix getA()
	{
		return A;
	}

	public void setA(RealMatrix a)
	{
		A = a;
	}

	@Override
	public RealMatrix hessian(RealVector x)
	{
		return A.add(A.transpose()); // (A+A^T)
	}

	@Override
	public RealVector gradient(RealVector x)
	{
		return A.add(A.transpose()).operate(x); // (A+A^T)x
	}

	@Override
	public double value(RealVector x)
	{
		return A.operate(x).dotProduct(x); // x^TAx
	}
}