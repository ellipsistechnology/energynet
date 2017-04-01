package ellipsis.energy.sandbox;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.HessianFunction;
import ellipsis.util.NewtonSolver;

public class Sandbox002
{
	public static void main(String[] args)
	{
		NewtonSolver ns = new NewtonSolver();
		ns.setFunction(new HessianFunction()
		{
			@Override
			public RealVector gradient(RealVector v)
			{
				// x^3 -6x^2+9x-4
				double x = v.getEntry(0);
				return new ArrayRealVector(new double[]{
						x*x*x - 6*x*x + 9*x - 4
				});
			}

			@Override
			public RealMatrix hessian(RealVector v)
			{
				//3x^2 -12x+9
				double x = v.getEntry(0);
				return new Array2DRowRealMatrix(new double[][]{
						{3*x*x - 12*x + 9}
				});
			}
		});
		
		RealVector x_star = ns.solve(new ArrayRealVector(new double[]{6}));
		System.out.println(6+": "+x_star.getEntry(0));
		
		System.out.println();
		x_star = ns.solve(new ArrayRealVector(new double[]{3.5}));
		System.out.println(3.5+": "+x_star.getEntry(0));
	}
}
