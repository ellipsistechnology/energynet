package ellipsis.util;

import static ellipsis.util.VectorHelper.vector;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class BarrierSolver
{
	private static class Minimiser implements TwiceDifferentiableFunction
	{
		double z;
		HessianFunction barrier;
		HessianFunction function;

		@Override
		public RealVector gradient(RealVector x)
		{
			// Gradient of f(x):
//Timer.getGlobalTimer("min.f.g").start();
			RealVector gradF = function.gradient(x);
//Timer.getGlobalTimer("min.f.g").stop();

			// Gradient of barrier b(x):
//Timer.getGlobalTimer("min.b.g").start();
			RealVector gradB = barrier.gradient(x);
//Timer.getGlobalTimer("min.b.g").stop();

			// Gradient; grad(x) = z*grad f(x) + grad b(x):
			RealVector grad = gradF.mapMultiply(z).add(gradB);
			
			return grad;
		}

		@Override
		public RealMatrix hessian(RealVector x)
		{
			// 2nd Gradient of f(x):
//Timer.getGlobalTimer("min.f.h").start();
			RealMatrix HF = function.hessian(x);
//Timer.getGlobalTimer("min.f.h").stop();
			
			// 2nd Gradient of barrier b(x):
//Timer.getGlobalTimer("min.b.h").start();
			RealMatrix HB = barrier.hessian(x);
//Timer.getGlobalTimer("min.b.h").stop();
			
			// Gradient; H = z*Hf(x) + Hb(x):
			RealMatrix H = HF.scalarMultiply(z).add(HB);
			
			return H;
		}

		@Override
		public double value(RealVector x)
		{
			if(!(barrier instanceof TwiceDifferentiableFunction && function instanceof TwiceDifferentiableFunction))
				return 0;
			return z*((TwiceDifferentiableFunction)function).value(x) + ((TwiceDifferentiableFunction)barrier).value(x);
		}
	}
	
	// Parameters:
	private double initialZ = 0.01;
	private double epsilon = 1e-6;
	private double multiplierZ = 10;
	private int constraintDimension;
	private double stepSize = 1;

	private TwiceDifferentiableFunction barrier; // defining constraints
	private HessianFunction function; // to minimise
	
	public RealVector solve(RealVector x0)
	{
		RealVector x = x0;
		
		// Create the function for combining constraints with costs:
		Minimiser minimiser = new Minimiser();
		minimiser.barrier = barrier;
		minimiser.function = function;
		
		GradientDecent solver = new NewtonSolver();
		solver.setFunction(minimiser);
		solver.setStep(stepSize);
		solver.setEpsilon(epsilon);
		solver.setMaxIterations(40);

		// Iterate, reducing barrier strength:
		double stepSize = this.stepSize; // Need to change this locally some times.
		double finalZ = constraintDimension/epsilon;
		int k = 1;
//System.out.println("k,z,x");
		for(double z = initialZ; z < finalZ; )
		{
			if(z*multiplierZ >= finalZ) // Final round, make it more precise.
				solver.setMaxIterations(100);
			
			minimiser.z = z;
			solver.setStep(stepSize/k);
			solver.setEpsilon(z*epsilon*10);
			
			// Minimise_{x} z*f(x) + b(x):
			RealVector nextX = solver.solve(x);
			
			// Ensure the result is valid:
			if(Double.isNaN(barrier.value(nextX))) // We've stepped over the barrier.
			{
				stepSize /= 10.0; // Reduce step size and try again.
			}
			else
			{
				x = nextX;
				stepSize = this.stepSize;
			}
//System.out.println(k+","+z+","+x.getEntry(0));
			
			++k;
			z *= multiplierZ;
		}
		
		return x;
	}
	
	
	//// Accessors ////

	public double getInitialZ()
	{
		return initialZ;
	}

	public void setInitialZ(double initialZ)
	{
		this.initialZ = initialZ;
	}

	public double getEpsilon()
	{
		return epsilon;
	}

	public void setEpsilon(double epsilon)
	{
		this.epsilon = epsilon;
	}

	public double getMultiplierZ()
	{
		return multiplierZ;
	}

	public void setMultiplierZ(double multiplierZ)
	{
		this.multiplierZ = multiplierZ;
	}

	public int getConstraintDimension()
	{
		return constraintDimension;
	}

	public void setConstraintDimension(int constraintDimension)
	{
		this.constraintDimension = constraintDimension;
	}

	public double getStepSize()
	{
		return stepSize;
	}

	public void setStepSize(double stepSize)
	{
		this.stepSize = stepSize;
	}

	public TwiceDifferentiableFunction getBarrier()
	{
		return barrier;
	}

	/**
	 * The barrier {@link TwiceDifferentiableFunction#value(RealVector)} function must return NaN for invalid states.
	 * @param barrier
	 */
	public void setBarrier(TwiceDifferentiableFunction barrier)
	{
		this.barrier = barrier;
	}

	public HessianFunction getFunction()
	{
		return function;
	}

	public void setFunction(HessianFunction function)
	{
		this.function = function;
	}
	
	
	//// Debug ////
	
	public static String showMinimiser(Minimiser min)
	{
		StringBuffer sb = new StringBuffer();
		sb.append("x,f(x),b(x),min(x),grad f(x),grad b(x),grad min(x),grad^2 f(x),grad^2 b(x),grad^2 min(x)\n");
		for(double x = -2.5; x <= 1.5; x += 0.01)
		{
			RealVector _x = vector(x);
			
			sb.append(x);
			sb.append(',');
			sb.append(((TwiceDifferentiableFunction)min.function).value(_x));
			sb.append(',');
			sb.append(((TwiceDifferentiableFunction)min.barrier).value(_x));
			sb.append(',');
			sb.append(min.value(_x));
			sb.append(',');
			sb.append(min.function.gradient(_x).getEntry(0));
			sb.append(',');
			sb.append(min.barrier.gradient(_x).getEntry(0));
			sb.append(',');
			sb.append(min.gradient(_x).getEntry(0));
			sb.append(',');
			sb.append(min.function.hessian(_x).getEntry(0,0));
			sb.append(',');
			sb.append(min.barrier.hessian(_x).getEntry(0,0));
			sb.append(',');
			sb.append(min.hessian(_x).getEntry(0,0));
			sb.append('\n');
		}
		return sb.toString();
	}
}