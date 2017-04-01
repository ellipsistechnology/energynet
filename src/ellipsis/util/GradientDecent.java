package ellipsis.util;

import static ellipsis.util.ArrayHelper.maxMagnitude;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

public class GradientDecent
{
	public static boolean verbose = false; 
	
	// Parameters:
	private int maxIterations = 20;
	private double epsilon = 1e-6;
	protected GradientFunction function;
	protected double stepSize = 1.0;//0.1;
	
	// Accessible for debug/test purposes only:
	public int k;
	public RealVector grad;
	
	public RealVector solve(RealVector x0)
	{
		RealVector x = x0;
		for(k = 0; k < maxIterations; ++k)
		{
			// Gradient at x:
			RealVector step = step(x);
			if(step == null)
				break;
			if(step.isNaN())
				Util.breakpoint();
			
			// Update estimate of minimal x; x = x - grad(x):
			x = x.subtract(step);
			if(x.isInfinite())
				Util.breakpoint();
			
			// If gradient close enough to zero then finish:
			if(maxMagnitude(grad.toArray()) < epsilon)
				break;
System.out.println(x.getEntry(0));
		}
		
		return new ArrayRealVector(x);
	}

	public RealVector step(RealVector x)
	{
		grad = function.gradient(x);
		RealVector step = grad.mapMultiply(stepSize);
		return step;
	}
	
	
	//// Accessors ////

	public int getMaxIterations()
	{
		return maxIterations;
	}

	public void setMaxIterations(int maxIterations)
	{
		this.maxIterations = maxIterations;
	}

	public double getEpsilon()
	{
		return epsilon;
	}

	public void setEpsilon(double epsilon)
	{
		this.epsilon = epsilon;
	}

	public double getStep()
	{
		return stepSize;
	}

	public void setStep(double step)
	{
		this.stepSize = step;
	}

	public GradientFunction getFunction()
	{
		return function;
	}

	public void setFunction(GradientFunction function)
	{
		this.function = function;
	}
}
