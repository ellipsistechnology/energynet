package ellipsis.energy.util.constraint;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;

import ellipsis.energy.util.TimeDependent;
import ellipsis.util.TwiceDifferentiableFunction;

/**
 * Implements the bounds on a control.
 * Can specify device maximum, or minimum.
 * Only considers active power.
 * @author bmillar
 *
 */
public class DeviceConstraint implements TwiceDifferentiableFunction, TimeDependent, ConvexMultivariateRealFunction
{
	private int i;
	private double[] max; // [time]
	private double[] min; // [time]
	private int t = -1;
	private int dim;

	public DeviceConstraint(int i, double[] min, double[] max, int dimension)
	{
		this.i = i;
		this.max = max;
		this.min = min;
		if(min != null && max != null)
			throw new RuntimeException("Only min or max may be specified.");
		this.dim = dimension;
	}

	@Override
	public RealMatrix hessian(RealVector u)
	{
		return null; // null is the same as [0]
	}

	@Override
	public RealVector gradient(RealVector u)
	{
		RealVector grad = new ArrayRealVector(u.getDimension());
		if(max != null)
		{
			grad.setEntry(i, 1);
		}
		else
		{
			grad.setEntry(i, -1);
		}

		return grad;
	}

	@Override
	public double value(RealVector u)
	{
		if(max != null)
			return u.getEntry(i) - max[t];
		else
			return min[t] - u.getEntry(i);
	}

	@Override
	public void setTime(int t)
	{
		this.t = t;
	}
	
	@Override
	public String toString()
	{
		return min == null ? 
				"i="+i+", max["+t+"]="+max[t] : 
				"i="+i+", min["+t+"]="+min[t];
	}
	
	
	//// JOptimizer function implementation ////

	@Override
	public double value(double[] u)
	{
		if(max != null)
			return u[i] - max[t];
		else
			return min[t] - u[i];
	}

	@Override
	public double[] gradient(double[] u)
	{
		double[] grad = new double[u.length];
		if(max != null)
		{
			grad[i] = 1;
		}
		else
		{
			grad[i] =  -1;
		}

		return grad;
	}

	@Override
	public double[][] hessian(double[] u)
	{
		return new double[u.length][u.length];
	}

	@Override
	public int getDim()
	{
		return dim;
	}
	
	
	//// Accessors ////

	public int getI()
	{
		return i;
	}

	public double[] getMax()
	{
		return max;
	}

	public double[] getMin()
	{
		return min;
	}
}