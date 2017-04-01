package ellipsis.energy.util.constraint;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;

import ellipsis.energy.util.TimeDependent;
import ellipsis.util.TwiceDifferentiableFunction;

public class ApparentPowerConstraint implements TwiceDifferentiableFunction, TimeDependent, ConvexMultivariateRealFunction
{
	// Read-only:
	private int i;
	private RealMatrix H;
	private double[][] H_array;
	private double[] sSquared; // [time]
	private int dim;
	
	// Variables:
	private int t;

	public ApparentPowerConstraint(int i, double[][] max, int dimension)
	{
		sSquared = new double[max.length];
		for(int t = 0; t < max.length; ++t)
		{
			sSquared[t] = max[t][i]*max[t][i];
			if(sSquared[t] == 0)
				sSquared[t] = 1e-9; // can't let this be zero
		}
		
		init(i, dimension);
	}

	public ApparentPowerConstraint(int i, double[] sSquared, int dimension)
	{
		this.sSquared = sSquared;
		init(i, dimension);
	}

	private void init(int i, int dimension)
	{
		this.i = i;
		this.dim = dimension;

		H_array = new double[dimension][dimension];
		H_array[i][i] = 2;
		int qi = i+dimension/2;
		H_array[qi][qi] = 2;
	}

	/**
	 * d^2c(u)/d^2P^2 = 2
	 * d^2c(u)/d^2Q^2 = 2
	 */
	@Override
	public RealMatrix hessian(RealVector u)
	{
		if(H == null)
		{
			int dimension = u.getDimension();
			H = new Array2DRowRealMatrix(dimension, dimension);
			H.setEntry(i, i, 2);
			int qi = i+u.getDimension()/2;
			H.setEntry(qi, qi, 2);
		}
		return H;
	}

	/**
	 * dc(u)/dP = 2P
	 * dc(u)/dQ = 2Q
	 */
	@Override
	public RealVector gradient(RealVector u)
	{
		double p = u.getEntry(i);
		double q = u.getEntry(i+u.getDimension()/2);
		RealVector grad = new ArrayRealVector(u.getDimension());
		grad.setEntry(i, 2*p);
		grad.setEntry(i+u.getDimension()/2, 2*q);

		return grad;
	}

	/**
	 * c(u) = P^2 + Q^2 - S^2 <= 0
	 * where S^2 is the maximum apparent power.
	 */
	@Override
	public double value(RealVector u)
	{
		double p = u.getEntry(i);
		double q = u.getEntry(i+u.getDimension()/2);
		double v = p*p + q*q - sSquared[t];

		return v;
	}
	
	@Override
	public String toString()
	{
		return "i="+i+", p^2+q^2<s_"+t+"^2="+sSquared[t];
	}
	
	
	//// JOptimizer implementation ////

	@Override
	public double value(double[] u)
	{
		double p = u[i];
		double q = u[i+dim/2];
		double v = p*p + q*q - sSquared[t];

		return v;
	}

	@Override
	public double[] gradient(double[] u)
	{
		double p = u[i];
		double q = u[i+dim/2];
		double[] grad = new double[dim];
		grad[i] = 2*p;
		grad[i+dim/2] = 2*q;

		return grad;
	}

	@Override
	public double[][] hessian(double[] u)
	{
		return H_array;
	}

	@Override
	public int getDim()
	{
		return dim;
	}
	
	
	//// Accessors ////
	
	@Override
	public void setTime(int t)
	{
		this.t = t;
	}

	public double[] getSSquared()
	{
		return sSquared;
	}
}