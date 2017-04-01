package ellipsis.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;

/**
 * A quadratic wrapper for {@link LinearLeastSquaresEstimator}.
 * @author bmillar
 *
 */
public class QuadraticEstimator extends LinearLeastSquaresEstimator implements TwiceDifferentiableFunction, SampleBasedEstimator
{
	protected List<Pair<RealVector, Double>> quadraticSamples = new ArrayList<>(); // for debug purposes only
	private RealMatrix defaultHessian;
	private RealVector defaultGradient;
	
	@Override
	public RealMatrix hessian(RealVector x)
	{
		if(dirty)
		{
			try
			{
				solve();
			}
			catch(SingularMatrixException sme)
			{
				return defaultHessian(x);
			}
		}
		else if(b == null)
		{
			return defaultHessian(x);
		}
		
		// If f(x) = a + c*x + d*x^2, 
		// then grad^2 f(x) = 2*diag(d):
		int dimension = b.getDimension()/2;
		RealVector d = b.getSubVector(0, dimension).mapMultiply(2.0);
		return MatrixUtils.createRealDiagonalMatrix(d.toArray());
	}

	private RealMatrix defaultHessian(RealVector x)
	{
		if(defaultHessian == null)
		{
			int dimension = x.getDimension();
			defaultHessian = new Array2DRowRealMatrix(dimension, dimension);
		}
		return defaultHessian;
	}

	@Override
	public RealVector gradient(RealVector x)
	{
		if(dirty)
		{
			try
			{
				solve();
			}
			catch(SingularMatrixException sme)
			{
				return defaultGradient(x);
			}
		}
		else if(b == null)
		{
			return defaultGradient(x);
		}
		
		// If f(x) = a + c*x + d*x^2, 
		// then grad f(x) = 2*diag(d)x + c:
		int dimension = b.getDimension()/2;
		RealVector d = b.getSubVector(0, dimension);
		RealVector c = b.getSubVector(dimension, dimension);
		return d.ebeMultiply(x).mapMultiply(2).add(c);
	}

	private RealVector defaultGradient(RealVector x)
	{
		if(defaultGradient == null)
		{
			defaultGradient = new ArrayRealVector(x.getDimension());
		}
		return defaultGradient;
	}

	@Override
	public double value(RealVector x)
	{
		return super.value(squared(x));
	}

	@Override
	public void addSample(RealVector x, double y)
	{
		super.addSample(squared(x), y);
		quadraticSamples.add(new Pair<RealVector, Double>(x, y));
		while(quadraticSamples.size() > maxSampleCount)
			quadraticSamples.remove(0);
	}
	
	@Override
	public List<Pair<RealVector, Double>> getSamplePoints()
	{
		return quadraticSamples;
	}
	
	/**
	 * Gives all parameters for estimating a general form quadratic function:
	 * b1x^2 + b2x, 
	 * where b1, b2 are vectors of the same dimension as x,
	 * and x^2 is the element by element square giving a vector with the same dimension as x.
	 * @param x
	 * @return
	 */
	private RealVector squared(RealVector x)
	{
		return x.ebeMultiply(x).append(x);
	}
}