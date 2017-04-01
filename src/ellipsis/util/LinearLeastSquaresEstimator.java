package ellipsis.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;

import com.mls.util.Timer;

public class LinearLeastSquaresEstimator implements TwiceDifferentiableFunction, SampleBasedEstimator
{
	private double defaultValue;
	
	// Lazy solving:
	boolean dirty = false;
	
	// Samples:
	private RealMatrix X;
	private RealVector y;
	private List<RealVector> xs = new ArrayList<>();
	private List<Double> ys = new ArrayList<>();
	private List<Pair<RealVector, Double>> samples = new ArrayList<>(); // for debug purposes only
	protected int maxSampleCount = Integer.MAX_VALUE; // defaults to no limit
	
	// Cooefficents:
	protected RealVector b;
	protected double a;
	
	public void solve()
	{
		dirty = false;
		
		// Setup matrices:
		updateSamples();
		
		// (X^TX)b=X^Ty:
		RealMatrix X2 = X.transpose().multiply(X);
		RealVector Xy = X.transpose().operate(y);
		DecompositionSolver solver = new LUDecomposition(X2).getSolver();
		RealVector ba = solver.solve(Xy);
		int dimension = ba.getDimension()-1;
		b = ba.getSubVector(0, dimension);
		a = ba.getEntry(dimension);
	}

	// TODO this method is somewhat inefficient
	public void updateSamples()
	{
Timer.getGlobalTimer("updateSamples").start();
		if(xs.isEmpty())
			throw new RuntimeException("Cannot update without samples; use addSample(x, y)");
		
		int oldRowCount;
		RealMatrix newX;
		
		// Setup sample vector and matrix if needed:
		if(X == null)
		{
			oldRowCount = 0;
			int dimension = xs.get(0).getDimension();
			newX = new Array2DRowRealMatrix(xs.size(), dimension+1);
			y = new ArrayRealVector((Double[])ys.toArray(new Double[ys.size()]));
		}
		// Expand sample vector and matrix for new samples:
		else
		{
			oldRowCount = X.getRowDimension();
			newX = new Array2DRowRealMatrix(oldRowCount+xs.size(), X.getColumnDimension());
			newX.setSubMatrix(X.getData(), 0, 0);
			
			// Append new values onto y:
			RealVector newy = new ArrayRealVector(y.getDimension()+ys.size());
			newy.setSubVector(0, y);
			for(int i = 0; i < ys.size(); ++i)
			{
				newy.setEntry(oldRowCount+i, ys.get(i));
			}
			y = newy;
		}
		
		// Append new value onto X:
		for(int i = 0; i < xs.size(); ++i)
		{
			newX.setRowVector(oldRowCount+i, xs.get(i).append(1.0));
		}
		X = newX;
		
		// Trim if there is a sample limit:
		if(X.getRowDimension() > maxSampleCount)
		{
			int newStart = X.getRowDimension()-maxSampleCount;
			int newEnd = X.getRowDimension()-1;
			X = X.getSubMatrix(newStart, newEnd, 0, X.getColumnDimension()-1);
			y = y.getSubVector(newStart, newEnd);
		}
		
		// Clear samples:
		xs.clear();
		ys.clear();
Timer.getGlobalTimer("updateSamples").stop();
	}
	
	@Override
	public double value(RealVector x)
	{
		if(dirty)
		{
			try
			{
				solve();
			}
			catch(SingularMatrixException sme)
			{
//System.out.println(getClass().getName()+": Singular!");
				return defaultValue;
			}
		}
		else if(b == null)
		{
			return defaultValue;
		}
//System.out.println(getClass().getName()+": "+b+" dot "+x+" = "+b.dotProduct(x));
		return b.dotProduct(x)+a;
	}
	
	public void addSample(RealVector x, double y)
	{
		dirty = true;
		
		xs.add(x);
		ys.add(y);
		
		samples.add(new Pair<RealVector, Double>(x, y));
		while(samples.size() > maxSampleCount)
			samples.remove(0);
	}
	
	public RealMatrix getX()
	{
		return X;
	}
	
	public RealVector gety()
	{
		return y;
	}

	@Override
	public RealMatrix hessian(RealVector x)
	{
		int dimension = x.getDimension();
		return new Array2DRowRealMatrix(dimension, dimension);
	}

	@Override
	public RealVector gradient(RealVector x)
	{
		return b;
	}

	@Override
	public void setDefaultValue(double def)
	{
		this.defaultValue = def;
	}

	@Override
	public List<Pair<RealVector, Double>> getSamplePoints()
	{
		return samples;
	}

	public int getMaxSampleCount()
	{
		return maxSampleCount;
	}

	public void setMaxSampleCount(int maxSampleCount)
	{
		this.maxSampleCount = maxSampleCount;
	}
	
	public RealVector getB()
	{
		return b;
	}
	
	public double getA()
	{
		return a;
	}
}