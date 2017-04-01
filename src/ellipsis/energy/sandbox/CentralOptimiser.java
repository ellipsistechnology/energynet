package ellipsis.energy.sandbox;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import com.joptimizer.optimizers.OptimizationResponse;

public class CentralOptimiser
{
	private ConvexMultivariateRealFunction cost;
	private ConvexMultivariateRealFunction[] inequalities;
	
	public ArrayRealVector solve(double[] initialX)
	{
		OptimizationRequest or = new OptimizationRequest();
		or.setF0(cost);
		or.setFi(inequalities);
		or.setInitialPoint(initialX);
	
		JOptimizer opt = new JOptimizer();
		opt.setOptimizationRequest(or);
		int returnCode;
		try
		{
			returnCode = opt.optimize();
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}
		if(returnCode != OptimizationResponse.SUCCESS)
			throw new RuntimeException("Failed to optimize.");
		double[] u = opt.getOptimizationResponse().getSolution();
		return new ArrayRealVector(u);
	}
	
	public ConvexMultivariateRealFunction getCost()
	{
		return cost;
	}

	public void setCost(ConvexMultivariateRealFunction cost)
	{
		this.cost = cost;
	}

	public ConvexMultivariateRealFunction[] getInequalities()
	{
		return inequalities;
	}

	public void setInequalities(ConvexMultivariateRealFunction[] inequalities)
	{
		this.inequalities = inequalities;
	}
	
	
	//// Test ////
	
	/**
	 * f(x) = 0.5x^TAx + b^Tx + c
	 * f'(x) = Ax + b
	 * f''(x) = A
	 * g(x) = 
	 * @param args
	 */
	public static void main(String[] args)
	{
		double[] initialX = new double[] {2.0, -1.0};
		CentralOptimiser opt = new CentralOptimiser();
		RealMatrix A = new Array2DRowRealMatrix(new double[][]{{1,0},{0,1}});
		RealVector b = new ArrayRealVector(new double[]{1,1});
		opt.setCost(new ConvexMultivariateRealFunction()
		{
			@Override
			public double value(double[] X)
			{
				ArrayRealVector x = new ArrayRealVector(X);
				return x.dotProduct(A.operate(x)) + b.dotProduct(x);
			}
			
			@Override
			public double[][] hessian(double[] X)
			{
				return A.getData();
			}
			
			@Override
			public double[] gradient(double[] X)
			{
				ArrayRealVector x = new ArrayRealVector(X);
				return A.operate(x).add(b).toArray();
			}
			
			@Override
			public int getDim()
			{
				return 2;
			}
		});
		ArrayRealVector solution = opt.solve(initialX);
		for (int i = 0; i < solution.getDimension(); i++)
		{
			System.out.println(i+","+solution.getEntry(i));
		}
	}
}
