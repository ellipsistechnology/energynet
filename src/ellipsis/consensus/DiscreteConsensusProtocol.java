package ellipsis.consensus;

import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.average;
import static ellipsis.util.VectorHelper.vector;

import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class DiscreteConsensusProtocol
{
	public static interface NeighbourFunction
	{
		double value(int i, int j);
	}
	
	public static interface IndexedValuedFunction
	{
		double value(int i, RealVector x);
	}
	
	public static void main(String[] args)
	{
		DiscreteConsensusProtocol prot = new DiscreteConsensusProtocol();
		int dimension = 3;
		prot.setBias((i, x) -> 0.0);
		double[][] w = new double[][]{
				{0,   1,   1},
				{1,   0,   1},
				{1,   1,   0}
		};
		prot.setWeights((i, j) -> w[i][j]);
		Random rand = new Random(1);
		RealVector x0 = vector(dimension, i -> rand.nextDouble());
		RealVector x_sol = prot.solve(x0);
		
		System.out.println("Average initial vector was "+average(x0));
		System.out.print("Initial vector is ");
		for (int i = 0; i < dimension; i++)
		{
			System.out.print(x0.getEntry(i)+",");
		}
		System.out.println();
		System.out.print("Solution is ");
		for (int i = 0; i < dimension; i++)
		{
			System.out.print(x_sol.getEntry(i)+",");
		}
		System.out.println();
	}
	
	private int iterations = 30;
	private NeighbourFunction weights;
	private IndexedValuedFunction bias;

	public RealVector solve(RealVector x0) // TODO Change to return list of vectors (a vector per agent) or a matrix.
	{
		RealVector x = new ArrayRealVector(x0);
		
		for(int k = 0; k < iterations; ++k)
		{
			x = step(x);
		}
		
		return x;
	}

	public double stepSize(int dimension)
	{
		double degree = 0;
		for(int i = 0; i < dimension; ++i)
		{
			double sum_i = 0;
			for(int j = 0; j < dimension; ++j)
			{
				if(i != j)
					sum_i += weights.value(i, j);
			}
			degree = Math.max(degree, sum_i);
		}
		double stepSize = 0.9/degree;
		return stepSize;
	}

	public RealVector step(RealVector x)
	{
System.out.println(x.getEntry(0)+","+x.getEntry(1)+","+x.getEntry(2));
		RealVector nextX = new ArrayRealVector(x);
		RealVector _x = x;
		int dimension = x.getDimension();
		double stepSize = stepSize(dimension);
		for(int i = 0; i < dimension; ++i)
		{
			double x_i = x.getEntry(i);
			final int _i = i;
			double u_i = sum(j -> weights.value(_i,j)*(_x.getEntry(j) - x_i), dimension);
			double b_i = bias.value(i, x);
			nextX.setEntry(i, x_i + stepSize*u_i + b_i);
		}
		x = nextX;
		return x;
	}
	
	
	//// Accessors ////

	public int getIterations()
	{
		return iterations;
	}

	public void setIterations(int iterations)
	{
		this.iterations = iterations;
	}
	
	public double weight(int i, int j)
	{
		return weights.value(i, j);
	}

	public NeighbourFunction getWeights()
	{
		return weights;
	}

	public void setWeights(NeighbourFunction weights)
	{
		this.weights = weights;
	}

	public IndexedValuedFunction getBias()
	{
		return bias;
	}

	public void setBias(IndexedValuedFunction bias)
	{
		this.bias = bias;
	}
	
	public double bias(int i, RealVector x)
	{
		return bias.value(i, x);
	}
}