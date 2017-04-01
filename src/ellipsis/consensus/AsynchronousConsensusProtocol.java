package ellipsis.consensus;

import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.average;
import static ellipsis.util.VectorHelper.vector;

import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/**
 * Implements the following asynchronous protocol/algorithm. 
 * For each local agent i:
 * 1. x_i(k+1) = x_i(k) + gamma*sum_neighbours(x_j(k) - x_i(k)) - w_i + b_i
 * 2. w_j := w_j + gamma*(x_j(k) - x_i(k)), for all j neighbours of i
 * 3. w_i := 0; 
 * 
 * @author bmillar
 *
 */
public class AsynchronousConsensusProtocol extends DiscreteConsensusProtocol
{
	/**
	 * f(x) = 0.5*x^TAx + b^Tx
	 * grad f(x) = Ax + b
	 *           = a^Tx+b
	 * Solution should be x = [-2.0 0.5 -3.0]^T with average -1.5.
	 */
	private static class Bias implements IndexedValuedFunction
	{
		double gamma = 0.8;
		double[] a = new double[]{0.5,  2.0, 0.3};
		double[] b = new double[]{1.0, -1.0, 0.9};
		RealVector x0;
		
		public Bias(RealVector x0)
		{
			this.x0 = new ArrayRealVector(x0);
		}
		
		@Override
		public double value(int i, RealVector x)
		{
			double x0_i = x0.getEntry(i);
			double step = -gamma*(a[i]*x0.getEntry(i) + b[i]);
			x0_i += step;
			x0.setEntry(i, x0_i);
			return step;
		}
		
		public RealVector getSolution()
		{
			return x0;
		}
	}

	static int dimension = 3;
	static Random rand = new Random(0);
	static RealVector x0 = vector(dimension, i -> rand.nextDouble());
	static Bias bias = new Bias(x0);
	public static void main(String[] args)
	{
		AsynchronousConsensusProtocol prot = new AsynchronousConsensusProtocol();
		prot.setIterations(100);
		double[][] w = new double[][]{
				{0,   1,   1},
				{1,   0,   1},
				{1,   1,   0}
		};
		prot.setWeights((i, j) -> w[i][j]);
		prot.setBias((n, v) -> 0.0);//FIXME bias);
		RealVector x_sol = prot.solve(x0);
		
		System.out.println("Average initial vector was "+average(x0));
		System.out.println("Average solution vector is "+average(x_sol));
		System.out.print("Initial vector is ");
		for (int i = 0; i < dimension; i++)
		{
			System.out.print(x0.getEntry(i)+",");
		}
		System.out.println();
		System.out.print("Solution is x^=");
		for (int i = 0; i < dimension; i++)
		{
			System.out.print(x_sol.getEntry(i)+",");
		}
		System.out.println();
		System.out.print("Solution is x*=");
		RealVector x = bias.getSolution();
		for (int i = 0; i < dimension; i++)
		{
			System.out.print(x.getEntry(i)+",");
		}
		System.out.println();
	}
	
	private RealVector w;
	
	@Override
	public RealVector solve(RealVector x0)
	{
		w = new ArrayRealVector(x0.getDimension());
		return super.solve(x0);
	}

	public RealVector step(RealVector x)
	{
		int dimension = x.getDimension();
		double stepSize = stepSize(dimension);
System.out.print(x.getEntry(0)+","+x.getEntry(1)+","+x.getEntry(2)+",");
System.out.print(w.getEntry(0)+","+w.getEntry(1)+","+w.getEntry(2)+",");
RealVector b = bias.getSolution();
System.out.print(b.getEntry(0)+","+b.getEntry(1)+","+b.getEntry(2));
		for(int i = 0; i < dimension; ++i)
		{
			if(rand.nextDouble() < 0.5) // 50/50 chance of update
				continue;

			double x_i = x.getEntry(i);
			
			// Prepare bias terms:
			double b_i = bias(i, x);
			
			// Update neighbours' w:
			double w_i = w.getEntry(i);
			for(int j = 0; j < dimension; ++j)
			{
				double w_j = w.getEntry(j);
				double x_j = x.getEntry(j);
				w_j += stepSize*weight(i, j)*(x_j - x_i);
				w.setEntry(j, w_j);
			}
			
			// Update x:
			final int _i = i;
			final double _x_i = x_i;
			double u_i = sum(j -> weight(_i,j)*(x.getEntry(j) - _x_i), dimension);
			x_i += stepSize*u_i - w_i + b_i;
			x.setEntry(i, x_i);
			
			// Reset local w:
			w.setEntry(i, 0.0);
		}
System.out.println();
		return x;
/*System.out.println(x.getEntry(0)+","+x.getEntry(1)+","+x.getEntry(2));
		RealVector nextX = new ArrayRealVector(x);
		RealVector _x = x;
		int dimension = x.getDimension();
		double stepSize = stepSize(dimension);
		for(int i = 0; i < dimension; ++i)
		{
			double x_i = x.getEntry(i);
			final int _i = i;
			double u_i = sum(j -> weight(_i,j)*(_x.getEntry(j) - x_i), dimension);
			double b_i = bias(i);
			nextX.setEntry(i, x_i + stepSize*u_i + b_i);
		}
		x = nextX;
		return x;*/
	}
}
