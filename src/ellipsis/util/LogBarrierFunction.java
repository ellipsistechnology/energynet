package ellipsis.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * Given a set of inequality constraints in the form c(x) <= 0 a barrier
 * function based on the natural logarithm is given: -sum{ln(-c(x))}.
 * @author bmillar
 *
 */
public class LogBarrierFunction implements BarrierFunction
{
	private List<TwiceDifferentiableFunction> constraints = new ArrayList<>();

//public static int count_g_constraint_v = 0;
	/**
	 * sum{(-1/c(x))*(grad(c(x)))}
	 */
	@Override
	public RealVector gradient(RealVector x)
	{
		RealVector sum = new ArrayRealVector(x.getDimension());
		for (TwiceDifferentiableFunction constraint : constraints)
		{
//Timer.getGlobalTimer("LBF.g.constraint.v_").start();
//++count_g_constraint_v;
			double c = constraint.value(x);
//Timer.getGlobalTimer("LBF.g.constraint.v_").stop();
//Timer.getGlobalTimer("LBF.g.constraint.g").start();
			RealVector gradc = constraint.gradient(x);
//Timer.getGlobalTimer("LBF.g.constraint.g").stop();
			RealVector summand = gradc.mapMultiply(-1/c);
			sum = sum.add(summand);
		}
		
		return sum;
	}

//public static int count_h_constraint_v = 0;
	/**
	 * sum{ 1/c^2(x))grad(c(x))(grad(c(x))^T + (-1/c(x))*(H(c(x))) }
	 */
	@Override
	public RealMatrix hessian(RealVector x)
	{
		int dimension = x.getDimension();
		RealMatrix sum = new Array2DRowRealMatrix(dimension, dimension);
		for (TwiceDifferentiableFunction constraint : constraints)
		{
//Timer.getGlobalTimer("LBF.h.constraint.v_").start();
//++count_h_constraint_v;
			double c = constraint.value(x);
//Timer.getGlobalTimer("LBF.h.constraint.v_").stop();
//Timer.getGlobalTimer("LBF.h.constraint.g").start();
			double[] gradc = constraint.gradient(x).toArray();
//Timer.getGlobalTimer("LBF.h.constraint.g").stop();
//Timer.getGlobalTimer("LBF.h.constraint.h").start();
			RealMatrix Hc = constraint.hessian(x);
//Timer.getGlobalTimer("LBF.h.constraint.h").stop();

//Timer.getGlobalTimer("LBF.h.comp").start();
			double[][] leftArray = new double[dimension][dimension];
			double c2 = c*c;
			for(int i = 0; i < dimension; ++i) // rows
			{
				for(int j = 0; j < dimension; ++j) // columns
				{
					leftArray[i][j] = gradc[i]*gradc[j]/c2;
				}
			}

			RealMatrix right = Hc == null ? null : Hc.scalarMultiply(-1/c);

			RealMatrix left = new Array2DRowRealMatrix(leftArray, false);
			RealMatrix summand = right == null ? left : left.add(right);
			sum = sum.add(summand);
//Timer.getGlobalTimer("LBF.h.comp").stop();
		}
		
		return sum;
	}

	@Override
	public double value(RealVector x)
	{
		double sum = 0;
		for (TwiceDifferentiableFunction constraint : constraints)
		{
			sum += Math.log(-constraint.value(x));
		}
		return -sum;
	}
	
	public List<TwiceDifferentiableFunction> getConstraints()
	{
		return constraints;
	}
	
	public void addConstraint(TwiceDifferentiableFunction con)
	{
		constraints.add(con);
	}
}