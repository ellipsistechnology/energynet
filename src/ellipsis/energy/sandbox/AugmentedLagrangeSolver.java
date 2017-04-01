package ellipsis.energy.sandbox;

import static ellipsis.util.VectorHelper.normSquared;
import static ellipsis.util.VectorHelper.vector;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.direct.BOBYQAOptimizer;

import com.mls.util.Util;

@SuppressWarnings("deprecation")
public class AugmentedLagrangeSolver
{
	public static interface ScalarFunction
	{
		double value(RealVector x);
	}
	
	public static interface VectorFunction
	{
		RealVector value(RealVector x);
	}
	
	public static interface MatrixFunction
	{
		RealMatrix value(RealVector x);
	}
	
	private ScalarFunction costFunction;
	private VectorFunction costGradient;
	private VectorFunction constraintFunction;
	private MatrixFunction constraintJacobian;
	private VectorFunction projection;
	private double beta = 2;
	private int iterations = 100;
	private int gradDecIterations = 10000;
	private double gradDecTargetGradient = 1e-3;
	private double gradDecStepSize = 1e0;
	private double gamma = 0.5;
	
	public RealVector solve(RealVector initialX)
	{
		double alpha = 1e-6;
		RealVector x = initialX;
		int dimension = x.getDimension()/4;
		RealVector lambda = vector(dimension*2, 0.0);
		double gradTarget = gradDecTargetGradient;
		for(int k = 0; k < iterations; ++k)
		{
if(k == 10)
	Util.nullop();
			RealVector xHat = min(x, lambda, alpha, gradTarget);
			x = projection.value(xHat);
			RealVector g = constraintFunction.value(x);
			
			{
				System.out.print(k+",");
				for(int i = 0; i < x.getDimension(); ++i)
					System.out.print(x.getEntry(i)+",");
				for(int i = 0; i < lambda.getDimension(); ++i)
					System.out.print(lambda.getEntry(i)+",");
				for(int i = 0; i < g.getDimension(); ++i)
					System.out.print(g.getEntry(i)+",");
				RealMatrix gradG = constraintJacobian.value(x);
				RealVector gradGLambda = gradG.operate(lambda);
				RealVector gradC = costGradient.value(x);
				RealVector gradAug = gradG.operate(g);
				RealVector gradL = gradC.add(gradGLambda).add(gradAug.mapMultiply(alpha));
				for(int i = 0; i < gradL.getDimension(); ++i)
					System.out.print(gradL.getEntry(i)+",");
				System.out.print(costFunction.value(x)+",");
				System.out.print(lagrangian(x, lambda, alpha)+",");
				for(int i = 0; i < dimension; ++i)
					System.out.print(",,,"); // skipping t, c and inf, a.k.a. gamma, alpha and influence
				System.out.print(",,"); // skipping ||grad|| and CPU Time
				System.out.print(gradTarget+",");
				System.out.print(gradL.getNorm()+",");
				System.out.println();
			}
			
			lambda = lambda.add(g.mapMultiply(alpha));
			alpha *= beta;
			gradTarget *= gamma;
		}
		return x;
	}
	
	private RealVector min(RealVector x0, RealVector lambda, double alpha, double targetGradient)
	{
		final String FINISHED = "FINISHED";
		int dimension = x0.getDimension();
		org.apache.commons.math3.optimization.direct.BOBYQAOptimizer opt = new BOBYQAOptimizer(2*dimension);
		final RealVector prevX[] = new RealVector[1];
		MultivariateFunction f = new MultivariateFunction()
		{
			double prevL;
			@Override
			public double value(double[] point)
			{
				RealVector x = new ArrayRealVector(point);
				RealVector g = constraintFunction.value(x);
				double L = costFunction.value(x) + lambda.dotProduct(g) + 0.5*alpha*normSquared(g);
//				if(prevX[0] != null)
//				{
//					double grad = Math.abs(prevL-L)/prevX[0].subtract(x).getNorm();
//					if(grad <= targetGradient)
//					{
//						throw new RuntimeException(FINISHED);
//					}
//				}
				prevX[0] = x;
//				prevL = L;
				return L;
			}
		};
		try
		{
			PointValuePair results = opt.optimize(1000, f, GoalType.MINIMIZE, x0.toArray());
			return vector(results.getKey());
		}
		catch(RuntimeException e)
		{
//			if(e.getMessage().equals(FINISHED))
				return prevX[0];
//			e.printStackTrace();
//			return prevX[0];
		}
	}
	
	@SuppressWarnings("unused")
	private RealVector min_old(RealVector x, RealVector lambda, double alpha, double targetGradient)
	{
//		{
//			RealVector _x = new ArrayRealVector(x);
//			double _f = -0.05;
//			_x.setEntry(10, _f);
////			_x.setEntry(11, _f);
//			for(double _e = 0.9; _e <= 1.10001; _e += 0.001)
//			{
//				_x.setEntry(7, _e);
////				_x.setEntry(8, _e);
//				RealVector _g = constraintFunction.value(_x);
//				double _L = costFunction.value(_x) + lambda.dotProduct(_g) + 0.5*alpha*normSquared(_g);
//				RealVector _alphaG = _g.mapMultiply(alpha); // alpha*g(x)
//				RealVector _gradC = costGradient.value(_x); // grad c(x)
//				RealMatrix _gradG = constraintJacobian.value(_x);
//				RealVector _grad = _gradC.add(_gradG.operate(lambda.add(_alphaG))); // grad c(x) + grad g(x)(lambda + alpha*g(x))
//				System.out.println(_e+","+_L+","+_grad.getEntry(7));
////				System.out.println(_e+","+_L+","+_grad.getEntry(8));
//			}
//		}
		
//boolean test = false;
		double gradNorm = Double.MAX_VALUE;
		for(int j = 0; j < gradDecIterations && gradNorm > targetGradient; ++j)
		{
//if(j == 2200)
//{
//	Util.nullop();
//	test = true;
//}
//if(test)
//{
//	RealVector _x = new ArrayRealVector(x);
//	for(double _p = 0.1; _p <= 0.30001; _p += 0.001)
//	{
//		_x.setEntry(2, _p);
//		RealVector _g = constraintFunction.value(_x);
//		double _L = costFunction.value(_x) + lambda.dotProduct(_g) + 0.5*alpha*normSquared(_g);
//		RealVector _alphaG = _g.mapMultiply(alpha); // alpha*g(x)
//		RealVector _gradC = costGradient.value(_x); // grad c(x)
//		RealMatrix _gradG = constraintJacobian.value(_x);
//		RealVector _grad = _gradC.add(_gradG.operate(lambda.add(_alphaG))); // grad c(x) + grad g(x)(lambda + alpha*g(x))
//		System.out.println(_p+","+_L+","+_grad.getEntry(2));
//	}
//}
			
			// Only do grad dec on one dimension at a time:
			int dimension = x.getDimension();
			RealVector grad = null;
			for(int i = dimension/2 /*FIXME DO NOT COMMIT*/; i < dimension; ++i)
			{
//int i = dimension/2 + 1; // FIXME DO NOT COMMIT
				double stepSize = gradDecStepSize;
				RealVector g = constraintFunction.value(x);
				RealVector alphaG = g.mapMultiply(alpha); // alpha*g(x)
				RealVector gradC = costGradient.value(x); // grad c(x)
				RealMatrix gradG = constraintJacobian.value(x);
				grad = gradC.add(gradG.operate(lambda).add(gradG.operate(alphaG))); // grad c(x) + grad g(x)(lambda + alpha*g(x))
				RealVector maskedGrad = new ArrayRealVector(grad);
				maskZeros(maskedGrad, i);
				double f = costFunction.value(x) + lambda.dotProduct(g) + 0.5*alpha*normSquared(g);
				double gradSquared = normSquared(maskedGrad);
				double nextF;
				RealVector nextX;
				do
				{
					nextX = x.subtract(maskedGrad.mapMultiply(stepSize));
					RealVector nextG = constraintFunction.value(nextX);
					nextF = costFunction.value(nextX) + lambda.dotProduct(nextG) + 0.5*alpha*normSquared(nextG);

//					System.out.println(stepSize+","+nextF+","+(f - stepSize*gradSquared)+","+(f - gamma*stepSize*gradSquared));
					stepSize *= gamma;
					
				} while(nextF > f - stepSize*gradSquared); // note that gamma is used one extra time above, so not needed here
				
				x = nextX;
			}
			 
			gradNorm = grad
					.getSubVector(dimension/2, dimension/2) // FIXME DO NOT COMMIT
					.getNorm();
			
//			if(j%10 == 0)
//			{
//				System.out.print(j+","+gradNorm+",");
//				for(int i = 0; i < x.getDimension(); ++i)
//					System.out.print(grad.getEntry(i)+",");
//				System.out.println(lagrangian(x, lambda, alpha));
//			}
		}
		
//		if(gradNorm > targetGradient)
//			System.err.println
////			throw new RuntimeException
//			("Gradient Decent Didn't Converge; gradient was "+gradNorm);
		
		return x;
	}

	public double lagrangian(RealVector x, RealVector lambda, double alpha)
	{
		RealVector g = constraintFunction.value(x);
		return costFunction.value(x) + lambda.dotProduct(g) + 0.5*alpha*normSquared(g);
	}

	private void maskZeros(RealVector v, int i)
	{
		int dimension = v.getDimension();
		for (int j = 0; j < dimension; j++)
		{
			if(j != i)
				v.setEntry(j, 0.0);
		}
	}
	
	
	//// Accessors ////


	public ScalarFunction getCostFunction()
	{
		return costFunction;
	}
	public void setCostFunction(ScalarFunction costFunction)
	{
		this.costFunction = costFunction;
	}
	public VectorFunction getCostGradient()
	{
		return costGradient;
	}
	public void setCostGradient(VectorFunction costGradient)
	{
		this.costGradient = costGradient;
	}
	public VectorFunction getConstraintFunction()
	{
		return constraintFunction;
	}
	public void setConstraintFunction(VectorFunction constraintFunction)
	{
		this.constraintFunction = constraintFunction;
	}
	public MatrixFunction getConstraintJacobian()
	{
		return constraintJacobian;
	}
	public void setConstraintJacobian(MatrixFunction constraintJacobian)
	{
		this.constraintJacobian = constraintJacobian;
	}

	public VectorFunction getProjection()
	{
		return projection;
	}

	public void setProjection(VectorFunction projection)
	{
		this.projection = projection;
	}

	public double getBeta()
	{
		return beta;
	}

	public void setBeta(double beta)
	{
		this.beta = beta;
	}

	public int getIterations()
	{
		return iterations;
	}

	public void setIterations(int iterations)
	{
		this.iterations = iterations;
	}

	public int getGradDecIterations()
	{
		return gradDecIterations;
	}

	public void setGradDecIterations(int gradDecIterations)
	{
		this.gradDecIterations = gradDecIterations;
	}

	public double getGradDecTargetGradient()
	{
		return gradDecTargetGradient;
	}

	public void setGradDecTargetGradient(double gradDecTargetGradient)
	{
		this.gradDecTargetGradient = gradDecTargetGradient;
	}

	public double getGradDecStepSize()
	{
		return gradDecStepSize;
	}

	public void setGradDecStepSize(double gradDecStepSize)
	{
		this.gradDecStepSize = gradDecStepSize;
	}
}
