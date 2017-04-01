package ellipsis.energy.sandbox;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.MatrixHelper.matrix;
import static ellipsis.util.VectorHelper.vector;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

public class Sandbox021_AugLagSolver
{
	public static Solution solve(Problem problem, ProblemConfiguration config)
	{
		Solution sol = new Solution(problem);
		sol.x = problem.project(config.initialState); // ensure we start in a valid state
		RealVector lambda = new ArrayRealVector(problem.getEqualityConstraintDimension(), 0.0);
		double epsilon = config.initialEpsilon;
		double alpha = config.initialAlpha;
		Set<Integer> stochasticIndeces = new HashSet<>();
		sol.storeDataPoint(config, lambda, alpha, epsilon, stochasticIndeces);
		for(int k = 0; k < config.iterations; ++k)
		{
if(k == 898)
	Util.nullop();
			if(config.stochastic)
			{
				stochasticIndeces = problem.nextStochasticIndeces(k);
			}

			// Minimise:
			sol.x = minimise(problem, config, sol.x, lambda, alpha, epsilon, stochasticIndeces);
			
//if(k < 1000)
{
			// Step dual variables:
			lambda = bound(lambda.add(problem.g(sol.x).mapMultiply(alpha)), config.lambdaBound);
			
			// Step penalty multiplier:
			if(problem.g(sol.x).getNorm() > config.equalityConstraintTarget)
				alpha = alpha*config.alphaMultiplier;
			
			// Step epsilon:
			epsilon *= config.epsilonMultiplier;
}
//else
//	Util.nullop();
			sol.storeDataPoint(config, lambda, alpha, epsilon, stochasticIndeces);
		}
		
		return sol;
	}
	
	private static double bound(double d, double bound)
	{
		double newD;
		if(Math.abs(d) > bound)
		{
			if(d < 0)
				newD = -bound;
			else
				newD = bound;
		}
		else
			newD = d;
		
		return newD;
	}

	private static RealVector bound(RealVector v, double bound)
	{
		int dimension = v.getDimension();
		RealVector newV = new ArrayRealVector(dimension);
		for (int i = 0; i < dimension; i++)
		{
			double v_i = v.getEntry(i);
			double newV_i = bound(v_i, bound);
			newV.setEntry(i, newV_i);
		}
		
		return newV;
	}

	private static RealVector minimise(Problem problem, ProblemConfiguration config, RealVector initialX, RealVector lambda, double alpha, double epsilon, Set<Integer> stochasticIndeces)
	{
		RealVector x = initialX;
		RealVector grad = lagrangeGradient(problem, config, x, lambda, alpha, stochasticIndeces);
		double improvement;
		double prevLag = lagrange(problem, x, lambda, alpha);
		do
		{
			// Gradient decent step: 
//			double stepSize = backtrack(problem, config, x, grad, lambda, alpha);
			double stepSize = trinarySearch(problem, config, x, grad, lambda, alpha);
//if(lagrange(problem, x.add(grad.mapMultiply(-stepSize2)), lambda, alpha) < lagrange(problem, x.add(grad.mapMultiply(-stepSize)), lambda, alpha))
//	Util.nullop();
			x = x.add(grad.mapMultiply(-stepSize));
			
			// Check for improvement:
			double lag = lagrange(problem, x, lambda, alpha);
//			if(prevLag < lag)
//				throw new RuntimeException("Minimisation failed to reduce Lagrange function: L(x0)="+prevLag+", L(x\')="+lag);
			
			// Project:
			x = problem.project(x);
			grad = lagrangeGradient(problem, config, x, lambda, alpha, stochasticIndeces);
			
			// Check for improvement:
			lag = lagrange(problem, x, lambda, alpha);
			improvement = prevLag - lag;
//			if(improvement < 0) // CAN'T GUARANTEE THIS since backtracking is approximate
//				return initialX;
//				throw new RuntimeException("Projected minimisation failed to reduce Lagrange function: L(x0)="+prevLag+", L(x\')="+lag);
			prevLag = lag;
			
		} while(grad.getNorm() > epsilon && improvement > epsilon);
		
		return x;
	}
	
	/*
	     StringBuffer sb = new StringBuffer();
	     RealVector xp = problem.project(x);
	     for(int i = 0; i < 12; ++i)
	     {
	         sb.append(x.getEntry(i) == xp.getEntry(i) ? 
	         	"x"+i+"\n" : 
	         	"[x"+i+x.getEntry(i)+"->"+xp.getEntry(i)+"]\n");
	     }
	     
	     StringBuffer sb = new StringBuffer();
	     RealVector _x = new ArrayRealVector(x);
	     for(double x4 = x.getEntry(4)-1e-4; x4 < x.getEntry(4)+1e-4; x4 += 1e-5)
	     {
	     	 _x.setEntry(4, x4);
		     for(double x5 = x.getEntry(5)-1e-4; x5 < x.getEntry(5)+1e-4; x5 += 1e-5)
		     {
		     	_x.setEntry(5, x5);
		     	sb.append(lagrange(problem, _x, lambda, alpha)+",");
		     }
		     sb.append("\n");
	     }
	 */
	
	private static double trinarySearch(Problem problem, ProblemConfiguration config, RealVector x, RealVector grad, RealVector lambda, double alpha)
	{
		double l0 = lagrange(problem, x, lambda, alpha);
		
		double t_left = 0;
		double t_right = 1.0;
		
		while(true)
		{
			// Lagrange values at outer points:
			double l_left = lagrange(problem, x.add(grad.mapMultiply(-t_left)), lambda, alpha);
			double l_right = lagrange(problem, x.add(grad.mapMultiply(-t_right)), lambda, alpha);
			
			// Stop if close enough to the minimum:
			if(t_left == t_right || l_left <= l0 && l_right <= l0 && Math.abs(l_left-l_right) < 1e-24)
				return (t_left+t_right)/2;
	
			// Evaluate intermittant points:
			double t1 = t_left + (t_right-t_left)/3;
			double t2 = t_left + 2*(t_right-t_left)/3;
			double l1 = lagrange(problem, x.add(grad.mapMultiply(-t1)), lambda, alpha);
			double l2 = lagrange(problem, x.add(grad.mapMultiply(-t2)), lambda, alpha);
			
			if(l_left <= l1) // min is between left and t1
			{
				t_right = t1;
			}
			else if(l_right < l2) // expand the search
			{
				t_left = t2;
				t_right = 2*t_right;
			}
			else if(l1 <= l_left && l1 <= l2) // min is between left and t2
			{
				t_right = t2;
			}
			else if(l2 <= l1 && l2 <= l_right) //min is between t1 and right
			{
				t_left = t1;
			}
			else
				throw new RuntimeException("Impossible state reached during trinary search.");
		}
	}

	/* plot(problem, x, lambda, alpha) */
	@SuppressWarnings("unused")
	private static String plot(Problem problem, RealVector x, RealVector lambda, double alpha)
	{
		StringBuffer sb = new StringBuffer();
		double x0 = x.getEntry(0);
		double x1 = x.getEntry(1);
		
		double inc = 0.03;
		double min = 0.0;
		double max = 0.60;
		for(double i = min; i < max; i += inc)
		{
			for(double j = min; j < max; j += inc)
			{
				sb.append(lagrange(problem, vector(x0+j, x1+i), lambda, alpha));
				sb.append(",");
			}
			sb.append("\n");
		}
		
		return sb.toString();
	}

	public static RealVector lagrangeGradient(Problem problem, ProblemConfiguration config, RealVector x, RealVector lambda, double alpha, Set<Integer> stochasticIndeces)
	{
		RealVector grad = lagrangeGradientRaw(problem, x, lambda, alpha);
		
		// If stochastic only use one agent (2-d step):
		if(config.stochastic && stochasticIndeces != null)
		{
			int dimension = grad.getDimension();
			RealVector stochasticGrad = vector(dimension, 0.0);
			for (int i = 0; i < dimension; i++)
			{
				if(stochasticIndeces.contains(i))
					stochasticGrad.setEntry(i, grad.getEntry(i));
			}
			return stochasticGrad;
		}
		else
		{
			return grad;
		}
	}

	public static RealVector lagrangeGradientRaw(Problem problem, RealVector x,
			RealVector lambda, double alpha)
	{
		RealVector costGradient = problem.costGradient(x);
		RealMatrix equalityConstraintGradient = problem.equalityConstraintGradient(x);
		RealVector g = problem.g(x);
		RealVector grad = costGradient.add(
		                  equalityConstraintGradient.operate(lambda.add(g.mapMultiply(alpha))));
		return grad;
	}

	private static double backtrack(Problem problem, ProblemConfiguration config, RealVector x, RealVector grad, RealVector lambda, double alpha)
	{
		double lagrange = lagrange(problem, x, lambda, alpha);
		double stepSize = 1;
		double grad2 = grad.dotProduct(grad);
		RealVector newX = x.add(grad.mapMultiply(-stepSize));
		while(lagrange(problem, newX, lambda, alpha) > lagrange - stepSize*config.stepScale*grad2)
		{
//System.out.println(stepSize+","+lagrange(problem, newX, lambda, alpha)+","+(lagrange - stepSize*grad2));
			stepSize *= config.stepMultiplier;
			newX = x.add(grad.mapMultiply(-stepSize));
		}

//for(double d = 0; d > -1e-8; d -= 1e-9)
//{
//	System.out.println(d+","+lagrange(problem, newX, lambda, alpha)+","+(lagrange - d*grad2));
//	newX = x.add(grad.mapMultiply(-d));
//}
		return stepSize;
	}

	public static double lagrange(Problem problem, RealVector x, RealVector lambda, double alpha)
	{
		RealVector g = problem.g(x);
		return problem.cost(x) 
			 + g.dotProduct(lambda) 
			 + alpha*g.dotProduct(g)/2;
	}
	
	public static class Solution
	{
		private Problem problem;
		
		public RealVector x;
		public ArrayList<RealVector> xs = new ArrayList<>();
		public ArrayList<Double> lagrangeValues = new ArrayList<>();
		public ArrayList<Double> costValues = new ArrayList<>();
		public ArrayList<RealVector> gradientValues = new ArrayList<>();
		public ArrayList<RealVector> gValues = new ArrayList<>();
		public ArrayList<Double> epsilonValues = new ArrayList<>();
		public ArrayList<Double> alphaValues = new ArrayList<>();
		public ArrayList<RealVector> lambdaValues = new ArrayList<>();
		
		public Solution(Problem problem)
		{
			this.problem = problem;
		}

		public void storeDataPoint(ProblemConfiguration config, RealVector lambda, double alpha, double epsilon, Set<Integer> stochasticIndeces)
		{
			xs.add(x);
			lagrangeValues.add(lagrange(problem, x, lambda, alpha));
			costValues.add(problem.cost(x));
			gradientValues.add(lagrangeGradientRaw(problem, x, lambda, alpha));
			gValues.add(problem.g(x));
			epsilonValues.add(epsilon);
			alphaValues.add(alpha);
			lambdaValues.add(lambda);
		}
		
		public int size()
		{
			return xs.size();
		}
	}

	public static class ProblemConfiguration
	{
		// Initialisation:
		public double initialAlpha;
		public double initialEpsilon;
		public RealVector initialState;
		
		// Iteration parameters:
		public int iterations;
		public double alphaMultiplier;
		public double epsilonMultiplier;
		public double lambdaBound;
		public double equalityConstraintTarget = 1e-3; // alpha increase will only proceed while ||g(x)|| > target
		public boolean stochastic = false;
		
		// Backtracking parameters:
		public double stepMultiplier; // stepSize *= stepMultiplier at each iteration
		public double stepScale; // scales the slope for approximating the minimum
	}

	public static interface Problem
	{
		int getEqualityConstraintDimension();

		default Set<Integer> nextStochasticIndeces(int k) { return null; }

		/**
		 * Projects the given point to a valid point according to a closed convex constraint set.
		 * @param x The initial state.
		 * @return The projected state.
		 */
		RealVector project(RealVector x);

		// Problem functions:
		double cost(RealVector x);
		RealVector g(RealVector x);

		// Function gradients:
		RealVector costGradient(RealVector x);
		RealMatrix equalityConstraintGradient(RealVector x);
	}
	
	
	//// TEST CASES ////
	
	public static void main(String[] args)
	{
		Problem p = testProblem001();
		ProblemConfiguration config = testConfig001();
		
		Solution sol = solve(p, config);
		
//		printHumanReadable(p, sol);
		printCSV(p, sol, 1000);
	}

	public static void printHumanReadable(Problem problem, Solution sol)
	{
		int i = 0;
		for (RealVector x : sol.xs)
		{
			System.out.println(
					"x=["+x.getEntry(0)+","+x.getEntry(1)+"] "
				  + "L="+sol.costValues.get(i)
				  + " grad_L="+sol.gradientValues.get(i).getNorm()
				  + " g(x)="+sol.gValues.get(i).getEntry(0)
				  + " epsilon="+sol.epsilonValues.get(i)
				  + " alpha="+sol.alphaValues.get(i));
			++i;
		}
	}
	
	public static void printCSV(Problem problem, Solution sol, int lineCount)
	{
		PrintStream out = System.out;
		int dimension = sol.x.getDimension();
		int gDimension = problem.getEqualityConstraintDimension();
		int logFrequency = sol.size()/lineCount;
		
		// Header:
		out.print("k,");
		for (int i = 0; i < dimension; i++)
		{
			out.print("x"+i+",");
		}
		out.print("c(x),L,");
		for (int i = 0; i < gDimension; i++)
		{
			out.print("g"+i+"(x),");
		}
		out.print("alpha,epsilon,");
		for (int i = 0; i < dimension; i++)
		{
			out.print("grad_x"+i+" L(x),");
		}
		for (int i = 0; i < gDimension; i++)
		{
			out.print("lambda"+i+",");
		}
		out.println();
		
		// Data:
		int k = 0;
		for (RealVector x : sol.xs)
		{
			if(k%logFrequency == 0)
			{
				out.print(k+",");
				for (int i = 0; i < dimension; i++)
				{
					out.print(x.getEntry(i)+",");
				}
				out.print(sol.costValues.get(k)+",");
				out.print(sol.lagrangeValues.get(k)+",");
				RealVector g = sol.gValues.get(k);
				for (int i = 0; i < gDimension; i++)
				{
					out.print(g.getEntry(i)+",");
				}
				out.print(sol.alphaValues.get(k)+","+sol.epsilonValues.get(k)+",");
				RealVector grad = sol.gradientValues.get(k);
				for (int i = 0; i < dimension; i++)
				{
					out.print(grad.getEntry(i)+",");
				}
				RealVector lambda = sol.lambdaValues.get(k);
				for (int i = 0; i < gDimension; i++)
				{
					out.print(lambda.getEntry(i)+",");
				}
				
				out.println();
			}
			
			++k;
		}
	}

	public static ProblemConfiguration testConfig001()
	{
		ProblemConfiguration config = new ProblemConfiguration();
		
		config.initialAlpha = 1e-0;
		config.alphaMultiplier = 1.05;
		config.initialEpsilon = 10;
		config.epsilonMultiplier = 0.99;
		config.iterations = 1000;
		config.initialState = vector(0.6, 0.3);//Math.random(), Math.random());
		config.stepMultiplier = 0.5;
		config.stepScale = 0.5;
		config.lambdaBound = 2.0;
		config.stochastic = true;

		return config;
	}

	public static Problem testProblem001()
	{
		return new Problem()
		{
			@Override
			public int getEqualityConstraintDimension()
			{
				return 1;
			}

			@Override
			public double cost(RealVector x)
			{
				double x0 = x.getEntry(0);
				double x1 = x.getEntry(1);
				return x0*x0 + x1*x1;
			}

			@Override
			public RealVector g(RealVector x)
			{
				return vector(x.getEntry(0) - x.getEntry(1) - 1.0);
			}

			@Override
			public RealVector costGradient(RealVector x)
			{
				return vector(2*x.getEntry(0), 2*x.getEntry(1));
			}

			@Override
			public RealMatrix equalityConstraintGradient(RealVector x)
			{
				return matrix(
					array(1.0),
					array(-1.0)
				);
			}

			@Override
			public RealVector project(RealVector x)
			{
				double x0 = x.getEntry(0);
				double x1 = x.getEntry(1);
				if(x0+x1+1 > 0)
				{
					return vector((-x1+x0-1)/2, (x1-x0-1)/2);
				}
				return x;
			}
			
			Random rand = new Random(0);
			@Override
			public Set<Integer> nextStochasticIndeces(int k)
			{
				Set<Integer> is = new HashSet<>();
				if(rand.nextBoolean())
					is.add(1);
				else
					is.add(0);
				return is;
			}
		};
	}
}
