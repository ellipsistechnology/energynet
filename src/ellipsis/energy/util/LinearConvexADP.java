package ellipsis.energy.util;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Logger;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import com.joptimizer.optimizers.OptimizationResponse;

import ellipsis.util.BarrierFunction;
import ellipsis.util.HessianFunction;
import ellipsis.util.SampleBasedEstimator;
import ellipsis.util.TwiceDifferentiableFunction;

/**
 * Implements a linear state transition function in the form f(x) = Ax + Bu + Cw,
 * and assumes that the cost function is convex.
 * @author bmillar
 *
 */
public abstract class LinearConvexADP extends ADP
{
	// Parameters:
	protected RealMatrix A, B, C;
	protected double[][] noiseMean; // [time][unit]
	protected double[][] noiseSD; // [time][unit]
	protected double[] zeroControl; // [unit]
//	protected Set<TwiceDifferentiableFunction> constraints;
	protected BarrierFunction barrier; // FIXME change this to be a set of ConvexMultivariateRealFunction
	
	// Variables:
	private Random rand = new Random(0);

	protected abstract RealVector randomControl(int t, RealVector x_t, RealVector w_t);
	
	@Override
	public RealVector f(int t, RealVector x_t, RealVector u_t, RealVector w_t)
	{
		return A.operate(x_t).add(B.operate(u_t)).add(C.operate(w_t));
	}

	private RealVector[] wZeros;
	/**
	 * @return noiseMean for all entries.
	 */
	@Override
	public RealVector wZero(int time)
	{
		if(wZeros == null)
		{
			wZeros = new RealVector[getT()+1];
			for (int t = 0; t < getT(); t++)
			{
				wZeros[t] = new ArrayRealVector(noiseMean[t]);
			}
			wZeros[getT()] = new ArrayRealVector(noiseMean[0].length);
		}
		return wZeros[time];
	}

	@Override
	public RealVector uZero()
	{
		return new ArrayRealVector(zeroControl);
	}

	@Override
	public RealVector randomVariations(int t)
	{
		int dimension = noiseMean[t].length;
		RealVector w_t = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
		{
			w_t.setEntry(i, noiseMean[t][i]+rand.nextGaussian()*noiseSD[t][i]);
		}
		
		return w_t;
	}
	
public static int converged = 0;
public static int notConverged = 0;
	@Override
	public RealVector estimateOptimalControl(int t, RealVector x_t, RealVector w_t)
	{
		// Constraints:
		List<TwiceDifferentiableFunction> constraints = barrier.getConstraints();
		ConvexMultivariateRealFunction[] inequalities = constraints.toArray(new ConvexMultivariateRealFunction[constraints.size()]);
		
		// Set times and states:
		for (TwiceDifferentiableFunction constraint : constraints)
		{
			if(constraint instanceof TimeDependent)
				((TimeDependent)constraint).setTime(t);
			if(constraint instanceof StateDependent)
				((StateDependent)constraint).setState(x_t);
		}
		
		if(costFunction instanceof TimeDependent)
			((TimeDependent)costFunction).setTime(t);
		if(costFunction instanceof StateDependent)
			((StateDependent)costFunction).setState(x_t);
		
		SampleBasedEstimator approxJ_t = approxJ[t];
		if(approxJ_t instanceof TimeDependent)
			((TimeDependent)approxJ_t).setTime(t);
		if(approxJ_t instanceof StateDependent)
			((StateDependent)approxJ_t).setState(x_t);
		
		// Objective:
		CostToGoFunction ctg = new CostToGoFunction();
		if(!(approxJ_t instanceof HessianFunction))
			throw new RuntimeException(approxJ_t.getClass()+" is not a HessianFunction at time "+t);
		ctg.approx = (HessianFunction) approxJ_t;
		ctg.t = t;
		ctg.x_t = x_t;
		if(!(costFunction instanceof HessianFunction))
			throw new RuntimeException(costFunction.getClass()+" is not a HessianFunction");
		ctg.g = (HessianFunction) costFunction;
		
		// Optimize:
		OptimizationRequest or = new OptimizationRequest();
		or.setF0(ctg);
		or.setFi(inequalities);
		double[] uZero = uZero().toArray();
		or.setInitialPoint(uZero);

		JOptimizer opt = new JOptimizer();
		opt.setOptimizationRequest(or);
		try
		{
			int returnCode = opt.optimize();
			if(returnCode != OptimizationResponse.SUCCESS)
				throw new RuntimeException("Failed to optimize.");
			double[] u = opt.getOptimizationResponse().getSolution();
			++converged;
			return new ArrayRealVector(u);
		}
		catch (Exception e)
		{
			Logger.getLogger(LinearConvexADP.class).warn("Failed to find optimal control", e);
			++notConverged;
			String msg = e.getMessage();
			if("initial point must be strictly feasible".equals(msg))
				return uZero();
//			if("Max iterations limit reached".equals(msg) || "KKT solution failed".equals(msg))
				return null;
//			else
//				throw new RuntimeException(e);
		}
	}

//	@Override
//	public RealVector estimateOptimalControl(int t, RealVector x_t, RealVector w_t)
//	{
//		// TODO Analytically find unconstrained optimum first, 
//		// and return it if no constraints are breached:
////Timer.getGlobalTimer("estimateOptimalControl").start();
//
//		// Barrier Method:
//		BarrierSolver bs = new BarrierSolver();
//		List<TwiceDifferentiableFunction> constraints = barrier.getConstraints();
//		
//		for (TwiceDifferentiableFunction constraint : constraints)
//		{
//			if(constraint instanceof TimeDependent)
//				((TimeDependent)constraint).setTime(t);
//		}
//		if(barrier instanceof TimeDependent)
//			((TimeDependent)barrier).setTime(t);
//		bs.setConstraintDimension(constraints.size()+1);
//		
//		bs.setBarrier(barrier);
//		bs.setStepSize(1.0);
//		bs.setInitialZ(10.0);
//		bs.setEpsilon(1e-5);
//		
//		CostToGoFunction ctg = new CostToGoFunction();
//		if(!(approxJ[t] instanceof HessianFunction))
//			throw new RuntimeException(approxJ[t].getClass()+" is not a HessianFunction at time "+t);
//		if(approxJ[t] instanceof TimeDependent)
//			((TimeDependent)approxJ[t]).setTime(t);
//		ctg.approx = (HessianFunction) approxJ[t];
//		ctg.t = t;
//		ctg.x_t = x_t;
//		if(!(costFunction instanceof HessianFunction))
//			throw new RuntimeException(costFunction.getClass()+" is not a HessianFunction");
//		if(costFunction instanceof TimeDependent)
//			((TimeDependent)costFunction).setTime(t);
//		ctg.g = (HessianFunction) costFunction;
//		bs.setFunction(ctg);
//		
//		RealVector uZero = uZero();
//		RealVector u = bs.solve(uZero);
//		
//		if(u.isNaN() || u.isInfinite())
//			throw new RuntimeException("Invalid control selected: "+u);
////Timer.getGlobalTimer("estimateOptimalControl").stop();
//		return u;
//	}

//	@Override
//	public RealVector estimateOptimalControl(int t, RealVector x_t, RealVector w_t)
//	{
//		final List<TwiceDifferentiableFunction> constraints = barrier.getConstraints();
//		
//		// Pass time and state to functions that need it:
//		for (TwiceDifferentiableFunction constraint : constraints)
//		{
//			if(constraint instanceof TimeDependent)
//				((TimeDependent)constraint).setTime(t);
//			if(constraint instanceof StateDependent)
//				((StateDependent)constraint).setState(x_t);
//		}
//		
//		if(approxJ[t] instanceof TimeDependent)
//			((TimeDependent)approxJ[t]).setTime(t);
//		if(approxJ[t] instanceof StateDependent)
//			((StateDependent)approxJ[t]).setState(x_t);
//
//		if(costFunction instanceof TimeDependent)
//			((TimeDependent)costFunction).setTime(t);
//		if(costFunction instanceof StateDependent)
//			((StateDependent)costFunction).setState(x_t);
//		
//		// Setup cost-go-go for minimising:
//		final CostToGoFunction ctg = new CostToGoFunction();
//		if(!(approxJ[t] instanceof HessianFunction))
//			throw new RuntimeException(approxJ[t].getClass()+" is not a HessianFunction at time "+t);
//		ctg.approx = (HessianFunction) approxJ[t];
//		ctg.t = t;
//		ctg.x_t = x_t;
//		if(!(costFunction instanceof HessianFunction))
//			throw new RuntimeException(costFunction.getClass()+" is not a HessianFunction");
//		ctg.g = (HessianFunction) costFunction;
//		
//		// Setup solver:
//		// TODO Implement this without wrapping CostToGoFunction.
//		Objective f = new Objective()
//		{
//			@Override
//			public double value(Matrix u)
//			{
//				return ctg.value(toVector(u));
//			}
//
//			@Override
//			public Matrix hessian(Matrix u)
//			{
//				return toMatrix(ctg.hessian(toVector(u)));
//			}
//
//			@Override
//			public Matrix gradient(Matrix u)
//			{
//				return toMatrix(ctg.gradient(toVector(u)));
//			}
//		};
//		Constraints F = new Constraints()
//		{
//			@Override
//			public Matrix value(Matrix x)
//			{
//				int m = constraints.size();
//				double[][] F = new double[m][1];
//				for(int i = 0; i < m; ++i)
//				{
//					F[i][0] = constraints.get(i).value(toVector(x));
//				}
//				
//				return new DenseMatrix(F);
//			}
//			
//			@Override
//			public Matrix derivative(Matrix x)
//			{
//				int m = constraints.size();
//				int n = x.getRowDimension();
//				double[][] DF = new double[m][n];
//				for(int i = 0; i < m; ++i)
//				{
//					DF[i] = constraints.get(i).gradient(toVector(x)).toArray();
//				}
//				
//				return new DenseMatrix(DF);
//			}
//
//			@Override
//			public Matrix hessian(int i, Matrix x)
//			{
//				RealMatrix hessian = constraints.get(i).hessian(toVector(x));
//				if(hessian == null)
//				{
//					int dim = x.getRowDimension();
//					return new DenseMatrix(dim, dim);
//				}
//				return toMatrix(hessian);
//			}
//
//			@Override
//			public int constrantCount()
//			{
//				return constraints.size();
//			}
//		};
//		Matrix u = toMatrix(uZero());
//		int dim = u.getRowDimension();
//		PrimalDualInteriorPointSolver solver = new PrimalDualInteriorPointSolver(f, F, zero(1, dim), zero(1, 1));
//		boolean converged = solver.solve(u);
//if(!converged && t == 9)
//	Util.nullop();
//		
//		return toVector(solver.getX());
//	}
//
//	private Matrix zero(int rows, int cols)
//	{
//		return new SparseMatrix(rows, cols);
//	}
//
//	private Matrix toMatrix(RealVector v)
//	{
//		int dimension = v.getDimension();
//		double[][] m = new double[dimension][1];
//		for(int i = 0; i < dimension; ++i)
//		{
//			m[i][0] = v.getEntry(i);
//		}
//		return new DenseMatrix(m);
//	}
//	
//	private RealVector toVector(Matrix x)
//	{
//		int dimension = x.getRowDimension();
//		double[][] xArray = x.getData();
//		double[] v = new double[dimension];
//		for(int i = 0; i < dimension; ++i)
//		{
//			v[i] = xArray[i][0];
//		}
//		return new ArrayRealVector(v);
//	}
//	
//	private Matrix toMatrix(RealMatrix rm)
//	{
//		double[][] m = new double[rm.getRowDimension()][rm.getColumnDimension()];
//		for(int i = 0; i < rm.getRowDimension(); ++i)
//		{
//			for(int j = 0; j < rm.getColumnDimension(); ++j)
//			{
//				m[i][j] = rm.getEntry(i, j);
//			}
//		}
//		return new DenseMatrix(m);
//	}

	@Override
	public RealVector combinedControl(int k, int t, RealVector x_t, RealVector w_t)
	{
		// Approximate optimal control:
		RealVector uOptimal;
//		if(explorationRate(k, t) < 1.0)
			uOptimal = estimateOptimalControl(t, x_t, w_t); // TODO improve efficiency by not calling this if only random will be used anyway
//		else
//			uOptimal = null;

		// Choose random control:
		RealVector uRandom = randomControl(t, x_t, w_t);
		
		// Combine optimal and random control: 
		if(uOptimal == null)
			return uRandom;
		RealVector uCombined = combineControls(k, t, uOptimal, uRandom);

		return uCombined;
	}
	
	@Override
	public void train(double[][] defaultBandwidth)
	{
		super.train(defaultBandwidth);

//		checkConstraints(); // FIXME would be better if this could stay here
	}

	@SuppressWarnings("unused")
	private void checkConstraints()
	{
		// Check constraints:
		RealVector[] x = new RealVector[T+1];
		RealVector[] w = new RealVector[T+1];
		RealVector[] u = schedule(null, x, w);
		for(int t = 0; t < T; ++t)
		{
			List<TwiceDifferentiableFunction> inequalities = barrier.getConstraints();
			for (TwiceDifferentiableFunction con : inequalities)
			{
				if(u[t] == null)
					continue; //throw new RuntimeException("null control at time "+t);
				if(con instanceof TimeDependent)
					((TimeDependent)con).setTime(t);
				if(con instanceof StateDependent)
					((StateDependent)con).setState(x[t]);
				if(((ConvexMultivariateRealFunction)con).value(u[t].toArray()) > 0)
					throw new RuntimeException("Invalid solution found for adp '"+this.toString()+"': t="+t+", x[t]="+x[t]+", u[t]="+u[t]+", w[t]="+w[t]);
			}
		}
	}
	
	public class CostToGoFunction implements TwiceDifferentiableFunction, ConvexMultivariateRealFunction
	{
		private HessianFunction g; // with respect to u
		private HessianFunction approx;
		private RealVector x_t;
		private int t;

		@Override
		public RealVector gradient(RealVector u)
		{
			RealVector xu = f(t, x_t, u, wZero(t)); // Post decision state.
			
			// Gradient of g(x,u):
			RealVector gradg = g.gradient(u);
			
			// Gradient of ~V(x(u)) with respect to u:
			// grad_u ~V(x(u)) = B^T[grad_x ~V(x(u))]
			// Given that x(u) = Ax + Bu
			RealVector gradV = B.transpose().operate(approx.gradient(xu));
			
			// Gradient is grad_u g(x,u) + grad_u ~V(x(u)):
			return gradg == null ? gradV : 
				   gradV == null ? gradg :
					               gradg.add(gradV);
		}

		@Override
		public RealMatrix hessian(RealVector u)
		{
			RealVector xu = f(t, x_t, u, wZero(t)); // Post decision state.
			
			// H(g(x,u)):
			RealMatrix Hg = g.hessian(u);
			
			// H(~V(x(u)):
			RealMatrix HV = B.transpose().multiply(approx.hessian(xu).multiply(B));
			
			// Hessian is H(g(x,u)) + H(~V(x(u))):
			return Hg == null ? HV :
				   HV == null ? Hg :
					            Hg.add(HV);
		}

		@Override
		public double value(RealVector u)
		{
			double cost = ((TwiceDifferentiableFunction)g).value(u);
			RealVector w_t = LinearConvexADP.this.wZero(t);
			RealVector x_u = LinearConvexADP.this.f(t, x_t, u, w_t);
			double ctg = ((TwiceDifferentiableFunction)approx).value(x_u);
			return cost + ctg;
		}
		
		
		//// J

		@Override
		public double value(double[] u)
		{
			double cost = ((ConvexMultivariateRealFunction)g).value(u);
			RealVector w_t = LinearConvexADP.this.wZero(t);
			RealVector x_u = LinearConvexADP.this.f(t, x_t, new ArrayRealVector(u), w_t);
			double ctg = ((TwiceDifferentiableFunction)approx).value(x_u);
			return cost + ctg;
		}

		@Override
		public double[] gradient(double[] u)
		{
			return gradient(new ArrayRealVector(u)).toArray();
		}

		@Override
		public double[][] hessian(double[] u)
		{
			return hessian(new ArrayRealVector(u)).getData();
		}

		@Override
		public int getDim()
		{
			return ((ConvexMultivariateRealFunction)g).getDim();
		}
	}
	
	
	//// Accessors ////

	public RealMatrix getA()
	{
		return A;
	}

	public void setA(RealMatrix a)
	{
		A = a;
	}

	public RealMatrix getB()
	{
		return B;
	}

	public void setB(RealMatrix b)
	{
		B = b;
	}

	public RealMatrix getC()
	{
		return C;
	}

	public void setC(RealMatrix c)
	{
		C = c;
	}

	public double[][] getNoiseMean()
	{
		return noiseMean;
	}

	public void setNoiseMean(double[][] noiseMean)
	{
		this.noiseMean = noiseMean;
	}

	public double[][] getNoiseSD()
	{
		return noiseSD;
	}

	public void setNoiseSD(double[][] noiseSD)
	{
		this.noiseSD = noiseSD;
	}

	public double[] getZeroControl()
	{
		return zeroControl;
	}

	public void setZeroControl(double[] zeroControl)
	{
		this.zeroControl = zeroControl;
	}

	public BarrierFunction getBarrier()
	{
		return barrier;
	}

	public void setBarrier(BarrierFunction barrier)
	{
		this.barrier = barrier;
	}
}