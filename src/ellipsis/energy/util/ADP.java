package ellipsis.energy.util;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.util.NonParametricApproximator;
import ellipsis.util.Pair;
import ellipsis.util.SampleBasedEstimator;
import ellipsis.util.VectorHelper;

public abstract class ADP implements ADPInterface
{
	private static int nameIndex = 0;
	public String name = "ADP-"+(nameIndex++);
	
	// Parameters:
	public int N_r = 20;
	public int T = 12;
	public double gamma = 1.0;
	private RealVector x_0;

	public CostFunction costFunction;
	
	// Variables:
	public SampleBasedEstimator[] approxJ;
	public Random rand = new Random(nameIndex);

	/**
	 * Trains the estimates of the cost-to-go functions for all times.
	 * @param defaultBandwidth The bandwidth used by the estimators.
	 */
	public void train(double[][] defaultBandwidth)
	{
		// Initialise estimators:
		approxJ = new SampleBasedEstimator[T];
		for(int t = 0; t < T; ++t)
		{
			approxJ[t] = newEstimator(defaultBandwidth[t]);
		}
		
		// Nothing to train if the state is empty:
		if(x_0.getDimension() == 0)
			return;

		// Training iterations:
		int trainingIterations = (T+1)*N_r;

//System.out.println("training ");
		for(int k = 0; k < trainingIterations; ++k)
		{
//System.out.print(k+",");
			// Make sample path for t \in [0, T):
			RealVector[] x = new RealVector[T+1]; // State trajectory.
			x[0] = x_0;
			RealVector[] u = new RealVector[T]; // Control policy.
			RealVector[] x_u = new RealVector[T]; // Post decision ('noiseless') state trajectory.
			for(int t = 0; t < T; ++t)
			{
if(k > 340 && t == 0 && !trained(k, t) && !beforeTraining(k, t))
{
//	System.out.println(((QuadraticEstimator)approxJ[9]).getB());
	Util.nullop();
}
				// Choose random variations:
				RealVector w_t = randomVariations(t);

				// Get optimal control estimate:
				u[t] = combinedControl(k, t, x[t], w_t);
				
				// Get post decision state and next state:
				RealVector wZero = wZero(t);
				x_u[t] = f(t, x[t], u[t], wZero);
				x[t+1] = f(t, x[t], u[t], w_t);
			} // end for t
			
			// Add observation:
			RealVector u_zero = uZero();
			double J = g(T, x[T], u_zero);
			for(int t = T-1; t >= 0; --t)
			{
				// Add sample:
				if(!trained(k, t) && !beforeTraining(k, t))
					addObservation(t, x_u, J);

				// Cost:
				double g = g(t, x[t], u[t]);
				
				// Cost-to-go:
				J = g + gamma*J;
			}
		} // end for k
//System.out.println(".");
	}

	/**
	 * Defaults to creating a {@link NonParametricApproximator}.
	 * @param defaultBandwidth 
	 * @return
	 */
	public SampleBasedEstimator newEstimator(double[] defaultBandwidth)
	{
		NonParametricApproximator npa = new NonParametricApproximator();
		npa.setBandwidth(defaultBandwidth);
		npa.setHistorySize(N_r);
//		npa.setCorrectInterpolation(true);
		return npa;
	}
	
	/**
	 * Cost at time t, when in state x, applying control u.
	 * @param t
	 * @param x
	 * @param u
	 * @return
	 */
	public double g(int t, RealVector x, RealVector u)
	{
		return costFunction.g(this, t, x, u);
	}

	public void addObservation(int t, RealVector[] x_u, double J)
	{
		approxJ[t].addSample(x_u[t], J);
//		BandwidthSelector bs = new BandwidthSelector(approxJ[t]);
//		bs.crossValidate();
	}

	public boolean trained(int k, int t)
	{
		return k > N_r*(T-t);
	}

	public double explorationRate(int k, int t)
	{
		if(beforeTraining(k, t))
			return 1; // all random before training
		if(trained(k, t))
			return 0; // only optimal after training
		
		double rate = T-t-k/(double)N_r; // linearly from random to optimal during training

		if(rate < 0 || rate > 1)
			throw new RuntimeException("Invalid rate "+rate);
		return rate;
	}

	protected boolean beforeTraining(int k, int t)
	{
		return k < N_r*(T-t-1);
	}

	/**
	 * Convenience method for combining a random control with an optimal control.
	 * They will be combined according to the value of {@link #explorationRate(int, int)}.
	 * @param k
	 * @param t
	 * @param uOptimal
	 * @param uRandom
	 * @return
	 */
	public RealVector combineControls(int k, int t, RealVector uOptimal, RealVector uRandom)
	{
		double rho = explorationRate(k, t);
		RealVector uCombined = uRandom.mapMultiply(rho).add(uOptimal.mapMultiply(1-rho));
		return uCombined;
	}

	/**
	 * Calculates the costs-to-go for the given control path.
	 * @param u The control path to apply.
	 * @return The costs-to-go.
	 */
	public double[] ctgs(RealVector[] u)
	{
		double[] ctgs = new double[T+1];
		
		RealVector[] x = new RealVector[T+1];
		x[0] = x_0;
		for(int t = 0; t < T; ++t)
 		{
			RealVector zero = wZero(t);
			x[t+1] = f(t, x[t], u[t], zero);
		}
		
		// Fill in actual costs to go:
		if(ctgs != null)
		{
			double ctg = g(T, x[T], null); // Terminating CTG.
			ctgs[T] = ctg;
			for(int t = T-1; t >= 0; --t)
			{
				double g = g(t, x[t], u[t]);
				ctg += g;
				ctgs[t] = ctg;
			}
		}
		
		return ctgs;
	}
	
	public RealVector[] schedule(double[] ctgs)
	{
		return schedule(ctgs, null, null);
	}

	/**
	 * 
	 * @param ctgs Array of length {@link #T}+1 to be filled with costs-to-go. May be null.
	 * @param x State array of length {@link #T}+1 to be filled with states. May be null.
	 * @return
	 */
	public RealVector[] schedule(double[] ctgs, RealVector[] x, RealVector[] w)
	{
		RealVector[] u = new RealVector[T];
		
		if(x == null)
			x = new RealVector[T+1];
		else if(x.length != T+1)
			throw new RuntimeException("State vector length incorrect. Was "+x.length+", should be "+(T+1));
		
		if(w == null)
			w = new RealVector[T+1];
		else if(w.length != T+1)
			throw new RuntimeException("Noise vector length incorrect. Was "+w.length+", should be "+(T+1));
		
		x[0] = x_0;
		
		for(int t = 0; t < T; ++t)
 		{
			w[t] = wZero(t);
			u[t] = estimateOptimalControl(t, x[t], w[t]);
//			if(u[t] == null)
//				throw new RuntimeException("null control at time "+t);
			RealVector u_ = u[t] == null ? uZero() : u[t]; // FIXME this should not happen
			x[t+1] = f(t, x[t], u_, w[t]);
		}
		
		// Fill in actual costs to go:
		if(ctgs != null)
		{
			double ctg = g(T, x[T], null); // Terminating CTG.
			ctgs[T] = ctg;
			for(int t = T-1; t >= 0; --t)
			{
				double g = g(t, x[t], u[t]);
				ctg += g;
				ctgs[t] = ctg;
			}
		}
		
		return u;
	}
	
	@Override
	public String toString()
	{
		return name;
	}
	
	@Override
	public int getT()
	{
		return T;
	}
	
	@Override
	public RealVector getX0()
	{
		return x_0;
	}
	
	@Override
	public void setNr(int i)
	{
		this.N_r = i;
	}
	
	@Override
	public void setCostFunction(CostFunction g)
	{
		this.costFunction = g;
	}
	
	@Override
	public void setX0(RealVector x_0)
	{
		this.x_0 = x_0;
	}
	
	@Override
	public CostFunction getCostFunction()
	{
		return costFunction;
	}
	
	public int getNr()
	{
		return N_r;
	}
	
	public double getGamma()
	{
		return gamma;
	}

	public String getName()
	{
		return name;
	}

	public void setname(String name)
	{
		this.name = name;
	}

	public void setT(int T)
	{
		this.T = T;
	}
	
	
	//// Debug ////
	
	public Object showEstimations = new Object()
	{
		public String toString() 
		{
try
{
			StringBuffer sb = new StringBuffer();
			int t = 9;

			// Heading:
			for(int i = 0; i < x_0.getDimension(); ++i)
			{
				sb.append("x");
				sb.append(i);
				sb.append(',');
			}
			sb.append("J_t,~J_t,(t="+t+")\n");
			
			// Samples and estimates:
			RealVector min = new ArrayRealVector(x_0.getDimension(), Double.MAX_VALUE);
			RealVector max = new ArrayRealVector(x_0.getDimension(), -Double.MAX_VALUE);
			List<Pair<RealVector, Double>> samples = approxJ[t].getSamplePoints();
			for (Pair<RealVector, Double> sample : samples)
			{
				RealVector x = sample.getKey();
				Double y = sample.getValue();
				double estimate = approxJ[t].value(x);
				for (int i = 0; i < x.getDimension(); i++)
				{
					double x_i = x.getEntry(i);
					sb.append(x_i);
					sb.append(',');
					
					if(x_i < min.getEntry(i))
						min.setEntry(i, x_i);
					else if(x_i > max.getEntry(i))
						max.setEntry(i, x_i);
				}
				sb.append(y);
				sb.append(',');
				sb.append(estimate);
				sb.append('\n');
			}
			
			if(x_0.getDimension() == 1)
			{
				double x_min = min.getEntry(0);
				double x_max = max.getEntry(0);
				double x_step = (x_max-x_min)/20.0;
				sb.append('\n');
				if(x_step != 0)
					for(double x = x_min; x <= x_max; x += x_step)
					{
						double estimate = approxJ[t].value(VectorHelper.vector(x));
						sb.append(x);
						sb.append(",,");
						sb.append(estimate);
						sb.append('\n');
					}
			}
			
			if(x_0.getDimension() == 2)
			{
				double x_min = min.getEntry(0);
				double x_max = max.getEntry(0);
				double x_step = (x_max-x_min)/20.0;
				sb.append('\n');
				if(x_step != 0)
					for(double x = x_min; x <= x_max; x += x_step)
					{
						double estimate = approxJ[t].value(VectorHelper.vector(0.0, x));
						sb.append("0,");
						sb.append(x);
						sb.append(",,");
						sb.append(estimate);
						sb.append('\n');
					}
			}
			
			return sb.toString();
} catch(Exception e) {e.printStackTrace(); return null; }
		}
	};
}