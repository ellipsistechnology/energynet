package ellipsis.energy.util;

import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public abstract class DiscreteADP extends ADP
{
	/**
	 * Called for each possible control when minimising (see {@link #argMin(int, RealVector, Set, RealVector)}). 
	 * Allows for control dependent noise.
	 * @param u The control being considered as possible optimal.
	 * @param w The noise to be modified.
	 */
	public abstract void prepareNoise(RealVector u, RealVector w);
	
	/**
	 * Indicates if a particular control and noise combination causes a constraint breach.
	 * Used by {@link #argMin(int, RealVector, Set, RealVector)} to discard candidate controls.
	 * @param t Time.
	 * @param u Control.
	 * @param w Noise.
	 * @return True if this control noise combination causes a constraint breach.
	 */
	public abstract boolean breached(int t, RealVector u, RealVector w);
	
	/**
	 * Provides a list of all possible controls.
	 * @param t
	 * @param x_t
	 * @return
	 */
	public abstract Set<RealVector> enumerateControls(int t, RealVector x_t);
	
	/**
	 * Estimates the optimal control from the set of viable controls given by {@link #enumerateControls(int, RealVector)}.
	 */
	@Override
	public RealVector combinedControl(int k, int t, RealVector x_t, RealVector w_t)
	{
		// Enumerate controls for t \in [0, T):
		Set<RealVector> U = enumerateControls(t, x_t);

		// Approximate optimal control:
		RealVector uOptimal = argMin(t, x_t, U, w_t);

		// Choose random control:
		RealVector uRandom = randomControl(U);
		
		// Combine optimal and random control: 
		RealVector uCombined = combineControls(k, t, uOptimal, uRandom);

		return uCombined;
	}

	/**
	 * Randomly selects a control from the given set of viable controls.
	 * @param U
	 * @return
	 */
	public RealVector randomControl(Set<RealVector> U)
	{
		RealVector uRandom = null;
		int r = rand.nextInt(U.size());
		int i = 0;
		for (RealVector _u : U)
		{
			if(i == r)
			{
				uRandom = _u;
				break;
			}
			++i;
		}
		return uRandom;
	}
	
	/**
	 * Searches each of the viable controls and returns the one that gives the minumal cost-to-go
	 * according to {@link #g(int, RealVector, RealVector)} and the estimators approxJ.
	 * @param t
	 * @param x_t
	 * @param U
	 * @param w_t
	 * @return
	 */
	public RealVector argMin(int t, RealVector x_t, Set<RealVector> U, RealVector w_t)
	{
		RealVector u_min = null;
		double ctg_min = Double.MAX_VALUE;
		RealVector u_max_breach = null;
		double ctg_max_breach = -Double.MAX_VALUE;

		for (RealVector u : U)
		{
			// Zero out unused portions of DG noise (for zero output DGs):
			RealVector w = new ArrayRealVector(w_t);
			prepareNoise(u, w);

			// Approximate CTG based on post decision state:
			RealVector x_u = f(t, x_t, u, w);
			double ctg_next = approxJ[t].value(x_u);
			double g = g(t, x_t, u);
			double ctg = g + gamma*ctg_next;
			
			if(Double.isNaN(ctg_next))
				throw new RuntimeException("ctg_next was NaN");

			// If constraints are breached then keep this solution separate:
			if(breached(t, u, w))
			{
				// This assumes a decreasing CTG function:
				if(ctg > ctg_max_breach)
				{
					ctg_max_breach = ctg;
					u_max_breach = u;
				}
			}
			else
			{
				if(ctg < ctg_min)
				{
					ctg_min = ctg;
					u_min = u;
				}
			}
		}

		if(u_min == null)
			return u_max_breach;
		else
			return u_min;
	}

	@Override
	public RealVector estimateOptimalControl(int t, RealVector x_t, RealVector w_t)
	{
		Set<RealVector> U = enumerateControls(t, x_t);
		RealVector u_t = argMin(t, x_t, U, w_t);
		return u_t;
	}
}