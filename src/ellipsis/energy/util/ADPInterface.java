package ellipsis.energy.util;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.SampleBasedEstimator;

public interface ADPInterface extends LocalController
{
	/**
	 * 
	 * @param ctgs Array of length {@link #T}+1 to be filled with costs-to-go. May be null.
	 * @param x State array of length {@link #T}+1 to be filled with states. May be null.
	 * @return
	 */
	public RealVector[] schedule(double[] ctgs, RealVector[] x, RealVector[] w);
	
	/**
	 * Next state.
	 * @param x_t
	 * @param u_t
	 * @param w_t
	 * @return
	 */
	public RealVector f(int t, RealVector x_t, RealVector u_t, RealVector w_t);
	
	/**
	 * Mean noise at time t.
	 * @param t
	 * @return
	 */
	public RealVector wZero(int t);
	
	/**
	 * A sample random variation for time t (w_t).
	 * @param t
	 * @return
	 */
	public RealVector randomVariations(int t);
	
	/**
	 * Estimates the optimal control for the given time, state and noise (minimises the cost-to-go).
	 * @param t
	 * @param x_t
	 * @param w_t
	 * @return
	 */
	public RealVector estimateOptimalControl(int t, RealVector x_t, RealVector w_t);
	
	/**
	 * Estimate the optimal control.
	 * @param k Training iteration.
	 * @param t Time.
	 * @param x_t State.
	 * @param w_t Noise.
	 * @return
	 */
	public RealVector combinedControl(int k, int t, RealVector x_t, RealVector w_t);
	
	/**
	 * Convenience method for the central controller to create the appropriate class
	 * for distributed ADPs.
	 * @return A new, default instance of the ADP with the same class as the subclass instance.
	 */
	public ADP makeNew();
	
	public SampleBasedEstimator newEstimator(double[] defaultBandwidth);
	
	/**
	 * Cost at time t, when in state x, applying control u.
	 * @param t
	 * @param x
	 * @param u
	 * @return
	 */
	public double g(int t, RealVector x, RealVector u);

	public void addObservation(int t, RealVector[] x_u, double J);

	public boolean trained(int k, int t);

	public double explorationRate(int k, int t);

	/**
	 * Convenience method for combining a random control with an optimal control.
	 * They will be combined according to the value of {@link #explorationRate(int, int)}.
	 * @param k
	 * @param t
	 * @param uOptimal
	 * @param uRandom
	 * @return
	 */
	public RealVector combineControls(int k, int t, RealVector uOptimal, RealVector uRandom);

	/**
	 * Calculates the costs-to-go for the given control path.
	 * @param u The control path to apply.
	 * @return The costs-to-go.
	 */
	public double[] ctgs(RealVector[] u);

	/**
	 * The state at time 0;
	 * @return
	 */
	public RealVector getX0();
	
	public void setX0(RealVector x_0);

	public void setNr(int i);

	public void setCostFunction(CostFunction g);

	public CostFunction getCostFunction();

	public int getNr();

	public double getGamma();

	public String getName();

	public void setname(String name);

	public void setT(int t);
}