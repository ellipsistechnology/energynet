package ellipsis.energy.util;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;

public abstract class CentralCoordinator
{
	static ThreadPool pool;

	// Parameters:
	private Set<Runnable> loopFinishedCallbacks = new HashSet<>();

	public double alpha = 1.0; // Defaults to no dampening.
	public double[] probabilityOfUpdate; // A probability for each decomposition - defaults to 1.
	public LocalController localControllerTemplate;// adpTemplate;
	public LocalController[] decomposition;//localADP;
	
	// Variables:
	protected Random rand_update = new Random(5);
	protected RealVector[] previousControlSchedule;
	public RealVector[][] schedule; // [subset][time]
	public double[][] ctgs; // [subset][time]
	
	
	//// Optimization ////
	
	/**
	 * Combines all subset controls, dampened by {@link #alpha} with respect to 
	 * the previous control schedule, into a control vector for each time.
	 * @return
	 */
	public abstract RealVector[] controlSchedule();
	
	/**
	 * Updates {@link #decomposition}s' values.
	 * @param u_0
	 */
	public abstract void update(RealVector[] u);
	
	protected abstract double[][] makeBandwidthArray(ADPInterface gridlessADP, double bw);

	public void loop(int iterations, double socBW)
	{
		for(int j = 0; j < iterations; ++j)
		{
System.out.print(j);
			iteration(socBW);
			for (Runnable callback : loopFinishedCallbacks)
			{
				callback.run();
			}
		}
	}

	/**
	 * Calls {@link #optimise(double, double, double, double, double, double, double)} then updates
	 * {@link #adpTemplate}.power_0, {@link #adpTemplate}.v_0, {@link #localADP}[*].power_0 and {@link #localADP}[*].v_0.
	 * @param bandwidth
	 * @param previousControl 
	 * @param dgPBW
	 * @param storagePBW
	 * @param loadPBW
	 * @param dgQBW
	 * @param storageQBW
	 * @param loadQBW
	 * @return Implemented control.
	 */
	public RealVector[] iteration(double bandwidth)
	{
		// Optimise:
		optimise(bandwidth);

		// Combine local controls into control vector and update:
		RealVector[] u = controlSchedule();
		update(u);
		
		// Remember control schedule:
		if(previousControlSchedule == null)
			previousControlSchedule = new RealVector[localControllerTemplate.getT()];
		for(int t = 0; t < previousControlSchedule.length; ++t)
		{
			if(u[t] != null)
				previousControlSchedule[t] = u[t];
			else if(previousControlSchedule[t] == null)
				previousControlSchedule[t] = localControllerTemplate.uZero();
		  //else leave old value
		}
		
		return u;
	}

	/**
	 * Trains each decompositions and gets the approximately optimal schedule.
	 * {@link #probabilityOfUpdate} is used to determine if each decomposition
	 * is actually trained and its schedule obtained.
	 * @param bandwidth Bandwidth for estimator.
	 */
	public void optimise(double bandwidth)
	{
		if(pool == null)
		{
			pool = new ThreadPool(Math.min(decomposition.length,4)); // Assuming 4 cores.
		}
		
		schedule = new RealVector[decomposition.length][];
		ctgs = new double[decomposition.length][localControllerTemplate.getT()+1];
		
		for(int i = 0; i < decomposition.length; ++i)
		{
			// Simulate delayed updated info:
			if(rand_update.nextDouble() > probabilityOfUpdate[i])
				continue;
			
			trainLocal(pool, bandwidth, i);
		}
		pool.waitForAll();
	}
	
	private void trainLocal(ThreadPool pool, double bw, final int i)
	{
		final double[][] bandwidth;
		if(decomposition[i] instanceof ADPInterface)
			bandwidth= makeBandwidthArray((ADPInterface) decomposition[i], bw);
		else
			bandwidth = null;
		
		pool.queueTask(new Runnable()
		{
			@Override
			public void run()
			{
				LocalController decomposition_i = decomposition[i];
				decomposition_i.train(bandwidth);
				schedule[i] = decomposition_i.schedule(ctgs[i]);
			}
		});
	}
	
	public boolean addLoopFinishedCallback(Runnable arg0)
	{
		return loopFinishedCallbacks.add(arg0);
	}
}