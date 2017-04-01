package ellipsis.energy.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;
import com.mls.util.Timer;

import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.energy.util.GridlessTransformerRegulator;
import ellipsis.test.TestHelper;
import ellipsis.util.Pair;
import ellipsis.util.VectorHelper;

public class GridTestHelper extends TestHelper
{
    protected static final double BASE_POWER = 100e6;
    protected static final double BASE_VOLTAGE = 1e3;
    protected static final double BASE_IMPEDANCE = BASE_VOLTAGE*BASE_VOLTAGE/BASE_POWER;

    /**
     * Uses the power vector for time t=0 as a template and multiplies it by the relative schedule
     * to produce a schedule.
     * @param power_0
     * @param offset
     * @param count
     * @param T
     * @param relativeSchedule
     * @return An array of scheduled powers; time against unit.
     */
	public static double[][] makeSchedule(RealVector power_0, int offset, int count, int T, double[] relativeSchedule)
	{
		double[][] schedule = new double[T][count]; // index 0 is P, index 1 is Q
		for(int i = 0; i < count; ++i)
		{
			double power_0_i = power_0.getEntry(offset+i);
			for(int t = 0; t < T; ++t)
			{
				schedule[t][i] = power_0_i*relativeSchedule[t];
			}
		}
		
		return schedule;
	}

	protected List<Pair<RealVector, Double>> dp(DiscreteGridlessADP adp)
	{
		long start = System.currentTimeMillis();
		
		@SuppressWarnings("unchecked")
		HashMap<RealVector, List<Pair<RealVector, Double>>>[] cache = new HashMap[adp.T+1]; // [0 T]
		List<Pair<RealVector, Double>> ctg = ctg(adp, cache, 0, adp.getX0());
		
		long time = System.currentTimeMillis()-start;
		if(verbose)
			System.out.println("dp() ran in "+time+"ms");
		
		return ctg;
	}
	
	/**
	 * Fills in state, control and noise vectors based on the given schedule.
	 * @param adp Used for x_0 and calculating noise terms.
	 * @param schedule Control and CTG pairs.
	 * @param x State vector to be filled in.
	 * @param u Control vector to be filled in.
	 * @param w Noise vector to be filled in.
	 */
	protected void fillStateControlNoise(GridlessADP adp, List<Pair<RealVector, Double>> schedule, RealVector x[], RealVector u[], RealVector w[])
	{
		int T = schedule.size()-1; // schedule includes a null control element for t=T, so T=size-1
		x[0] = adp.getX0();
		for(int t = 0; t < T; ++t)
		{
			Pair<RealVector, Double> pair = schedule.get(t);
			u[t] = pair.getKey();
			w[t] = adp.wZero(t);
			x[t+1] = adp.f(t, x[t], u[t], w[t]);
		}
	}
	
	/**
	 * Fills in state and noise vectors based on the given control vector.
	 * @param u The controls that define the schedule.
	 * @param adp Used for x_0 and calculating noise terms.
	 * @param x State vector to be filled in.
	 * @param w Noise vector to be filled in.
	 */
	protected void fillStateAndNoiseFromControl(RealVector u[], GridlessADP adp, RealVector x[], RealVector w[])
	{
		int T = u.length;
		x[0] = adp.getX0();
		for(int t = 0; t < T; ++t)
		{
			w[t] = adp.wZero(t);
			x[t+1] = adp.f(t, x[t], u[t], w[t]);
		}
	}

//static int maxChacheSize = 0;
	private List<Pair<RealVector, Double>> ctg(
			DiscreteGridlessADP adp, 
			HashMap<RealVector, List<Pair<RealVector, Double>>>[] cache, 
			int t, 
			RealVector x_t)
	{
		// Check cache:
		if(cache[t] == null)
		{
			cache[t] = new HashMap<>();//(int)Math.pow(t+2, x_t.getDimension())); // Initial capacity set to number of possible states.
		}
		else if(cache[t].get(x_t) != null)
		{
			return cache[t].get(x_t);
		}
		
		// Terminating cost:
		if(t == adp.T)
		{
			// Terminating cost is g(x_t, .):
			double g = adp.g(t, x_t, null);
			
			// Terminating schedule is (., g(.)):
			List<Pair<RealVector, Double>> schedule = new ArrayList<>();
			schedule.add(new Pair<RealVector, Double>(null, g));
			
			// Cache it:
			cache[t].put(x_t, schedule);
//if(maxChacheSize < cache[t].size())
//{
//	maxChacheSize = cache[t].size();
//	System.out.println("Cache size: "+maxChacheSize);
//}
			
			// Return terminating schedule:
			return schedule;
		}
		
		// Find min u \in U:
		RealVector w_noNoise = adp.wZero(t);
		Set<RealVector> U = adp.enumerateControls(t, x_t);
		
		List<Pair<RealVector, Double>> schedule_next_min = null;
		RealVector u_min = null;
		double ctg_min = Double.MAX_VALUE;

		for(RealVector u : U)
		{
			// If constraints are breached then ignore this solution:
			RealVector v = adp.voltagesFromControlAndNoise(t, u, w_noNoise);
			RealVector v_abs = v.getSubVector(v.getDimension()/2, v.getDimension()/2);
			if(v_abs.getMaxValue() > 1.05 || v_abs.getMinValue() < 0.95)
				continue;
			
			//
			RealVector x_next = adp.f(t, x_t, u, w_noNoise);
			double g = adp.g(t, x_t, u);

			List<Pair<RealVector,Double>> schedule_next = ctg(adp, cache, t+1, x_next); // t+1 scheduled controls and ctgs
			double ctg_next = schedule_next.get(0).getValue();
			double ctg = g + ctg_next;
			if(ctg < ctg_min)
			{
				ctg_min = ctg;
				u_min = u;
				schedule_next_min = schedule_next;
			}
		}
		
		// Set schedule value for t:
		List<Pair<RealVector, Double>> schedule = new ArrayList<>();
		schedule.add(new Pair<RealVector, Double>(u_min, ctg_min));
		schedule.addAll(schedule_next_min);
		
		// Add cache value for this state at this time:
		cache[t].put(x_t, schedule);
//if(maxChacheSize < cache[t].size())
//{
//	maxChacheSize = cache[t].size();
//	System.out.println("Cache size: "+maxChacheSize);
//}
		
		return schedule;
	}

	/**
	 * Randomly chooses a control path and builds a CTG array.
	 * Warning: Assumes {@link DiscreteGridlessADP} by default.
	 * @param adp
	 * @param ctgs The CTGs according to the random path (filled in by {@link #randomPath(DiscreteGridlessADP, double[])}).
	 * @return The random control path.
	 */
	public RealVector[] randomPath(GridlessADP gadp, double[] ctgs, Random rand)
	{
		DiscreteGridlessADP adp = (DiscreteGridlessADP)gadp;
		
		RealVector[] u = new RealVector[adp.T];
		RealVector[] x = new RealVector[adp.T+1];
		x[0] = adp.getX0();
		for(int t = 0; t < adp.T; ++t)
 		{
			Set<RealVector> U = adp.enumerateControls(t, x[t]);
			
			int r = /*GridlessADP.*/rand.nextInt(U.size());
			int i = 0;
			for (RealVector _u : U)
			{
				if(i == r)
				{
					u[t] = _u;
					break;
				}
				++i;
			}
			
			RealVector zero = adp.wZero(t);
			x[t+1] = adp.f(t, x[t], u[t], zero);
		}
		
		// Fill in actual costs to go:
		if(ctgs != null)
		{
			double ctg = adp.g(adp.getT(), x[adp.getT()], null); // Terminating CTG.
			ctgs[adp.getT()] = ctg;
			for(int t = adp.getT()-1; t >= 0; --t)
			{
				double g = adp.g(t, x[t], u[t]);
				ctg += g;
				ctgs[t] = ctg;
			}
		}
		
		return u;
	}

	protected double percentageError(
			List<Pair<RealVector, Double>> schedule,
			double[] approximateCtgs, 
			int t)
	{
		Double ctg = schedule.get(t).getValue();
		return Math.abs(ctg - approximateCtgs[t])/ctg;
	}

	protected void printScheduleComparisonAndVoltages(
			DiscreteGridlessADP adp,
			List<Pair<RealVector, Double>> schedule, 
			double[] approximateCtgs,
			RealVector[] approximateSchedule)
	{
		if(verbose)
		{
			System.out.println("\ncorrect,,,approximate");
			for(int t = 0; t < adp.T; ++t)
			{
				System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
			}

			// Voltages:
			printVoltageProfiles(adp, schedule);
		}
	}

	protected void printVoltageProfiles(DiscreteGridlessADP adp, List<Pair<RealVector, Double>> schedule)
	{
		System.out.println("\nVoltages:");
		RealVector x_t = adp.getX0();
		for(int t = 0; t < adp.T; ++t)
		{
			RealVector u_t = schedule.get(t).getKey();
			RealVector wZero = adp.wZero(t);
			
			// Print voltages:
			System.out.print(t+",");
			RealVector v_t = adp.voltagesFromControlAndNoise(t, u_t, wZero);
			for(int i = 0; i < v_t.getDimension(); ++i)
			{
				System.out.print(v_t.getEntry(i)+",");
			}
			System.out.println();

			// Next state:
			x_t = adp.f(t, x_t, u_t, wZero);
		}
	}
	
	protected void printAverageRandomPath(DiscreteGridlessADP adp, double[] randomCtgSum)
	{
		Random rand = new Random(0);
		for(int i = 0; i < iterations; ++i)
		{
			double[] randomCtg = new double[adp.T+1];
			randomPath(adp, randomCtg, rand);
			for (int t = 0; t < randomCtg.length; t++)
			{
				randomCtgSum[t] += randomCtg[t];
			}
		}
		if(verbose)
			System.out.println("Average random path cost-to-go:");
		for (int t = 0; t < randomCtgSum.length; t++)
		{
			randomCtgSum[t] /= iterations;
			if(verbose)
				System.out.println(t+","+randomCtgSum[t]);
		}
	}
	
	
	//// Test cases and logging ////

	protected List<Pair<RealVector, Double>> optimalSchedule;
	protected double[][] iterativeADPCTGs;
	protected double[][] iterativeOLTCCTGs;
	protected RealVector[][] iterativeADPControls;
	protected RealVector[][][] iterativeADPExternalVoltages; // [central iteration][decomposition][time]
	protected RealVector[][] iterativeADPV_0; // [central iteration][decomposition]
	protected RealVector[][] iterativeADPCentralVoltages; // [central iterations][time]
	protected double[][] iterativeOLTCTapPosition; // [central iterations][time]

	protected static RealVector[] randomPath_averageX;
	protected static RealVector[] randomPath_minX;
	protected static RealVector[] randomPath_maxX;

	public double[] randomCase(final GridlessADP adp)
	{
		// Setup min and max:
		int T = adp.getT();
		int dimension = adp.getX0().getDimension();
		randomPath_minX = new RealVector[T];
		randomPath_maxX = new RealVector[T];
		randomPath_averageX = new RealVector[T];
		for(int t = 0; t < T; ++t)
		{
			randomPath_minX[t] = new ArrayRealVector(dimension, Double.MAX_VALUE);
			randomPath_maxX[t] = new ArrayRealVector(dimension, -Double.MAX_VALUE);
			randomPath_averageX[t] = new ArrayRealVector(dimension, 0);
		}
		
		// Run test in another thread:
		final double[] averageCTGs = new double[adp.getT()+1];
		ThreadPool.getInstance().queueTask(new Runnable(){
			public void run() 
			{
				System.out.println("Random comparison started.");
				Random rand = new Random(0);
				int randomIterations = 1000;
				for(int i = 0; i < randomIterations; ++i)
				{
					double[] randomCTGs = new double[adp.getT()+1];
					randomPath(adp, randomCTGs, rand);
					for (int t = 0; t < randomCTGs.length; t++)
					{
						averageCTGs[t] += randomCTGs[t];
					}
				}
				for (int t = 0; t < averageCTGs.length; t++)
				{
					averageCTGs[t] /= randomIterations;
				}
				for (int t = 0; t < adp.getT(); t++)
				{
					randomPath_averageX[t] = randomPath_averageX[t].mapDivide(randomIterations);
				}
				System.out.println("Random comparison complete.");
			};
		});
		return averageCTGs;
	}

	public void adpCase(final GridlessADP adp, final double[][] bw)
	{
		ThreadPool.getInstance().queueTask(new Runnable(){
			public void run() 
			{
				System.out.println("Full network ADP comparison started.");
				Timer.getGlobalTimer("adpCase").start();
				adp.train(bw);
				Timer.getGlobalTimer("adpCase").stop();
				System.out.println("Full network ADP comparison complete.");
			}
		});
	}

	/**
	 * Runs a DP using {@link DiscreteGridlessADP#wZero(int)} at each step and stores
	 * the restult in {@link #optimalSchedule}.
	 * @param adp
	 * @param dpScheduleFileName
	 */
	public void dpCase(final DiscreteGridlessADP adp, final String dpScheduleFileName)
	{
		ThreadPool.getInstance().queueTask(new Runnable()
		{
			@SuppressWarnings({ "unchecked", "rawtypes" })
			public void run() 
			{
				System.out.println("DP comparison started.");
				
				optimalSchedule = (List<Pair<RealVector, Double>>) convertFromSerializable((List<Pair>) read(dpScheduleFileName));
				if(optimalSchedule == null)
					optimalSchedule = dp(adp);
				save(convertToSerializable(optimalSchedule), dpScheduleFileName);
				
				System.out.println("DP comparison complete.");
			}
		});
	}

	public void logCTGs(int T, double[] averageCTGs, int iterations, double[][] iterativeADPCTGs, double[] adpCTGs)
	{
		logCTGs(T, averageCTGs, iterations, iterativeADPCTGs, null, adpCTGs);
	}
	
	public void logCTGs(int T, double[] averageCTGs, int iterations, double[][] iterativeADPCTGs, double[][] iterativeOLTCCTGs, double[] adpCTGs)
	{
		System.out.println("Time,Rand,DP,Full ADP,Semi-dist...");
		for(int t = 0; t <= T; ++t)
		{
			System.out.print(t);
			System.out.print(",");
			
			// Random/average:
			if(averageCTGs != null && averageCTGs.length >= t)
				System.out.print(averageCTGs[t]);
			System.out.print(",");
			
			// DP:
			if(optimalSchedule != null && optimalSchedule.size() >= t)
				System.out.print(optimalSchedule.get(t).getValue());
			System.out.print(",");
			
			// Full network ADP:
			if(adpCTGs != null && adpCTGs.length >= t)
				System.out.print(adpCTGs[t]);
			System.out.print(",");
			
			// Semi-distributed ADP:
			if(	   iterativeADPCTGs != null 
				&& iterativeADPCTGs.length >= iterations 
				&& iterativeADPCTGs[0] != null 
				&& iterativeADPCTGs[0].length >= t)
			{
				for(int i = 0; i < iterations; ++i)
				{
					System.out.print(iterativeADPCTGs[i][t]);
					System.out.print(",");
				}
			}
			
			System.out.println();
		}
		
		// OLTC CTG:
		System.out.println("Time,OLTC CTG...");
		for(int t = 0; t <= T; ++t)
		{
			System.out.print(t);
			System.out.print(",");
			
			// Semi-distributed ADP:
			if(	   iterativeOLTCCTGs != null 
				&& iterativeOLTCCTGs.length >= iterations 
				&& iterativeOLTCCTGs[0] != null 
				&& iterativeOLTCCTGs[0].length >= t)
			{
				for(int i = 0; i < iterations; ++i)
				{
					System.out.print(iterativeOLTCCTGs[i][t]);
					System.out.print(",");
				}
			}
			
			System.out.println();
		}
	}

	/**
	 * Runs the semi-distributed test case and stores the results in {@link #iterativeADPControls},
	 * {@link #iterativeADPCTGs}, and {@link #iterativeADPExternalVoltages}.
	 * @param cc
	 * @param bandwidth
	 * @param iterations
	 */
	public void semiDistributedCase(final GridlessCentralCoordinator cc, final double bandwidth, final int iterations, final double probabilityOfUpdate)
	{
		ThreadPool.getInstance().queueTask(new Runnable(){
			public void run() 
			{
				Timer.getGlobalTimer("semiDistributedCase").start();
				// Log iteration counter to correspond with times logged below:
				{
					StringBuffer sb = new StringBuffer("ADP loop started:,Total Time,");
					for(int k = 0; k < iterations; ++k)
					{
						sb.append(k);
						sb.append(',');
					}
					System.out.println(sb.toString());
				}
				
				cc.decompose();
				cc.probabilityOfUpdate = new double[cc.decomposition.length];
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					cc.probabilityOfUpdate[i] = probabilityOfUpdate;
				}
				
				final long averageTimesPerLocalOptimisation[] = new long[iterations];
				cc.addLoopFinishedCallback(new Runnable()
				{
					int k = 0; // iteration
					long time = System.currentTimeMillis();
					
					@Override
					public void run()
					{
						// Gather data for logging:
						iterativeADPControls[k] = cc.controlSchedule();
						iterativeADPCTGs[k] = cc.gridlessADP().ctgs(iterativeADPControls[k]);
	
						iterativeADPExternalVoltages[k] = new RealVector[cc.decomposition.length][];
						iterativeADPV_0[k] = new RealVector[cc.decomposition.length];
						for (Integer d : cc.gridlessADP().getGridlessData().controlNames.keySet())
						{
							iterativeADPExternalVoltages[k][d] = cc.decomposition[d].getDeltaVExternal();
							iterativeADPV_0[k][d] = cc.decomposition[d].getV0();
						}
						
						long now = System.currentTimeMillis();
						averageTimesPerLocalOptimisation[k] = (now-time)/cc.decomposition.length;
						time = now;
						
						// OLTC Tap Position and CTG:
						int T = cc.localControllerTemplate.getT();
						if(cc.centralTemplate() != null && cc.centralTemplate().getTransformerCount() > 0)
						{
							// Tap:
							int transformerStartIndex = cc.centralTemplate().transformerPowerOffset();
							GridlessTransformerRegulator transformer = (GridlessTransformerRegulator)cc.decomposition[transformerStartIndex]; // First OLTC only.
							iterativeOLTCTapPosition[k] = transformer.copyRatioSchedule();
							
							// CTG:
							iterativeOLTCCTGs[k] = new double[T+1];
							iterativeOLTCCTGs[k][T] = 0;
							for(int t = T-1; t >= 0; --t)
								iterativeOLTCCTGs[k][t] = transformer.cost(t, iterativeOLTCTapPosition[k][t]) + iterativeOLTCCTGs[k][t+1];
						}
						
						// Voltages:
						for(int t = 0; t < T; ++t)
						{
							RealVector v = ((GridlessADP)cc.localControllerTemplate).voltagesFromControlAndNoise(t, iterativeADPControls[k][t], cc.gridlessADP().wZero(t));
							int dim = v.getDimension();
							iterativeADPCentralVoltages[k][t] = v.getSubVector(dim/2, dim/2); // Get magnitude only.
						}
						
						++k;
					}
				});
				
				long startLoop = System.currentTimeMillis();
				cc.loop(iterations, bandwidth);
				long loopTime = System.currentTimeMillis() - startLoop;
				
				// Log average time per decomposition for each iteration:
				{
					StringBuffer sb = new StringBuffer("ADP loop complete:,");
					sb.append(loopTime);
					sb.append(',');
					for(int k = 0; k < iterations; ++k)
					{
						sb.append(averageTimesPerLocalOptimisation[k]);
						sb.append(',');
					}
					System.out.println(sb.toString());
				}
				
				Timer.getGlobalTimer("semiDistributedCase").stop();
			}
		});
	}
}
