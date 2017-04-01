package ellipsis.energy.casestudy;

import static ellipsis.energy.casestudy.TestCaseLogger.logCCConfig;
import static ellipsis.energy.casestudy.TestCaseLogger.logControlStateNoiseVoltage;
import static ellipsis.energy.casestudy.TestCaseLogger.logDataWithNamesAndTimes;
import static ellipsis.energy.casestudy.TestCaseLogger.out;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeDG;

import java.lang.reflect.Method;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.mls.util.ThreadPool;
import com.mls.util.Timer;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.Load;
import ellipsis.energy.test.GridTestHelper;
import ellipsis.energy.test.IEEE13BusGrid;
import ellipsis.energy.util.ADP;
import ellipsis.energy.util.CentralTemplateGridlessADP;
import ellipsis.energy.util.CostFunction;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.energy.util.GridlessData;
import ellipsis.energy.util.LCGridlessADP;
import ellipsis.energy.util.LinearConvexADP;
import ellipsis.energy.util.StateDependent;
import ellipsis.energy.util.TimeDependent;
import ellipsis.energy.util.constraint.ApparentPowerConstraint;
import ellipsis.energy.util.constraint.DeviceConstraint;
import ellipsis.energy.util.constraint.VoltageConstraint;
import ellipsis.util.BarrierFunction;
import ellipsis.util.GroupBarrierFunction;
import ellipsis.util.LogBarrierFunction;
import ellipsis.util.MetricsAspect;
import ellipsis.util.QuadraticEstimator;
import ellipsis.util.SampleBasedEstimator;
import ellipsis.util.TwiceDifferentiableFunction;
import ellipsis.util.VectorHelper;

public class CS05_PQVoltageReg extends GridTestHelper
{
	public static final boolean groupAveragePowerOutput = false;

	public static class PQVoltageRegADP extends CentralTemplateGridlessADP
	{
		private static final double alpha = 1.0/12;
		
//		public PQVoltageRegADP(int dimension)
//		{
//			
//		}
		
		public void setDimension(int dimension)
		{
			costFunction = new PQVoltageRegCost(this, dimension);
		}
		
		/**
		 *   x_t = [ \hat{P}_{DG,t} ]    u_t = [ P_{DG,t} ]    w_t = [ \Delta P_{DG,t} ]
			                                                         [ P_{L,t}         ]
			                                                         [ Q_{L,t}         ]

			 x_{t+1} = f(x_t, u_t, w_t) = Ax_t + Bu_t + Cw_t
			         =>   A = I(1-alpha)
			              B = diag([alpha/\bar{P_{DG,i}} | i=1...])
			              C = B U 0 U 0 // where U is union (by column)
			 where \hat{P}_{DG,t} is the estimated average DG output scaled by DG rating,
			 and \bar{P_{DG,i}} is the DG rating (e.g. 3.5kW).

		 * @param x
		 * @param u
		 * @param w
		 * @return
		 */
		@Override
		public RealVector f(int t, RealVector x_t, RealVector u_t, RealVector w_t)
		{
			// Lazily initialise A, B, C:
			initABC(x_t);
			
			return super.f(t, x_t, u_t, w_t);
		}
		
		private synchronized void initABC(RealVector x_t)
		{
			if(A == null)
			{
				int dimension = x_t.getDimension();
				A = MatrixUtils.createRealIdentityMatrix(dimension).scalarMultiply(1-alpha);
				B = new Array2DRowRealMatrix(dimension, controlDimension());
				C = new Array2DRowRealMatrix(dimension, gridlessData.dgCount+2*gridlessData.loadCount); // DG noise and loads (P&Q)
				for(int i = 0; i < gridlessData.dgCount; ++i)
				{
					double b_ii = alpha/gridlessData.rating_DG[i];
					int j;
					if(groupAveragePowerOutput)
					{
						b_ii /= gridlessData.dgCount;
						j = 0;
					}
					else
					{
						j =  i;
					}
					B.setEntry(j, i, b_ii);
					C.setEntry(j, i, b_ii);
				}
			}
		}

		@Override
		public RealVector getX0()
		{
			if(groupAveragePowerOutput)
				return VectorHelper.vector(1.0);
			else
				return super.getX0();
		}
		
		@Override
		public ADP makeNew()
		{
			PQVoltageRegADP reg = new PQVoltageRegADP();
			reg.setDimension(0);
			return reg;
		}
		
		@Override
		public SampleBasedEstimator newEstimator(double[] defaultBandwidth)
		{
			QuadraticEstimator qe = new QuadraticEstimator();
			qe.setMaxSampleCount(60);
			return qe;
		}
		
		/*@Override
		public double explorationRate(int k, int t)
		{
			if(beforeTraining(k, t) || approxJ[t].getSamplePoints().size() < 5)
				return 1; // all random before training
			else
				return 0;
//			if(trained(k, t))
//				return 0; // only optimal after training
//			
//			double rate = T-t-k/(double)N_r; // linearly from random to optimal during training
//
//			if(rate < 0 || rate > 1)
//				throw new RuntimeException("Invalid rate "+rate);
//			return rate;
		}*/
	}
	
	public static class PQVoltageRegCost implements CostFunction, TwiceDifferentiableFunction, TimeDependent, StateDependent, ConvexMultivariateRealFunction
	{
		private int t;
		private ADP adp;
		private RealVector x_t;
		int dim;
		
		public PQVoltageRegCost(ADP adp, int dimension)
		{
			this.adp = adp;
			this.dim = dimension;
		}

		@Override
		public double g(ADP adp, int t, RealVector x, RealVector u)
		{
			RealVector one = VectorHelper.one(x.getDimension());
			RealVector oneMinusX = one.subtract(x);
			return oneMinusX.dotProduct(oneMinusX); // (1-x)^2
		}

		@Override
		public RealVector gradient(RealVector u)
		{
			// g(x,u) is constant with respect to u:
			return null; // null means [0]
		}

		@Override
		public RealMatrix hessian(RealVector u)
		{
			// g(x,u) is constant with respect to u:
			return null; // null means [0]
		}

		@Override
		public double value(RealVector u)
		{
			return g(adp, t, x_t, u);
		}
		
		@Override
		public void setTime(int t)
		{
			this.t = t;
		}

		@Override
		public void setState(RealVector x)
		{
			this.x_t = x;
		}
		
		
		//// JOptimizer ////

		@Override
		public double value(double[] u)
		{
			return g(adp, t, x_t, new ArrayRealVector(u));
		}

		@Override
		public double[] gradient(double[] u)
		{
			return new double[dim];
		}

		@Override
		public double[][] hessian(double[] u)
		{
			return new double[dim][dim];
		}

		@Override
		public int getDim()
		{
			return dim;
		}
	}
	
	public static class PQVoltageRegCentralController extends GridlessCentralCoordinator
	{
		@Override
		public void fillInForecastAndInitialState(GridlessADP adp_B, List<Integer> sortedB)
		{
			super.fillInForecastAndInitialState(adp_B, sortedB);

			GridlessData gridlessData = ((PQVoltageRegADP)localControllerTemplate).getGridlessData();
			GridlessData gridlessData_B = ((PQVoltageRegADP)adp_B).getGridlessData();
			int dgCount = gridlessData.dgCount;
			int dgCount_B = gridlessData_B.dgCount;
			int loadCount_B = gridlessData_B.loadCount;
			int loadCount = gridlessData.loadCount;
			int controlDimension = gridlessADP().controlDimension();
			int controlDimension_B = adp_B.controlDimension();
			
			// Set cost function with correct dimension:
			adp_B.setCostFunction(new PQVoltageRegCost((ADP)adp_B, adp_B.controlDimension()));
			
			// Set x_0s:
			ArrayRealVector x_0 = new ArrayRealVector(adp_B.getGridlessData().dgCount, 1.0);
			adp_B.setX0(x_0);
			
			// Set DG ratings and zeroControl:
			int i_B = 0;
			gridlessData_B.rating_DG = new double[adp_B.getGridlessData().dgCount];
			double[] zeroControl = ((PQVoltageRegADP)localControllerTemplate).getZeroControl();
			double[] zeroControl_B = new double[controlDimension_B];
			for (Integer i : sortedB)
			{
				if(i >= dgCount)
					break;
				
				gridlessData_B.rating_DG[i_B] = gridlessData.rating_DG[i];
				zeroControl_B[i_B] = zeroControl[i];
				zeroControl_B[i_B+controlDimension_B/2] = zeroControl[i+controlDimension/2];
				
				++i_B;
			}
			((PQVoltageRegADP)adp_B).setZeroControl(zeroControl_B);
			
			// Set noise mean and SD:
			// w = [Delta P_DG  P_L  Q_L]
			double[][] noiseMean = ((LinearConvexADP)localControllerTemplate).getNoiseMean();
			double[][] noiseSD = ((LinearConvexADP)localControllerTemplate).getNoiseSD();
			int T = localControllerTemplate.getT();
			double[][] noiseMean_B = new double[T][dgCount_B+2*loadCount_B];
			double[][] noiseSD_B = new double[T][dgCount_B+2*loadCount_B];
			i_B = 0;
			for (Integer i : sortedB)
			{
				// Delta P_DG = 0 for zero mean.
				// P_L & Q_L:
				if(i >= gridlessADP().getGridlessData().dgCount)
				{
					for(int t = 0; t < T; ++t)
					{
						noiseMean_B[t][i_B] = noiseMean[t][i]; // P
						noiseSD_B[t][i_B] = noiseSD[t][i];
						noiseMean_B[t][i_B+loadCount_B] = noiseMean[t][i+loadCount]; // Q
						noiseSD_B[t][i_B+loadCount_B] = noiseSD[t][i+loadCount];
					}
				}
				
				++i_B;
			}
			((LinearConvexADP)adp_B).setNoiseMean(noiseMean_B);
			((LinearConvexADP)adp_B).setNoiseSD(noiseSD_B);
			
			// Copy constraints:
			GroupBarrierFunction group = (GroupBarrierFunction) ((LinearConvexADP)gridlessADP()).getBarrier();
			List<BarrierFunction> barriers = group.getBarriers();
			GroupBarrierFunction group_B = new GroupBarrierFunction();
			i_B = 0;
			for (Integer i : sortedB)
			{
				group_B.addBarrier(copy(barriers.get(i), adp_B, i_B));
				++i_B;
			}
			((LinearConvexADP)adp_B).setBarrier(group_B);
		}
		
		private BarrierFunction copy(BarrierFunction barrierFunction, GridlessADP adp_B, int i_B)
		{
			LogBarrierFunction newBarrier = new LogBarrierFunction();
			int dim = adp_B.controlDimension();
			for (TwiceDifferentiableFunction constraint : barrierFunction.getConstraints())
			{
				if(constraint instanceof DeviceConstraint)
				{
					DeviceConstraint dc = (DeviceConstraint)constraint;
					DeviceConstraint newDC = new DeviceConstraint(i_B, dc.getMin(), dc.getMax(), dim);
					newBarrier.addConstraint(newDC);
					
				}
				else if(constraint instanceof ApparentPowerConstraint)
				{
					ApparentPowerConstraint apc = (ApparentPowerConstraint)constraint;
					ApparentPowerConstraint newAPC = new ApparentPowerConstraint(i_B, apc.getSSquared(), dim);
					newBarrier.addConstraint(newAPC);
				}
				else if(constraint instanceof VoltageConstraint)
				{
					VoltageConstraint vc = (VoltageConstraint)constraint;
					VoltageConstraint newVC = new VoltageConstraint(adp_B, dim);
					newVC.setA(vc.getA());
					newVC.setB(vc.getB());
					newVC.setI(i_B);
					newBarrier.addConstraint(newVC);
				}
			}
			return newBarrier;
		}

		@Override
		public double[][] makeBandwidthArray(int storageCount, int dgCount, int loadCount, double socBW, double dgBW, double storageBW, double loadBW)
		{
			double[][] bw = new double[localControllerTemplate.getT()][dgCount];
			for(int t = 0; t < localControllerTemplate.getT(); ++t)
				for(int i = 0; i < dgCount; ++i)
				{
					bw[t][i] = (t+1)*dgBW; // t = 0 -> T-1 gives bw = gdBW -> T*dgBW
				}
			return bw;
		}
		
		@Override
		public void makeGridlessADP(List<List<Integer>> Bs)
		{
			super.makeGridlessADP(Bs);
			
			for(int b = 0; b < decomposition.length; ++b)
			{
				if(decomposition[b] instanceof LCGridlessADP)
					((LCGridlessADP)decomposition[b]).finaliseData();
			}
		}
	}
	
	public static void main(String[] args)
	{
		assertAssertsOn();
		CS05_PQVoltageReg cs05 = new CS05_PQVoltageReg();
		cs05.parseArgs(args);
		cs05.run();
	}
	
	public void run()
	{
		final GridlessCentralCoordinator cc = makeCentralController();
		int iterations = 20;
		double probabilityOfUpdate = 1;
		double h = 2.0/12;
		
		logCCConfig(cc, iterations, probabilityOfUpdate, h);
		
		// Random comparison:
		double[] averageCTGs = randomCase(cc.gridlessADP());

//		System.out.println("-");
//		System.out.println("Random path range:");
//		System.out.print("t,");
//		for(int i = 0; i < randomPath_minX[0].getDimension(); ++i)
//		{
//			System.out.print("x["+i+"]av,min,max,");
//		}
//		System.out.println();
//		for(int t = 0; t < randomPath_minX.length; ++t)
//		{
//			System.out.print(t);
//			System.out.print(',');
//			for(int i = 0; i < randomPath_minX[0].getDimension(); ++i)
//			{
//				System.out.print(randomPath_averageX[t].getEntry(i));
//				System.out.print(',');
//				System.out.print(randomPath_minX[t].getEntry(i));
//				System.out.print(',');
//				System.out.print(randomPath_maxX[t].getEntry(i));
//				System.out.print(',');
//			}
//			System.out.println();
//		}

		// Full network ADP results:
		boolean disableADP = false;
		GridlessADP adp = makeGrid(false); // no transformer
		if(!disableADP)
		{
			double[][] bw = cc.makeBandwidthArray(0, adp.getGridlessData().dgCount, 0, 0, h, 0, 0);
			adpCase(adp, bw);
		}
		else
		{
			out.println("Full network ADP comparison not run.\n");
		}
		
		// Semi-distributed ADP results:
		iterativeADPCTGs = new double[iterations][]; // [loop iteration][time]
		iterativeOLTCCTGs = new double[iterations][]; // [loop iteration][time]
		iterativeADPControls = new RealVector[iterations][]; // [loop iteration][time]
		iterativeADPExternalVoltages = new RealVector[iterations][][]; // [loop iteration][decomposition][time]
		iterativeADPV_0 = new RealVector[iterations][]; // [loop iteration][decomposition]
		iterativeADPCentralVoltages = new RealVector[iterations][cc.localControllerTemplate.getT()];
		iterativeOLTCTapPosition = new double[iterations][];
		semiDistributedCase(cc, h, iterations, probabilityOfUpdate);

		ThreadPool.getInstance().waitForAll();
		
System.out.println("\n\n"+showResults());
		
		// Get full ADP schedule:
		double[] adpCTGs = new double[adp.getT()+1];
		RealVector[] adpStates = disableADP ? null : new RealVector[adp.getT()+1];
		RealVector[] adpNoise = disableADP ? null : new RealVector[adp.getT()+1];
		RealVector[] adpControls = disableADP ? null : adp.schedule(adpCTGs, adpStates, adpNoise);
		
		
		//// Logging: ////
		
		log(cc, iterations, averageCTGs, adp, adpCTGs, adpStates, adpNoise, adpControls);
		
		System.exit(0);
	}

	private void log(final GridlessCentralCoordinator cc, int iterations,
			double[] averageCTGs, GridlessADP adp, double[] adpCTGs,
			RealVector[] adpStates, RealVector[] adpNoise,
			RealVector[] adpControls)
	{
		out.println("Converged=,"+LinearConvexADP.converged+
				", non-converged=,"+LinearConvexADP.notConverged);
		
		out.println("Costs-to-go:");
		logCTGs(cc.localControllerTemplate.getT(), averageCTGs, iterations, iterativeADPCTGs, iterativeOLTCCTGs, adpCTGs);
		
		out.println(",Epsilon,Delta v at pilot t=1 k=[0.."+iterations+"]...");
		for(int i = 0; i < cc.decomposition.length; ++i)
		{
			Integer pilotIndex = cc.pilots.get(i);
			Integer pilotIndex_local = cc.pilots_i.get(i);
			String pilotName = cc.gridlessADP().getGridlessData().controlNames.get(pilotIndex);
			if(pilotName == null)
				continue;
			out.print(pilotName);
			out.print(",");
			out.print(cc.epsilon(i));
			out.print(",");
			
			for(int k = 0; k < iterativeADPExternalVoltages.length; ++k)
			{
				RealVector[] extDeltaV = iterativeADPExternalVoltages[k][i];
				out.print(extDeltaV[1].getEntry(pilotIndex_local)); // Pilot bus at time 1.
				out.print(",");
			}
			
			out.println();
		}
		
		out.println("-");
		out.println("**** ADP ****");
		out.println(
				"Times (ms):,ADP:,"+Timer.getGlobalTimer("adpCase").getTimeElapsed()
				+",Semi-dist:,"+Timer.getGlobalTimer("semiDistributedCase").getTimeElapsed()
//				+",estimateOptimalControl:,"+Timer.getGlobalTimer("estimateOptimalControl").getTimeElapsed()
//				+",VC.value:,"+Timer.getGlobalTimer("VC.value").getTimeElapsed()
//				+",powerFromControlAndNoise:,"+Timer.getGlobalTimer("powerFromControlAndNoise").getTimeElapsed()
//				+",sensitivities_i:,"+Timer.getGlobalTimer("sensitivities_i").getTimeElapsed()
//				+",v_abs_i:,"+Timer.getGlobalTimer("v_abs_i").getTimeElapsed()
				
//				+",min.f.g:,"+Timer.getGlobalTimer("min.f.g").getTimeElapsed()
//				+",min.b.g:,"+Timer.getGlobalTimer("min.b.g").getTimeElapsed()
//				+",min.f.h:,"+Timer.getGlobalTimer("min.f.h").getTimeElapsed()
//				+",min.b.h:,"+Timer.getGlobalTimer("min.b.h").getTimeElapsed()
//				+",VC.v:,"+Timer.getGlobalTimer("VC.v").getTimeElapsed()
//				+",VC.g:,"+Timer.getGlobalTimer("VC.g").getTimeElapsed()
//				+",VC.h:,"+Timer.getGlobalTimer("VC.h").getTimeElapsed()
//				+",count VC.v:,"+VoltageConstraint.count_v
//				+",count VC.g:,"+VoltageConstraint.count_g
//				+",LBF.h.comp:,"+Timer.getGlobalTimer("LBF.h.comp").getTimeElapsed()
//				+",LBF.g.constraint.v:,"+Timer.getGlobalTimer("LBF.g.constraint.v_").getTimeElapsed()
//				+",LBF.g.constraint.g:,"+Timer.getGlobalTimer("LBF.g.constraint.g").getTimeElapsed()
//				+",LBF.h.constraint.v:,"+Timer.getGlobalTimer("LBF.h.constraint.v_").getTimeElapsed()
//				+",LBF.h.constraint.g:,"+Timer.getGlobalTimer("LBF.h.constraint.g").getTimeElapsed()
//				+",LBF.h.constraint.h:,"+Timer.getGlobalTimer("LBF.h.constraint.h").getTimeElapsed()
//				+",LBF.count_g_constraint_v:,"+LogBarrierFunction.count_g_constraint_v
//				+",LBF.count_h_constraint_v:,"+LogBarrierFunction.count_h_constraint_v
				);

		{
			logControlStateNoiseVoltage(
					adpStates,
					adpControls,
					adpNoise,
					groupAveragePowerOutput ? new HashMap<Integer, String>() : adp.getGridlessData().controlNames,
					adp.getGridlessData().stateNames,
					adp.getGridlessData().noiseNames,
					adp.getGridlessData().voltageNames,
					adp);
		}
		
		out.println("-");
		out.println("**** Semi-distributed ADP ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			out.println("---- Iteration "+k+" ----");
			
			// OLTC:
			out.print("time:,");
			for(int t = 0; t < cc.centralTemplate().getT(); ++t)
			{
				out.print(t);
				out.print(",");
			}
			out.println();
			out.print("OLTC,");
			for(int t = 0; t < cc.centralTemplate().getT(); ++t)
			{
				if(iterativeOLTCTapPosition[k] != null)
					out.print(iterativeOLTCTapPosition[k][t]);
				out.print(",");
			}
			out.println();
			
			// State, control and noise:
			RealVector[] x = new RealVector[cc.gridlessADP().getT()+1];
			RealVector[] w = new RealVector[cc.gridlessADP().getT()];
			RealVector[] v = iterativeADPCentralVoltages[k];
			fillStateAndNoiseFromControl(iterativeADPControls[k], cc.gridlessADP(), x, w);
			logControlStateNoiseVoltage(
					x,
					iterativeADPControls[k],
					w,
					v,
					cc.gridlessADP().getGridlessData().controlNames,
					cc.gridlessADP().getGridlessData().stateNames,
					cc.gridlessADP().getGridlessData().noiseNames,
					cc.gridlessADP().getGridlessData().voltageNames,
					cc.gridlessADP());
		}

		// Log external delta voltages:
		out.println("-");
		out.println("**** External Delta V ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			out.println("---- Iteration "+k+" ----");
			for(int d = 0; d < cc.decomposition.length; ++d)
			{
				out.println("---- Decomposition "+d+" ----");
				RealVector[] extDeltaV = iterativeADPExternalVoltages[k][d];
				if(extDeltaV != null && cc.decomposition(d) != null)
					logDataWithNamesAndTimes(cc.gridlessADP().getT(), extDeltaV, cc.decomposition(d).getGridlessData().controlNames);
			}
		}

		// Log v_0:
		out.println("-");
		out.println("**** V_0 ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			out.println("---- Iteration "+k+" ----");
			for(int d = 0; d < cc.decomposition.length; ++d)
			{
				GridlessADP decomposition = cc.decomposition(d);
				if(decomposition == null)
					continue;
				Map<Integer, String> names = decomposition.getGridlessData().controlNames;
				out.println("---- Decomposition "+d+" ----");
				for (Integer i : names.keySet())
				{
					// First column (unit name):
					out.print(names.get(i));
					out.print(',');
					out.println(iterativeADPV_0[k][d].getEntry(i));
				}
			}
		}
	}
	
	/**
	 * Randomly chooses a control path and builds a CTG array. Paths may not be viable.
	 * @param adp
	 * @param ctgs The CTGs according to the random path (filled in by {@link #randomPath(GridlessADP, double[])}).
	 * @return The random control path.
	 */
	@Override
	public RealVector[] randomPath(GridlessADP adp, double[] ctgs, Random rand)
	{
		RealVector[] u = new RealVector[adp.getT()];
		RealVector[] x = new RealVector[adp.getT()+1];
		int controlDimension = adp.controlDimension();
		x[0] = adp.getX0();
		for(int t = 0; t < adp.getT(); ++t)
 		{
			// Remember min, max and avergae state:
			for(int i = 0; i < x[t].getDimension(); ++i)
			{
				double x_t_i = x[t].getEntry(i);
				if(x_t_i < randomPath_minX[t].getEntry(i))
				{
					randomPath_minX[t].setEntry(i, x_t_i);
				}
				else if(x_t_i > randomPath_maxX[t].getEntry(i))
				{
					randomPath_maxX[t].setEntry(i, x_t_i);
				}
			}
			randomPath_averageX[t] = randomPath_averageX[t].add(x[t]);
			
			// Random control:
			u[t] = new ArrayRealVector(controlDimension);
			for(int i = 0; i < adp.getGridlessData().dgCount; ++i) // can't include transformer here
			{
				double totalPower = rand.nextDouble()*adp.getGridlessData().meanP_DG[t][i];
				double p = rand.nextDouble()*totalPower;
				double q = Math.sqrt(totalPower*totalPower - p*p);
				u[t].setEntry(i, p);
				u[t].setEntry(controlDimension/2+i, q);
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

	public GridlessCentralCoordinator makeCentralController()
	{
		GridlessCentralCoordinator cc = new PQVoltageRegCentralController();
		GridlessADP adpTemplate = makeGrid(/*includeTransformer*/true);
		cc.C_max = 2;
		cc.U_max = 4;
		cc.alpha = 0.5; // FIXME alpha
		
		// Set the template (must be called after setting the transformer):
		cc.setADPTemplate(adpTemplate);
		
		return cc;
	}

	Map<Object, Object> results = new HashMap<>();
	
	public String showResults()
	{
		try
		{
			StringBuffer sb = new StringBuffer();
			sb.append("method,time,count\n");
			Method[] methods = PQVoltageRegADP.class.getMethods();
			for (int i = 0; i < methods.length; i++)
			{
				String id = MetricsAspect.id(PQVoltageRegADP.class, methods[i]);
				sb.append(id);
				sb.append(',');
				sb.append(Timer.getGlobalTimer(id).getTimeElapsed());
				sb.append(',');
				sb.append(results.get(id));
				sb.append('\n');
			}
			return sb.toString();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	public GridlessADP makeGrid(boolean includeTransformer)
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();

		int dgCount = 6;
		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
		
		// Remove transformer at bus 650:
		if(!includeTransformer)
		{
			Bus from = grid.xfm650.getFromBus();
			Bus to = grid.xfm650.getToBus();
			Line line650 = new Line();
			line650.setName("L650");
			line650.setFromBus(from);
			line650.setToBus(to);
			line650.setImpedencePerMetre(grid.xfm650.getImpedance());
			line650.setLength(grid.xfm650.getLength());
			from.removeChild(grid.xfm650);
			to.removeChild(grid.xfm650);
			if(grid.getLine("L650") == null || grid.getTransformer(grid.xfm650.getName()) != null)
				throw new RuntimeException("Failed to replace transformer 650 with a line.");
		}
		
		// Remove transformer at bus 633:
		Bus from = grid.xfm633_634.getFromBus();
		Bus to = grid.xfm633_634.getToBus();
		Line line633_634 = new Line();
		line633_634.setName("L633-634");
		line633_634.setFromBus(from);
		line633_634.setToBus(to);
		line633_634.setImpedencePerMetre(grid.xfm633_634.getImpedance());
		line633_634.setLength(grid.xfm633_634.getLength());
		from.removeChild(grid.xfm633_634);
		to.removeChild(grid.xfm633_634);
		if(grid.getLine("L633-634") == null || grid.getTransformer(grid.xfm633_634.getName()) != null)
			throw new RuntimeException("Failed to replace transformer 633-634 with a line.");
		
		reduceLoads(grid);
		
//		GridDisplay.showInFrame(grid, IEEE13BusGrid.BASE_POWER, IEEE13BusGrid.BASE_VOLTAGE)
//		           .setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		int contolDimension = dgCount*2;
//		PQVoltageRegADP reg = new PQVoltageRegADP();
		PQVoltageRegADP reg = MetricsAspect.newInterceptor(PQVoltageRegADP.class, results);
		reg.setDimension(contolDimension);
		CentralTemplateGridlessADP adp = reg;
		adp.setT(12);
		makeGridlessADP(grid, adp);
		
		if(adp.getGridlessData().dgCount != dgCount)
			throw new RuntimeException("dgCount incorrect: "+dgCount+", should be "+adp.getGridlessData().dgCount);

		adp.setNr(30);
		
		// State is accumulated time spent turned off:
		ArrayRealVector x_0 = new ArrayRealVector(groupAveragePowerOutput ? 1 : dgCount, 1.0);
		adp.setX0(x_0); // Default all to 1 at time 0.
		
		// Set DG sizes:
		double[] ratings = new double[dgCount];
		double[][] meanP_DG = adp.getGridlessData().meanP_DG;
		for(int t = 0; t < adp.getT(); ++t)
		{
			for(int i = 0; i < dgCount; ++i)
			{
				if(meanP_DG[t][i] > ratings[i])
					ratings[i] = meanP_DG[t][i];
			}
		}
		((PQVoltageRegADP)adp).getGridlessData().rating_DG = ratings;
		
		// Set state names to be DG:
		int i = 0;
		for (DistributedSource dg : grid.getDistributedSources(false))
		{
			adp.getGridlessData().stateNames.put(i, dg.getName());
			++i;
		}
		
		// Setup noise:
		int loadCount = adp.getGridlessData().loadCount;
		double[][] noiseMean = new double[adp.getT()][dgCount+loadCount*2];
		double[][] noiseSD = new double[adp.getT()][dgCount+loadCount*2];
		for(int t = 0; t < adp.getT(); ++t)
		{
			for(i = 0; i < dgCount; ++i)
			{
				noiseMean[t][i] = 0;
				noiseSD[t][i] = adp.getGridlessData().sdP_DG[t][i];
			}
			for(i = 0; i < loadCount; ++i)
			{
				noiseMean[t][dgCount+i] = adp.getGridlessData().meanP_L[t][i];
				noiseSD[t][dgCount+i] = adp.getGridlessData().sdP_L[t][i];
			}
			for(i = 0; i < loadCount; ++i)
			{
				noiseMean[t][dgCount+loadCount+i] = adp.getGridlessData().meanQ_L[t][i];
				noiseSD[t][dgCount+loadCount+i] = adp.getGridlessData().sdQ_L[t][i];
			}
		}
		((PQVoltageRegADP)adp).setNoiseMean(noiseMean);
		((PQVoltageRegADP)adp).setNoiseSD(noiseSD);
		
		// Setup zero controls:
		double[] zeroControl = new double[adp.controlDimension()];
		for (int j = 0; j < zeroControl.length; j++)
		{
			zeroControl[j] = 1e-12; // A little hack since zero is technically not a viable control for the barrier solver.
		}
		((PQVoltageRegADP)adp).setZeroControl(zeroControl);
		
		// Constraints:
		GroupBarrierFunction barrier = new GroupBarrierFunction();
		((PQVoltageRegADP)adp).setBarrier(barrier);
		
		int unitCount = dgCount+loadCount;
		for(i = 0; i < unitCount; ++i)
		{
			LogBarrierFunction barrier_i = new LogBarrierFunction();
			
			// Voltage constraints:
			VoltageConstraint lowerBound = new VoltageConstraint(adp, contolDimension);
			lowerBound.setA(-1.0);
			lowerBound.setB(0.95);
			lowerBound.setI(i);
			barrier_i.addConstraint(lowerBound);

			VoltageConstraint upperBound = new VoltageConstraint(adp, contolDimension);
			upperBound.setA(1.0);
			upperBound.setB(-1.05);
			upperBound.setI(i);
			barrier_i.addConstraint(upperBound);
			
			// DG constraints:
			if(i < dgCount)
			{
				// P > 0:
				DeviceConstraint dc_min = new DeviceConstraint(i, new double[adp.getT()], null, contolDimension); // >= 0
				barrier_i.addConstraint(dc_min);
				
				// P^2 + Q^2 <= S^2:
				ApparentPowerConstraint apc = new ApparentPowerConstraint(i, meanP_DG, contolDimension);
				barrier_i.addConstraint(apc);
			}
			
			// Add constraints for this unit:
			barrier.addBarrier(barrier_i);
		}
		
		((LCGridlessADP)adp).finaliseData();
		
		return adp;
	}

	private void reduceLoads(IEEE13BusGrid grid)
	{
		for(Load load : grid.getLoads())
		{
			load.setLoad(load.getLoad().multiply(0.4));
		}
	}
	
	public static void makeGridlessADP(IEEE13BusGrid grid, CentralTemplateGridlessADP adp)
	{
		GridlessADPFactory.fromGridToGridless(grid, adp, IEEE13BusGrid.BASE_VOLTAGE, IEEE13BusGrid.BASE_POWER, /*failIfNotConvergent=*/true);

		// Setup schedule:
		adp.getGridlessData().setMeanP_DG(new double[] {1.6,1.6,1.55,1.5,1.4,1.25,0.85,0.55,0.2,0.1,0,0});
		adp.getGridlessData().setsdP_DG(new double[] {0.02, 0.02, 0.02, 0.02, 0.015, 0.012, 0.01, 0.005, 0.002, 0.01, 0.0, 0.0});

		// Set load schedule:
		double[] relativeLoadSchedule = new double[]{0.5, 0.7, 0.8, 0.7, 0.9, 1.2, 1.4, 1.2, 0.7, 0.5, 0.4, 0.3};
		
		int pOffset = adp.loadPowerOffset();
		int qOffset = pOffset + adp.unitCount();
		
		adp.getGridlessData().meanP_L = makeSchedule(adp.getGridlessData().power_0, pOffset, adp.getGridlessData().loadCount, adp.getT(), relativeLoadSchedule);
		adp.getGridlessData().meanQ_L = makeSchedule(adp.getGridlessData().power_0, qOffset, adp.getGridlessData().loadCount, adp.getT(), relativeLoadSchedule);
		
		adp.getGridlessData().setsdP_L(new double[] {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
		adp.getGridlessData().setsdQ_L(new double[] {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
		
		adp.setNr(100);
	}
	
	public void showScenario(CentralTemplateGridlessADP adp)
	{
		double[][] pls = adp.getGridlessData().meanP_L;
		double[][] qls = adp.getGridlessData().meanQ_L;
		double[][] pdgs = adp.getGridlessData().meanP_DG;
		for (int t = 0; t < pls.length; t++) 
		{
			System.out.print(t);
			System.out.print('&');
			System.out.print(average(pls[t]));
			System.out.print('&');
			System.out.print(average(qls[t]));
			System.out.print('&');
			System.out.print(average(pdgs[t]));
			System.out.println("\\\\");
		}
	}

	private String average(double[] ds)
	{
		DecimalFormat format = new DecimalFormat("#.##");
		double d = 0;
		for(int i = 0; i < ds.length; ++d)
			d += ds[i];
		return format.format(d/ds.length);
	}
}