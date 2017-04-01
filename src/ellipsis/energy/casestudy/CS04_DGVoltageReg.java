package ellipsis.energy.casestudy;

import static ellipsis.energy.casestudy.TestCaseLogger.logCCConfig;
import static ellipsis.energy.casestudy.TestCaseLogger.logControlStateNoiseVoltage;
import static ellipsis.energy.casestudy.TestCaseLogger.logDataWithNamesAndTimes;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeDG;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeGridlessADP;
import static ellipsis.util.VectorHelper.vector;

import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.test.GridTestHelper;
import ellipsis.energy.test.IEEE13BusGrid;
import ellipsis.energy.util.ADP;
import ellipsis.energy.util.CostFunction;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.util.QuadraticEstimator;
import ellipsis.util.SampleBasedEstimator;
import ellipsis.util.VectorHelper;

/**
 * Tests the case of using DG to regulate voltage while aiming to minimise time
 * spent with DG off.
 * @author bmillar
 */
public class CS04_DGVoltageReg extends GridTestHelper
{
	public static class DGVoltageRegADP extends DiscreteGridlessADP
	{
		private static final double alpha = 1.0/12;
		
		private RealVector A, B, C;
		
		public double rating_DG[]; // [item]
		
		public DGVoltageRegADP()
		{
			costFunction = new DGVoltageRegCost();
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
			int dimension = x_t.getDimension();

			/*
			ArrayRealVector mean_DG = new ArrayRealVector(meanP_DG[t]);
			ArrayRealVector one = new ArrayRealVector(dimension, 1);
			RealVector timeOn = new ArrayRealVector(dimension);
			for(int i = 0 ; i < dimension; ++i)
			{
				double mean = mean_DG.getEntry(i);
				if(mean > 0)
					timeOn.setEntry(i, u_t.getEntry(i)/mean);
				else
					timeOn.setEntry(i, 1.0);
			}
			RealVector off = one.subtract(timeOn); // 1 - P_{DG,t}/P^+_{DG,t}
			
			return x_t.add(off);
			*/
			
			if(A == null)
			{
				A = vector(dimension, 1-alpha);
				B = new ArrayRealVector(dimension);
				for(int i = 0; i < dimension; ++i)
				{
					B.setEntry(i, alpha/rating_DG[i]);
				}
				C = B;
			}
			
			return 
					A.ebeMultiply(x_t).add(
					B.ebeMultiply(u_t)).add(
					C.ebeMultiply(w_t.getSubVector(0, dimension)));
		}
		
		@Override
		public ADP makeNew()
		{
			return new DGVoltageRegADP();
		}
		
		@Override
		public SampleBasedEstimator newEstimator(double[] defaultBandwidth)
		{
			QuadraticEstimator qe = new QuadraticEstimator();
			qe.setMaxSampleCount(60);
			return qe;
		}
	}
	
	public static class DGVoltageRegCentralController extends GridlessCentralCoordinator
	{
		@Override
		public void fillInForecastAndInitialState(GridlessADP adp_B, List<Integer> sortedB)
		{
			super.fillInForecastAndInitialState(adp_B, sortedB);
			
			// Set x_0s:
			adp_B.setX0(new ArrayRealVector(adp_B.getGridlessData().dgCount)); // Default all to 0 at time 0.
			
			// Set DG ratings:
			int i_new = 0;
			((DGVoltageRegADP)adp_B).rating_DG = new double[adp_B.getGridlessData().dgCount];
			for (Integer i : sortedB)
			{
				if(i >= gridlessADP().getGridlessData().dgCount)
					break;
				
				((DGVoltageRegADP)adp_B).rating_DG[i_new] = ((DGVoltageRegADP)localControllerTemplate).rating_DG[i];
				++i_new;
			}
		}
		
		@Override
		public double[][] makeBandwidthArray(int storageCount, int dgCount, int loadCount, double socBW, double dgBW, double storageBW, double loadBW)
		{
			int T = localControllerTemplate.getT();
			double[][] bw = new double[T][dgCount];
			for(int t = 0; t < T; ++t)
				for(int i = 0; i < dgCount; ++i)
				{
					bw[t][i] = (t+1)*dgBW; // t = 0 -> T-1 gives bw = gdBW -> T*dgBW
				}
			return bw;
		}
	}
	
	public static class DGVoltageRegCost implements CostFunction
	{
		@Override
		public double g(ADP adp, int t, RealVector x, RealVector u)
		{
			RealVector one = VectorHelper.one(x.getDimension());
			RealVector oneMinusX = one.subtract(x);
			return oneMinusX.dotProduct(oneMinusX); // (1-x)^2
//			double g = 0;
//			int dimension = x.getDimension();
//			for(int i = 0; i < dimension; ++i)
//			{
//				g += Math.exp(-x.getEntry(i));
//			}
//			return g;
		}
	}
	
	public static void main(String[] args)
	{
		assertAssertsOn();
		CS04_DGVoltageReg cs04 = new CS04_DGVoltageReg();
		cs04.parseArgs(args);
		cs04.run();
	}

	/**
	 * The test.
	 */
	public void run()
	{
		GridlessCentralCoordinator cc = new DGVoltageRegCentralController();
		cc.localControllerTemplate = makeGrid();
		((GridlessADP)cc.localControllerTemplate).setNr(60);
		cc.C_max = 2;
		cc.U_max = 4;
		cc.alpha = 0.2;
		int iterations = 20;
		double probabilityOfUpdate = 1;
		double h = 2.0/12;
		
		logCCConfig(cc, iterations, probabilityOfUpdate, h);
		
		// Random comparison:
		double[] averageCTGs = randomCase(cc.gridlessADP());
		
		// DP comparison:
		boolean disableDP = true;
		if(!disableDP)
			dpCase((DiscreteGridlessADP)cc.gridlessADP(), "/opt/energynet/cs04/dpSchedule");

		// Full network ADP results:
		boolean disableADP = true;
		double[][] bw = cc.makeBandwidthArray(0, cc.gridlessADP().getGridlessData().dgCount, 0, 0, h, 0, 0);
		if(!disableADP)
		{
			adpCase((DiscreteGridlessADP)cc.gridlessADP(), bw);
		}
		else
		{
			System.out.println("Full network ADP comparison not run.\n");
		}
		
		// Semi-distributed ADP results:
		iterativeADPCTGs = new double[iterations][]; // [loop iteration][time]
		iterativeADPControls = new RealVector[iterations][]; // [loop iteration][time]
		iterativeADPExternalVoltages = new RealVector[iterations][][]; // [loop iteration][decomposition][time]
		semiDistributedCase(cc, h, iterations, probabilityOfUpdate);

		// Get full ADP schedule:
		int T = cc.localControllerTemplate.getT();
		ThreadPool.getInstance().waitForAll();
		double[] adpCTGs = new double[T+1];
		RealVector[] adpStates = disableADP ? null : new RealVector[T+1];
		RealVector[] adpNoise = disableADP ? null : new RealVector[T+1];
		RealVector[] adpControls = disableADP ? null : cc.gridlessADP().schedule(adpCTGs, adpStates, adpNoise);

		
		//// Logging: ////
		
		System.out.println("Costs-to-go:");
		logCTGs(T, averageCTGs, iterations, iterativeADPCTGs, adpCTGs);
		
		System.out.println(",Epsilon,Delta v at pilot t=1 k=[0.."+iterations+"]...");
		for(int i = 0; i < cc.decomposition.length; ++i)
		{
			Integer pilotIndex = cc.pilots.get(i);
			Integer pilotIndex_local = cc.pilots_i.get(i);
			String pilotName = cc.gridlessADP().getGridlessData().controlNames.get(pilotIndex);
			System.out.print(pilotName);
			System.out.print(",");
			System.out.print(cc.epsilon(i));
			System.out.print(",");
			
			for(int k = 0; k < iterativeADPExternalVoltages.length; ++k)
			{
				System.out.print(iterativeADPExternalVoltages[k][i][1].getEntry(pilotIndex_local)); // Pilot bus at time 1.
				System.out.print(",");
			}
			
			System.out.println();
		}

		System.out.println("-");
		System.out.println("**** ADP ****");
		{
			logControlStateNoiseVoltage(
					adpStates,
					adpControls,
					adpNoise,
					cc.gridlessADP().getGridlessData().controlNames,
					cc.gridlessADP().getGridlessData().stateNames,
					cc.gridlessADP().getGridlessData().noiseNames,
					cc.gridlessADP().getGridlessData().voltageNames,
					cc.gridlessADP());
		}
		
		if(!disableDP)
		{
			System.out.println("-");
			System.out.println("**** DP ****");
			{
				RealVector[] x = new RealVector[cc.gridlessADP().getT()+1];
				RealVector[] u = new RealVector[cc.gridlessADP().getT()];
				RealVector[] w = new RealVector[cc.gridlessADP().getT()];
				fillStateControlNoise(cc.gridlessADP(), optimalSchedule, x, u, w);
				logControlStateNoiseVoltage(
						x,
						u,
						w,
						cc.gridlessADP().getGridlessData().controlNames,
						cc.gridlessADP().getGridlessData().stateNames,
						cc.gridlessADP().getGridlessData().noiseNames,
						cc.gridlessADP().getGridlessData().voltageNames,
						cc.gridlessADP());
			}
		}
		
		System.out.println("-");
		System.out.println("**** Semi-distributed ADP ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			System.out.println("---- Iteration "+k+" ----");
			RealVector[] x = new RealVector[cc.gridlessADP().getT()+1];
			RealVector[] w = new RealVector[cc.gridlessADP().getT()];
			fillStateAndNoiseFromControl(iterativeADPControls[k], cc.gridlessADP(), x, w);
			logControlStateNoiseVoltage(
					x,
					iterativeADPControls[k],
					w,
					cc.gridlessADP().getGridlessData().controlNames,
					cc.gridlessADP().getGridlessData().stateNames,
					cc.gridlessADP().getGridlessData().noiseNames,
					cc.gridlessADP().getGridlessData().voltageNames,
					cc.gridlessADP());
		}

		// Log external delta voltages:
		System.out.println("-");
		System.out.println("**** External Delta V ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			System.out.println("---- Iteration "+k+" ----");
			for(int d = 0; d < cc.decomposition.length; ++d)
			{
				System.out.println("---- Decomposition "+d+" ----");
				logDataWithNamesAndTimes(cc.gridlessADP().getT(), iterativeADPExternalVoltages[k][d], cc.decomposition(d).getGridlessData().controlNames);
			}
		}
		
		System.exit(0);
	}
	
	public GridlessADP makeGrid()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
//		GridDisplay.showInFrame(grid, IEEE13BusGrid.BASE_POWER, IEEE13BusGrid.BASE_VOLTAGE);
		
//		grid.bus632.addChild(makeStorage("S632"));
//		grid.bus671.addChild(makeStorage("S671"));
//		grid.bus684.addChild(makeStorage("S684"));
//		grid.bus645.addChild(makeStorage("S645"));

		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
//		grid.bus634.addChild(makeDG("DG7"));
//		grid.bus671.addChild(makeDG("DG8"));
		
		GridlessADP adp = makeGridlessADP(grid, new DGVoltageRegADP());
		
		// State is accumulated time spent turned off:
		adp.setX0(new ArrayRealVector(adp.getGridlessData().dgCount)); // Default all to 0 at time 0.
		
		// Set DG sizes:
		double[] ratings = new double[adp.getGridlessData().dgCount];
		for(int t = 0; t < adp.getT(); ++t)
		{
			for(int i = 0; i < adp.getGridlessData().dgCount; ++i)
			{
				if(adp.getGridlessData().meanP_DG[t][i] > ratings[i])
					ratings[i] = adp.getGridlessData().meanP_DG[t][i];
			}
		}
		((DGVoltageRegADP)adp).rating_DG = ratings;
		
		// Set state names to be DG:
		int i = 0;
		for (DistributedSource dg : grid.getDistributedSources(false))
		{
			adp.getGridlessData().stateNames.put(i, dg.getName());
			++i;
		}
		
		return adp;
	}
}
