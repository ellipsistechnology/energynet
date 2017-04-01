package ellipsis.energy.casestudy;

import static ellipsis.energy.casestudy.TestCaseLogger.logCCConfig;
import static ellipsis.energy.casestudy.TestCaseLogger.logControlStateNoiseVoltage;
import static ellipsis.energy.casestudy.TestCaseLogger.logDataWithNamesAndTimes;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeDG;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeGridlessADP;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeStorage;

import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;

import ellipsis.energy.test.GridTestHelper;
import ellipsis.energy.test.IEEE13BusGrid;
import ellipsis.energy.util.ADP;
import ellipsis.energy.util.CostFunction;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.util.Pair;

/**
 * IEEE 13 bus network with 2 storage and 6 DG.
 * Test case compares:
 * 	Random,
 *  DP,
 *  Full network ADP,
 *  Semi-distributed ADP.
 * 
 * @author bmillar
 *
 */
public class CS01 extends GridTestHelper
{
	public static class CS01CostFunction implements CostFunction
	{
		@Override
		public double g(ADP adp, int t, RealVector x, RealVector u)
		{
			DiscreteGridlessADP gadp = (DiscreteGridlessADP)adp;
			
			// Terminating case with no control:
			if(u == null)
				return 0;
			
			// Normal case:
			RealVector w = adp.wZero(t);
			RealVector PQ = gadp.powerFromControlAndNoise(u, w);
			
			// Minimisation actually minimises change in current: deltaI = Y(A)\Lambda_{A,B}PQ(B):
			RealVector argAbs = gadp.getGridlessData().slackSensitivities.operate(PQ);
			Complex deltaV = new Complex(argAbs.getEntry(1), argAbs.getEntry(0));
			Complex deltaI = gadp.getGridlessData().slackAdmittance.multiply(deltaV);
			
			return -Math.signum(deltaI.getReal())*deltaI.abs();
		}
	}

	public static void main(String[] args)
	{
		assertAssertsOn();
		CS01 cs01 = new CS01();
		cs01.parseArgs(args);
		cs01.run();
	}

	/**
	 * The test.
	 */
	public void run()
	{
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = makeGrid();
		((GridlessADP)cc.localControllerTemplate).setNr(60);
		((GridlessADP)cc.localControllerTemplate).setCostFunction(new CS01CostFunction());
		cc.C_max = 2;
		cc.U_max = 4;
		cc.alpha = 0.8;
		int iterations = 20;
		double probabilityOfUpdate = 1;

		double socBW = 0.5;
		logCCConfig(cc, iterations, probabilityOfUpdate, socBW);
		
		// Random comparison:
		double[] averageCTGs = randomCase((DiscreteGridlessADP)cc.gridlessADP());
		
		// DP comparison:
		dpCase((DiscreteGridlessADP)cc.gridlessADP(), "/opt/energynet/cs01/dpSchedule");
		
		// Full network ADP results:
		boolean disableADP = false;
		if(!disableADP)
		{
			adpCase((DiscreteGridlessADP)cc.gridlessADP(), cc.makeBandwidthArray(cc.gridlessADP().getGridlessData().storageCount, 0, 0, socBW, 0, 0, 0));
		}
		else
		{
			System.out.println("Full network ADP comparison not run.\n");
		}
		
		// Semi-distributed ADP results:
		iterativeADPCTGs = new double[iterations][]; // [loop iteration][time]
		iterativeADPControls = new RealVector[iterations][]; // [loop iteration][time]
		iterativeADPExternalVoltages = new RealVector[iterations][][]; // [loop iteration][decomposition][time]
		semiDistributedCase(cc, socBW, iterations, probabilityOfUpdate);

		// Get full ADP schedule:
		ThreadPool.getInstance().waitForAll();
		int T = cc.localControllerTemplate.getT();
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
		
		System.out.println("-");
		System.out.println("**** DP ****");
		{
			RealVector[] x = new RealVector[T+1];
			RealVector[] u = new RealVector[T];
			RealVector[] w = new RealVector[T];
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
		
		System.out.println("-");
		System.out.println("**** Semi-distributed ADP ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			System.out.println("---- Iteration "+k+" ----");
			RealVector[] x = new RealVector[T+1];
			RealVector[] w = new RealVector[T];
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
				logDataWithNamesAndTimes(T, iterativeADPExternalVoltages[k][d], cc.decomposition(d).getGridlessData().controlNames);
			}
		}
	}

	@SuppressWarnings("unused")
	private RealVector[] controls(List<Pair<RealVector, Double>> optimalSchedule)
	{
		RealVector[] u = new RealVector[optimalSchedule.size()];
		for (int i = 0; i < u.length; i++)
		{
			u[i] = optimalSchedule.get(i).getKey();
		}
		return u;
	}

	public GridlessADP makeGrid()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
//		GridDisplay.showInFrame(grid, IEEE13BusGrid.BASE_POWER, IEEE13BusGrid.BASE_VOLTAGE);
		
		grid.bus632.addChild(makeStorage("S632"));
		grid.bus671.addChild(makeStorage("S671"));
		grid.bus684.addChild(makeStorage("S684"));
		grid.bus645.addChild(makeStorage("S645"));

		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
//		grid.bus634.addChild(makeDG("DG7"));
//		grid.bus671.addChild(makeDG("DG8"));
		
		GridlessADP adp = makeGridlessADP(grid);
		
		return adp;
	}
}
