package ellipsis.energy.casestudy;

import static ellipsis.energy.casestudy.TestCaseLogger.logCCConfig;
import static ellipsis.energy.casestudy.TestCaseLogger.logControlStateNoiseVoltage;

import java.io.IOException;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;
import com.mls.util.Timer;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.casestudy.CS01.CS01CostFunction;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.grid.GridGenerator;
import ellipsis.energy.test.GridTestHelper;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.energy.util.log.LoggedGridlessADP;
import ellipsis.energy.util.log.LoggedGridlessCentralController;

public class CS03_LargeNets extends GridTestHelper
{
	public static final boolean ESTIMATE_PDFS = false;
	public static final boolean ESTIMATE_PDFS_testUncontrolledDGs = false;
	public static final boolean ESTIMATE_PDFS_testPassiveGrid = false;
	public static final boolean SHOW_IN_FRAME = false;
	public static final boolean TEST_MANY = true;
	public static final boolean ALLOW_INCREMENTAL_CONSOLE_OUTPUT = false;
	
	public static final double BASE_VOLTAGE = 11e3;
	public static void main(String[] args)
	{
		/*
		  4. Network with controlled storage:
		    1. Use uncontrolled DG case - check cost, voltages, and level of output.
		  5. Network with controlled DG and storage:
		    1. Check cost, voltages, and level of output. 
		  6. Semi-distributed problem:
		    1. Network with controlled DG,
		    2. Network with controlled storage,
		    3. Network with controlled DG and storage.
		 */

		CS03_LargeNets cs03 = new CS03_LargeNets();
		if(ESTIMATE_PDFS)
		{
			cs03.estimatePDFs();
		}
		else if(TEST_MANY)
		{
			int testsPerCase = 5;
			int seed = 36;
			for(int i = 0; i < 100; ++i)
			{
				int minSlackLines = 4 + i/testsPerCase;
				int maxSlackLines = minSlackLines + 10;
				cs03.testControlledDGsAndStorage(seed, minSlackLines, maxSlackLines);
				++seed;
			}
		}
		else
		{
//			cs03.testPassiveGrid(0);
//			cs03.testUncontrolledDGs(0);
//			cs03.testControlledDGs(0);
//			cs03.testControlledDGsAndStorage(0, 4, 8);
			cs03.testControlledDGsAndStorage(0, 6, 10);
//			cs03.testControlledDGsAndStorage(0, 10, 14);
//			cs03.testControlledDGsAndStorage(0, 15, 20);
//			cs03.testControlledDGsAndStorage(0, 25, 35);
		}
	}

	/**
	 * 
	 */
	@SuppressWarnings("unused")
	private void estimatePDFs()
	{
		if(ESTIMATE_PDFS_testPassiveGrid)
		{
			System.out.println("Passive Grid:\nMin V, Max V");
			for(int seed = 0; seed < 10000; ++seed)
			{
				testPassiveGrid(seed);
			}
		}
		 
		if(ESTIMATE_PDFS_testPassiveGrid && ESTIMATE_PDFS_testUncontrolledDGs)
		{
			System.out.println("Press enter to continue...");
			try
			{
				System.in.read();
			} catch (IOException e)
			{
				e.printStackTrace();
			}
		}
		 
		if(ESTIMATE_PDFS_testUncontrolledDGs)
		{
			System.out.println("Uncontrolled DG:\nMin V, Max V");
			for(int seed = 0; seed < 10000; ++seed)
			{
				testUncontrolledDGs(seed);
			}
		}
	}
	
	/**
	 * 3. Network with controlled DG:
		    1. Start with network exhibiting excessive voltages - check cost function in this case.
		    2. Add control in order to stop voltage breaches - check cost, voltages, and level of output. Note that this is not a DP problem (no dynamics).
	 */
	@SuppressWarnings("unused")
	private void testControlledDGs(int seed)
	{
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.C_max = 2;
		cc.U_max = 4;
		cc.alpha = 0.5;
		int iterations = 20;

		GridGenerator gen = makeGeneratorForControlledDG();
		Grid grid = gen.createGrid("Controlled DG", seed);
		
		GridDisplay.showInFrame(grid, BASE_POWER, BASE_VOLTAGE);
		
		testCC(cc, grid, iterations);
	}

	@SuppressWarnings("unused")
	public void testCC(GridlessCentralCoordinator cc, Grid grid, int iterations)
	{
		cc.localControllerTemplate = GridlessADPFactory.fromGridToGridless(grid, new LoggedGridlessADP(), BASE_VOLTAGE, BASE_POWER, true);
		((GridlessADP)cc.localControllerTemplate).setNr(60);
		((GridlessADP)cc.localControllerTemplate).setCostFunction(new CS01CostFunction());
		setupADP(cc.gridlessADP());

		double socBW = 0.5;
		logCCConfig(cc, iterations, 1, socBW);
		
		// Random comparison:
		double[] averageCTGs = null;//randomCase(adp);
		
		// DP comparison:
//		dpCase(adp, "/opt/energynet/cs03/dpSchedule");
		
		// Full network ADP results:
//		adpCase(adp, socBW);
		
		// Semi-distributed ADP results:
		iterativeADPCTGs = new double[iterations][]; // [loop iteration][time]
		iterativeADPControls = new RealVector[iterations][]; // [loop iteration][time]
		iterativeADPExternalVoltages = new RealVector[iterations][][]; // [loop iteration][decomposition][time]
		semiDistributedCase(cc, socBW, iterations, 1);
		
		// Wait for results:
		ThreadPool.getInstance().waitForAll();
		
		// Logging:
//		double[] adpCTGs = new double[cc.adp.T+1];
//		RealVector[] adpStates = new RealVector[cc.adp.T+1];
//		RealVector[] adpNoise = new RealVector[cc.adp.T+1];
//		/*RealVector[] adpControls = */cc.adp.schedule(adpCTGs, adpStates, adpNoise);
		int T = cc.localControllerTemplate.getT();
		logCTGs(T, averageCTGs, iterations, iterativeADPCTGs, null);//adpCTGs);

		int t = 9;
		System.out.println(",Epsilon,Delta v at pilot t="+t+" k=[0.."+iterations+"]...");
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
				System.out.print(iterativeADPExternalVoltages[k][i][t].getEntry(pilotIndex_local)); // Pilot bus at time t.
				System.out.print(",");
			}
			
			System.out.println();
		}
		
		System.out.println("-");
		System.out.println("**** Semi-distributed ADP ****");
		for(int k = 0; k < iterativeADPControls.length; ++k) // iterations
		{
			System.out.print("---- Iteration "+k+" ----");
			if(ALLOW_INCREMENTAL_CONSOLE_OUTPUT && grid.getUnitCount() > 200)
			{
				System.out.print("(press enter to continue)");
				try
				{
					System.in.read();
				} catch (IOException e)
				{
					e.printStackTrace();
				}
			}
			else
			{
				System.out.println();
			}
			
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
	}
	
	private void testControlledDGsAndStorage(int seed, int minSlackLines, int maxSlackLines)
	{
		GridlessCentralCoordinator cc = new LoggedGridlessCentralController();
		cc.C_max = 2;
		cc.U_max = 4;
		cc.alpha = 0.1;
		int iterations = 20;

		GridGenerator gen = makeGeneratorForControlledDGAndStorage(minSlackLines, maxSlackLines);
		Grid grid = gen.createGrid("Controlled DG & Storage", seed);

		if(SHOW_IN_FRAME)
			GridDisplay.showInFrame(grid, BASE_POWER, BASE_VOLTAGE);
		
		Timer timer = Timer.startNewTimer();
		testCC(cc, grid, iterations);
		System.out.println("testCC time was "+timer.stop());
	}

	private void setupADP(GridlessADP adp)
	{
		// Setup schedule:
		adp.getGridlessData().setMeanP_DG(
//				new double[] {1.4, 1.55, 1.6, 1.6, 1.5, 1.25, 0.85, 0.55, 0.2, 0.1, 0, 0});
//				new double[] {0.875,	0.96875,	1,	1,	0.9375,	0.78125,	0.53125,	0.34375,	0.125,	0.0625,	0,	0});
//				new double[] {0.4375,	0.484375,	0.5,	0.5,	0.46875,	0.390625,	0.265625,	0.171875,	0.0625,	0.03125,	0,	0});
//				new double[] {0.004375,	0.00484375,	0.005,	0.005,	0.0046875,	0.00390625,	0.00265625,	0.00171875,	0.000625,	0.0003125,	0,	0});
				new double[] {0.013125,	0.01453125,	0.015,	0.015,	0.0140625,	0.01171875,	0.00796875,	0.00515625,	0.001875,	0.0009375,	0,	0});


		adp.getGridlessData().setsdP_DG(new double[] {0.02, 0.02, 0.02, 0.02, 0.015, 0.012, 0.01, 0.005, 0.002, 0.01, 0.0, 0.0});

		// Set load schedule:
		double[] relativeLoadSchedule = new double[]{1.0, 1.0, 0.8, 0.7, 0.9, 1.2, 1.4, 1.2, 0.7, 0.5, 0.4, 0.3};
		
		int pOffset = adp.getGridlessData().dgCount + adp.getGridlessData().storageCount; // Past P_DGs, P_Ss
		int qOffset = pOffset + adp.getGridlessData().loadCount + adp.getGridlessData().dgCount + adp.getGridlessData().storageCount; // Past P_DGs, P_Ss, P_Ls, Q_DGs, Q_Ss
		
		adp.getGridlessData().meanP_L = makeSchedule(adp.getGridlessData().power_0, pOffset, adp.getGridlessData().loadCount, adp.getT(), relativeLoadSchedule);
		adp.getGridlessData().meanQ_L = makeSchedule(adp.getGridlessData().power_0, qOffset, adp.getGridlessData().loadCount, adp.getT(), relativeLoadSchedule);
		
		adp.getGridlessData().setsdP_L(new double[] {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
		adp.getGridlessData().setsdQ_L(new double[] {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
		
		adp.setNr(100);
		
		adp.setCostFunction(new CS01CostFunction());
	}

	/**
	 * 2. Network with uncontrolled DG:
		    2. Operating with excess DG output at some times; voltage rises above limits.
	 */
	private void testUncontrolledDGs(int seed)
	{
		GridGenerator gen = makeGeneratorForUncontrolledDG();
		Grid grid = gen.createGrid("Uncontrolled DG", seed);
		
		if(!ESTIMATE_PDFS_testUncontrolledDGs)
		{
			GridDisplay.showInFrame(grid, BASE_POWER, BASE_VOLTAGE);
			System.out.println("Uncontrolled DG:\nMin V, Max V");
		}
		
		analyseGrid(grid);
	}

	/**
	 * 1. Network with no DG or storage.
		    1. Operating normally, always within voltage constraints.
	 */
	private void testPassiveGrid(int seed)
	{
		GridGenerator gen = makeGeneratorForPassiveNetwork();
		Grid grid = gen.createGrid("Passive Network", seed);
		
		if(!ESTIMATE_PDFS_testPassiveGrid)
		{
			GridDisplay.showInFrame(grid, BASE_POWER, BASE_VOLTAGE);
			System.out.println("Passive Network:\nMin V, Max V");
		}

		analyseGrid(grid);
	}
	
	private void analyseGrid(Grid grid)
	{
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(BASE_POWER);
        lfa.setBaseVoltage(BASE_VOLTAGE);
        lfa.setIterations(100);
        lfa.setTargetError(1e-6);
        
        AnalysisResults results = lfa.analyse();
        Complex minV = Complex.INF;
        Complex maxV = Complex.ZERO;
        for (Bus bus : grid.getBusses())
		{
        	String name = bus.getName();
			Complex v = results.getBusVoltage(name);
			if(v.abs() > maxV.abs())
				maxV = v;
			if(v.abs() < minV.abs())
				minV = v;
		}

        System.out.println(minV.abs()+","+maxV.abs());
	}
	
	private GridGenerator makeGeneratorForPassiveNetwork()
	{
		GridGenerator gen = new GridGenerator();
		gen.slackVoltage = BASE_VOLTAGE;
		
		// PV:
		gen.pvDensity = 0;
		
		// Storage:
		gen.storageDensity = 0;
		
		// Loads:
		gen.minInitialLoad = 120e3;
		gen.maxInitialLoad = 1200e3;
		
		// Busses:
		gen.maxBusDepth = 4;
		
		// Lines:
		gen.minLineResistance = 0.1;
		gen.maxLineResistance = 0.7;
		gen.minLineInductance = 0.01;
		gen.maxLineInductance = 0.07;
		gen.minLineLength = 100e-3; // km
		gen.maxLineLength = 1000e-3; // km
		
		return gen;
	}

	private GridGenerator makeGeneratorForUncontrolledDG()
	{
		GridGenerator gen = makeGeneratorForControlledDG();
		
		// Consider the case where load has dipped:
		gen.minInitialLoad = gen.minInitialLoad*0.3;
		gen.maxInitialLoad = gen.maxInitialLoad*0.3;
		
		return gen;
	}

	public GridGenerator makeGeneratorForControlledDG()
	{
		GridGenerator gen = makeGeneratorForPassiveNetwork();
		
		// DG:
		gen.pvDensity = 0.9;
		gen.minPvCapacity = 1500e3;
		gen.maxPvCapacity = 1500e3;

//		// Loads:
//		gen.minInitialLoad = 500e3;
//		gen.maxInitialLoad = 500e3;
		
		return gen;
	}

	public GridGenerator makeGeneratorForControlledDGAndStorage(int minSlackLines, int maxSlackLines)
	{
		GridGenerator gen = makeGeneratorForPassiveNetwork();
		
		// Network size:
		gen.maxLineCountFromSlack = maxSlackLines;
		gen.minLineCountFromSlack = minSlackLines;
		
		// DG:
		gen.pvDensity = 0.9;
		gen.minPvCapacity = 1500e3;
		gen.maxPvCapacity = 1500e3;

		// Storage:
		gen.storageDensity = 0.3;
		gen.minStorageCapacity = 1000e3;
		gen.maxStorageCapacity = 2000e3;
		gen.minStorageChargeRate = 500e3;
		gen.maxStorageChargeRate = 1000e3;
		
		return gen;
	}
}
