package ellipsis.energy.test;

import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Timer;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;
import ellipsis.util.Pair;

public class GridlessADPIEEE13BusGridTest extends GridTestHelper
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		GridlessADPIEEE13BusGridTest tests = new GridlessADPIEEE13BusGridTest();
		tests.parseArgs(args);
		tests.test_IEEE13BusGrid();
		tests.test_IEEE13BusGrid_4DG2S();
		tests.test_IEEE13BusGrid_6DG2S();
		
		tests.printRandomPath();
	}

	private void printRandomPath()
	{
		if(verbose && iterations > 0)
		{
			System.out.println();
			GridlessADP adp = makeGridlessADPWith2DG2S();
			double[] randomCtgSum = new double[adp.getT()+1];
			printAverageRandomPath((DiscreteGridlessADP)adp, randomCtgSum);
		}
	}

	private void test_IEEE13BusGrid()
	{
		startTest();
		
		GridlessADP adp = makeGridlessADPWith2DG2S();
		
		// Deterministic comparison:
		if(verbose)
			System.out.println();
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)adp);
		
		// Train approximations:
		double[][] h = new double[adp.getT()][adp.getGridlessData().storageCount];//+2*(adp.dgCount+adp.storageCount+adp.loadCount)];
		for(int t = 0; t < adp.getT(); ++t)
			for(int i = 0; i < adp.getGridlessData().storageCount; ++i)
				h[t][i] = 0.5; // SOC
//		for(int i = 1; i < h.length; ++i)
//			h[i] = 0.5; // DG, storage, loads
		
		if(verbose)
			Timer.getTimer("train").start();
		
		adp.train(h);
		
		if(verbose)
			System.out.println("GridlessADP.train() ran in "+Timer.getTimer("train").stop()+"ms");
		
		// Get approximate schedule:
		double[] approximateCtgs = new double[adp.getT()+1];
		RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
		
		// Verbose output:
		printScheduleComparisonAndVoltages((DiscreteGridlessADP)adp, schedule, approximateCtgs, approximateSchedule);

		// Check that it's the same as the deterministic case:
		int t = 0;
		Double ctg_t = schedule.get(t).getValue();
		assertTrue( equals(ctg_t, approximateCtgs[t], 1) , "test_IEEE13BusGrid(): Schedule incorrect at time "+t+", cost-to-go is "+ctg_t+", approx. was "+approximateCtgs[t]);
		
		endTest("test_IEEE13BusGrid()");
	}

	private void test_IEEE13BusGrid_4DG2S()
	{
		startTest();
		
		// Setup grid with storage and DG:
		IEEE13BusGrid grid = makeGridWith2DG2S();
		grid.bus684.addChild(makeDG("DG3"));
		grid.bus645.addChild(makeDG("DG4"));
		
		// Make GridlessADP from grid:
		GridlessADP adp = makeGridlessADP(grid);
		
		// Deterministic comparison:
		if(verbose)
			System.out.println();
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)adp);
		
		// Train approximations:
		long start = System.currentTimeMillis();
		
		double[][] h = new double[adp.getT()][adp.getGridlessData().storageCount];//+2*(adp.dgCount+adp.storageCount+adp.loadCount)];
		for(int t = 0; t < adp.getT(); ++t)
			for(int i = 0; i < adp.getGridlessData().storageCount; ++i)
				h[t][i] = 0.5; // SOC
//		for(int i = 1; i < h.length; ++i)
//			h[i] = 0.5; // DG, storage, loads
		adp.train(h);
		
		long time = System.currentTimeMillis() - start;
		if(verbose)
			System.out.println("GridlessADP.train() ran in "+time+"ms");
		
		// Get approximate schedule:
		double[] approximateCtgs = new double[adp.getT()+1];
		RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
		
		// Verbose output:
		printScheduleComparisonAndVoltages((DiscreteGridlessADP)adp, schedule, approximateCtgs, approximateSchedule);

		// Check that it's the same as the deterministic case:
		int t = 0;
		Double ctg_t = schedule.get(t).getValue();
		assertTrue( equals(ctg_t, approximateCtgs[t], 1) , "test_IEEE13BusGrid_4DG4S(): Schedule incorrect at time "+t+", cost-to-go is "+ctg_t+", approx. was "+approximateCtgs[t]);
		
		endTest("test_IEEE13BusGrid_4DG4S()");
	}

	private void test_IEEE13BusGrid_6DG2S()
	{
		startTest();
		
		GridlessADP adp = makeGridlessADPWith6DG2S(false);
		
		// Deterministic comparison:
		if(verbose)
			System.out.println();
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)adp);
		
		// Train approximations:
		long start = System.currentTimeMillis();
		
		double[][] h = new double[adp.getT()][adp.getGridlessData().storageCount];//+2*(adp.dgCount+adp.storageCount+adp.loadCount)];
		for(int t = 0; t < adp.getT(); ++t)
			for(int i = 0; i < adp.getGridlessData().storageCount; ++i)
				h[t][i] = 0.5; // SOC
//		for(int i = 1; i < h.length; ++i)
//			h[i] = 0.5; // DG, storage, loads
		adp.train(h);
		
		long time = System.currentTimeMillis() - start;
		if(verbose)
			System.out.println("GridlessADP.train() ran in "+time+"ms");
		
		// Get approximate schedule:
		double[] approximateCtgs = new double[adp.getT()+1];
		RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
		
		// Verbose output:
		printScheduleComparisonAndVoltages((DiscreteGridlessADP)adp, schedule, approximateCtgs, approximateSchedule);

		// Check that it's the same as the deterministic case:
		int t = 0;
		Double ctg_t = schedule.get(t).getValue();
		assertTrue( equals(ctg_t, approximateCtgs[t], 1) , "test_IEEE13BusGrid_6DG4S(): Schedule incorrect at time "+t+", cost-to-go is "+ctg_t+", approx. was "+approximateCtgs[t]);
		
		endTest("test_IEEE13BusGrid_6DG4S()");
	}

	public static GridlessADP makeGridlessADPWith6DG2S(boolean showInFrame)
	{
		// Setup grid with storage and DG:
		IEEE13BusGrid grid = makeGridWith6DG2S();

		if(showInFrame)
			GridDisplay.showInFrame(grid, BASE_POWER, BASE_VOLTAGE);
		
		// Make GridlessADP from grid:
		GridlessADP adp = makeGridlessADP(grid);
		
		return adp;
	}

	public static IEEE13BusGrid makeGridWith6DG2S()
	{
		IEEE13BusGrid grid = makeGridWith2DG2S();
		
		grid.bus684.addChild(makeDG("DG3"));
		grid.bus645.addChild(makeDG("DG4"));
		grid.bus611.addChild(makeDG("DG5"));
		grid.bus646.addChild(makeDG("DG6"));
		
		return grid;
	}

	public static GridlessADP makeGridlessADPWith2DG2S()
	{
		// Setup grid with storage and DG:
		IEEE13BusGrid grid = makeGridWith2DG2S();

		// Make GridlessADP from grid:
		GridlessADP adp = makeGridlessADP(grid);
		
		return adp;
	}

	public static GridlessADP makeGridlessADP(IEEE13BusGrid grid)
	{
		return makeGridlessADP(grid, new DiscreteGridlessADP());
	}
	
	public static GridlessADP makeGridlessADP(IEEE13BusGrid grid, GridlessADP adp)
	{
		GridlessADPFactory.fromGridToGridless(grid, adp, IEEE13BusGrid.BASE_VOLTAGE, IEEE13BusGrid.BASE_POWER, /*failIfNotConvergent=*/true);

		// Setup schedule:
		adp.getGridlessData().setMeanP_DG(new double[] {1.4,1.55,1.6,1.6,1.5,1.25,0.85,0.55,0.2,0.1,0,0});
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
		return adp;
	}

	protected static IEEE13BusGrid makeGridWith2DG2S()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		
		ControllableDemand s1 = makeStorage("S1");
		grid.bus632.addChild(s1);

		ControllableDemand s2 = makeStorage("S2");
		grid.bus671.addChild(s2);

		DistributedSource dg1 = makeDG("DG1");
		grid.bus680.addChild(dg1);

		DistributedSource dg2 = makeDG("DG2");
		grid.bus675.addChild(dg2);
		
		return grid;
	}

	public static DistributedSource makeDG(String name)
	{
		DistributedSource dg = new DistributedSource();
		dg.setName(name);
		dg.setPmax(1400e3);
		dg.setPowerOutput(1400e3, 0.0);
		return dg;
	}

	public static ControllableDemand makeStorage(String name)
	{
		ControllableDemand s1 = new ControllableDemand();
		s1.setName(name);
		s1.setMaxCapacity(1000e3);
		s1.setCapacity(500e3);
		s1.setMaxChargeRate(400e3);
		s1.setMaxDischargeRate(400e3);
		s1.setChargeRate(new Complex(0.0));
		return s1;
	}
}
