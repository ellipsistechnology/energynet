package ellipsis.energy.test;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.Load;
import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;
import ellipsis.util.Pair;
import ellipsis.util.VectorHelper;

public class GridlessADPTests2 extends GridTestHelper
{
    protected static final double BASE_POWER = 100e3;
    protected static final double BASE_VOLTAGE = 4.16e3;
    protected static final double BASE_IMPEDANCE = BASE_VOLTAGE*BASE_VOLTAGE/BASE_POWER;
    
	public static void main(String[] args)
	{
		assertAssertsOn();
		
		GridlessADPTests2 tests = new GridlessADPTests2();
		tests.parseArgs(args);
		
		tests.test_train_3Feeder4BusSystem(0);
		tests.test_train_3Feeder7BusSystem();
		tests.test_train_3Feeder10BusSystem();

		if(tests.iterations > 0)
		{
			double averageError = 0;
			for(int i = 0; i < tests.iterations; ++i)
			{
//				GridlessADP.rand = new Random(i);
				averageError += tests.test_train_3Feeder4BusSystem(i);
				averageError += tests.test_train_3Feeder7BusSystem();
				averageError += tests.test_train_3Feeder10BusSystem();
			}
			averageError /= tests.iterations*3;
			System.out.println("Average error = "+averageError);
		}
	}
	
	static final double RES_PER_KM = 0.7;
	static final double IND_PER_KM = 0.07;
	
	// All of the following are in p.u.:
	final static double L_P_0 = -170e3/BASE_POWER; // kW
	final static double L_Q_0 = -125e3/BASE_POWER; // kVAr
	
	final static double SOC_0 = 200e3/BASE_POWER; // kWh
	final static double SOC_MAX = 400e3/BASE_POWER; // kWh
	final static double S_CHARGE = 100e3/BASE_POWER; // kW
	final static double S_DISCHARGE = 125e3/BASE_POWER; // kW
	final static double S_P_0 = -S_CHARGE;
	
	final static double DG_MAX = 500e3/BASE_POWER; // kW
	final static double DG_P_0 = DG_MAX;

	private double test_train_3Feeder7BusSystem()
	{
		startTest();
		
		Grid grid = gridWith3Feeders();
		make7BusSystem(grid);
		double averageError = testGrid(grid, "Test 3 feeder, 7 bus system:");
		
		endTest("test_train_3Feeder7BusSystem()");
		
		return averageError;
	}

	private double test_train_3Feeder10BusSystem()
	{
		startTest();
		
		Grid grid = gridWith3Feeders();
		make10BusSystem(grid);

		double averageError = testGrid(grid, "Test 3 feeder, 10 bus system:");
		
		endTest("test_train_3Feeder10BusSystem()");
		
		return averageError;
	}

	private void make10BusSystem(Grid grid)
	{
		make7BusSystem(grid);
		
		addNewBextBus(grid, grid.getBus("6"), 9, 5);
		addNewBextBus(grid, grid.getBus("7"), 10, 6);
		addNewBextBus(grid, grid.getBus("8"), 11, 7);
		
		// Add DG to bus 11 added above:
		DistributedSource dg = new DistributedSource();
		dg.setName("DG2");
		dg.setPmax(DG_MAX*BASE_POWER);
		dg.setPowerOutput(DG_P_0*BASE_POWER, 0);
		grid.getBus("11").addChild(dg);
	}

	private double testGrid(Grid grid, String testName)
	{
		double localAverageError = 0;
		int errorCount = 0;
		
		// Make GridlessADP from grid:
		GridlessADP adp = fromGrid(grid);
		
		// Test counts:
		testCounts((DiscreteGridlessADP)adp);
		
		// Test powerFromControlAndNoise():
		RealVector u_0 = adp.getGridlessData().power_0.getSubVector(0, adp.getGridlessData().dgCount+adp.getGridlessData().storageCount); // power_0 is [P_DG P_S P_L Q_DG Q_S Q_L], u is [P_DG P_S]
		RealVector wZero = adp.wZero(0);
		RealVector Sfromx = adp.powerFromControlAndNoise(u_0, wZero);
		RealVector S = adp.getGridlessData().power_0;
		assertTrue(Sfromx.equals(S), testName+" powerFromControlAndNoise(u_0, w) returned incorrect vector "+Sfromx+", should be "+S);
		
		// Check correct deltaV_external:
		int voltageDimension = adp.getGridlessData().power_0.getDimension();
		for(int t = 0; t <= adp.getT(); ++t)
		{
			assertTrue(adp.getGridlessData().deltaV_external[t].getDimension() == voltageDimension, "Check correct deltaV_external: Incorrect dimension of "+adp.getGridlessData().deltaV_external[t].getDimension()+", should be "+voltageDimension);
			for(int i = 0; i < voltageDimension; ++i)
			{
				assertTrue(adp.getGridlessData().deltaV_external[t].getEntry(i) == 0, "Check correct deltaV_external: Entry "+i+" was not 0.0");
			}
		}
		
		// Deterministic comparison:
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)adp);
		
		// Train ADP:
		double[][] bw = new double[adp.getT()][1];
		for(int t = 0; t < adp.getT(); ++t)
			bw[t] = new double[]{0.5};
		adp.train(bw);
		
		// Get approximate schedule:
		double[] approximateCtgs = new double[adp.getT()+1];
		RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
		
		// Test that approximate schedule costs match GridlessADP.scheduleFromControls():
		double[] calculatedCtgs = adp.ctgs(approximateSchedule);
		assertTrue(calculatedCtgs.length == adp.getT()+1, "scheduleFromControls() returned array of length "+calculatedCtgs.length+", should be "+(adp.getT()+1));
		for(int t = 0; t < adp.getT(); ++t)
		{
			assertTrue(calculatedCtgs[t] == approximateCtgs[t], "scheduleFromControls(): CTG mismatch at time "+t+", scheduleFromControls() gave "+calculatedCtgs[t]+", schedule() gave "+approximateCtgs[t]);
		}
		
		// Output schedule:
		printScheduleComparisonAndVoltages((DiscreteGridlessADP)adp, schedule, approximateCtgs, approximateSchedule);

		// Check that it's the same as the deterministic case:
		int t = 0;
		Double ctg_t = schedule.get(t).getValue();
		assertTrue( equals(ctg_t, approximateCtgs[t], 1.5) , testName+" Schedule incorrect at time "+t+", cost-to-go is "+ctg_t+", approx. was "+approximateCtgs[t]);

		// Calculate average error:
		for(t = 0; t < adp.getT(); ++t)
		{
			localAverageError += percentageError(schedule, approximateCtgs, t);
			++errorCount;
		}
		double averageError = localAverageError/errorCount;
		return averageError;
	}

	private void make7BusSystem(Grid grid)
	{
		addNewBextBus(grid, grid.getBus("3"), 6, 2);
		addNewBextBus(grid, grid.getBus("4"), 7, 3);
		addNewBextBus(grid, grid.getBus("5"), 8, 4);
	}

	private void addNewBextBus(Grid grid, Bus bus, int nextBusNumber, int nextLoadNumber)
	{
		// New line and new bus:
		Line line = new Line();
		Bus nextBus = new Bus(grid);
		
		// Set up line:
		line.setName(bus.getName()+"-"+nextBusNumber);
		line.setFromBus(bus);
		line.setToBus(nextBus);
		line.setResistancePerMetre(RES_PER_KM);
		line.setInductancePerMetre(IND_PER_KM);
		line.setLength(1.0);
		
		// Setup bus:
		nextBus.setName(""+nextBusNumber);
		grid.add(nextBus);
		Load load = new Load();
		load.setName("L"+nextLoadNumber);
		load.setLoad(-L_P_0*BASE_POWER, -L_Q_0*BASE_POWER);
		nextBus.addChild(load);
		bus.addChild(line);
	}

	private double test_train_3Feeder4BusSystem(int seed)
	{
		startTest();
		double localAverageError = 0;
		int errorCount = 0;
		
		// Grid with 4 busses - one root bus and 3 1-bus feeders:
		Grid grid = gridWith3Feeders();
		
		// Make GridlessADP from grid:
		GridlessADP adp = fromGrid(grid);
		((DiscreteGridlessADP)adp).rand = new Random(seed);
		
		// Test counts:
		testCounts((DiscreteGridlessADP)adp);
		
		// Random comparison:
		double[] randomCtgSum = new double[adp.getT()];
		if(iterations > 0)
		{
			printAverageRandomPath((DiscreteGridlessADP)adp, randomCtgSum);
		}
		
		// Deterministic comparison:
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)adp);
		
		// Train ADP:
		double[][] bw = new double[adp.getT()][1];
		for(int t = 0; t < adp.getT(); ++t)
			bw[t] = new double[]{0.5};
		adp.train(bw);
		
		// Get approximate schedule:
		double[] approximateCtgs = new double[adp.getT()+1];
		RealVector[] approximateSchedule = adp.schedule(approximateCtgs);
		
		// Output schedule:
		if(verbose)
		{
			System.out.println("\ncorrect,,,approximate");
			for(int t = 0; t < adp.getT(); ++t)
			{
				System.out.println(VectorHelper.printVector(schedule.get(t).getKey())+schedule.get(t).getValue()+","+VectorHelper.printVector(approximateSchedule[t])+approximateCtgs[t]);
			}
		}

		// Check that it's the same as the deterministic case:
		int t = 0;
		Double ctg_t = schedule.get(t).getValue();
		assertTrue( equals(ctg_t, approximateCtgs[t], 1.5) , "Test 3 feeder, 4 bus system: Schedule incorrect at time "+t+", cost-to-go is "+ctg_t+", approx. was "+approximateCtgs[t]);
		
		// Compare each time step with average case:
		for(t = 0; t < adp.getT(); ++t)
		{
			if(iterations > 0)
				assertTrue( approximateCtgs[t] < randomCtgSum[t] , "Test 3 feeder, 4 bus system: Schedule is worse than random average at time "+t);
			localAverageError += percentageError(schedule, approximateCtgs, t);
			++errorCount;
		}
		
		endTest("test_train_3Feeder4BusSystem()");
		
		return localAverageError/errorCount;
	}

	private void testCounts(DiscreteGridlessADP adp)
	{
		int unitCount = adp.getGridlessData().dgCount+adp.getGridlessData().storageCount+adp.getGridlessData().loadCount;
		int stateDimension = adp.getGridlessData().storageCount;
		int controlableCount = adp.getGridlessData().dgCount+adp.getGridlessData().storageCount;
		assertTrue(adp.stateDimension() == stateDimension, "Test counts: Incorrect state dimension, was "+adp.stateDimension()+", should be "+stateDimension);
		assertTrue(adp.stateDimension() == adp.getX0().getDimension(), "Test counts: Incorrect x_0 dimension, was "+adp.getX0().getDimension()+", stateDimension() = "+adp.stateDimension()+", should be "+stateDimension);
		assertTrue(adp.voltageDimension() == 2*unitCount, "Test counts: Incorrect voltage dimension, was "+adp.voltageDimension()+", should be "+2*unitCount);
		assertTrue(adp.powerDimension() == 2*unitCount, "Test counts: Incorrect power dimension, was "+adp.powerDimension()+", should be "+2*unitCount);
		assertTrue(adp.unitCount() == unitCount, "Test counts: Incorrect unit count, was "+adp.unitCount()+", should be "+adp.unitCount());
		assertTrue(adp.controlDimension() == controlableCount, "Test counts: Incorrect control dimension, was "+adp.controlDimension()+", should be "+controlableCount);
	}

	private GridlessADP fromGrid(Grid grid)
	{
		GridlessADP adp = GridlessADPFactory.fromGrid(grid, BASE_VOLTAGE, BASE_POWER, /*failIfNotConvergent=*/true);

		// Setup schedule:
		adp.getGridlessData().setMeanP_DG(new double[] {0.3, 0.25, 0.22, 0.2, 0.18, 0.14, 0.1, 0.05, 0.02, 0.0, 0.0, 0.0});
		adp.getGridlessData().setsdP_DG(new double[] {0.02, 0.02, 0.02, 0.02, 0.015, 0.012, 0.01, 0.005, 0.002, 0.01, 0.0, 0.0});
		
		adp.getGridlessData().setMeanP_L(new double[] {L_P_0/*-0.3*/, -0.45, -0.4, -0.3, -0.25, -0.35, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1});
		adp.getGridlessData().setsdP_L(new double[] {0.03, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.01});
		
		adp.getGridlessData().setMeanQ_L(new double[] {L_Q_0/*-0.1*/, -0.15, -0.15, -0.1, -0.1, -0.1, -0.15, -0.1, -0.1, -0.05, -0.05, -0.05});
		adp.getGridlessData().setsdQ_L(new double[] {0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01});
		
		adp.setNr(100);
		
		return adp;
	}

	private Grid gridWith3Feeders()
	{
		Grid grid = Grid.grid().
				Bus("1").
				SlackSource("Slack", BASE_VOLTAGE, 0, 0, 0).
				Line("1-2", 0.600, RES_PER_KM, IND_PER_KM).
					Bus("2").
						Line("2-3", 0.150, RES_PER_KM, IND_PER_KM).
							Bus("3").
								Load("L1").
							terminate().
						terminate().
						Line("2-4", 0.150, RES_PER_KM, IND_PER_KM).
							Bus("4").
								DistributedSource("DG1").
							terminate().
						terminate().
						Line("2-5", 0.6, RES_PER_KM, IND_PER_KM).
							Bus("5").
								ControllableDemand("S").
							terminate().
						terminate().
					terminate().
				terminate().
			terminate().
		grid();
		
		// x_0:
		grid.getLoad("L1").setLoad(-L_P_0*BASE_POWER, -L_Q_0*BASE_POWER);
		
		grid.getControllableDemand("S").setMaxCapacity(SOC_MAX*BASE_POWER);
		grid.getControllableDemand("S").setCapacity(SOC_0*BASE_POWER);
		grid.getControllableDemand("S").setMaxChargeRate(S_CHARGE*BASE_POWER);
		grid.getControllableDemand("S").setMaxDischargeRate(S_DISCHARGE*BASE_POWER);
		grid.getControllableDemand("S").setChargeRate(-S_P_0*BASE_POWER, 0);
		
		grid.getDistributedSource("DG1").setPmax(DG_MAX*BASE_POWER);
		grid.getDistributedSource("DG1").setPowerOutput(DG_P_0*BASE_POWER, 0);
		
		return grid;
	}
}
