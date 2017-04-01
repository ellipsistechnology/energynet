package ellipsis.energy.test;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessADPFactory;

public class GridlessADPFactoryTest extends GridTestHelper
{ 
	private static final double L_P_0 = -0.45;
	private static final double L_Q_0 = -0.15;
	
	private static final double SOC_0 = 0.5;
	private static final double S_MAX = 2;
	private static final double S_MAX_CHARGE = -0.15;
	private static final double S_MAX_DISCHARGE = 0.25;
	private static final double S_P_0 = -0.15;
	
	private static final double DG_MAX = 0.2;
	private static final double DG_P_0 = 0.2;
	
	public static void main(String[] args)
	{
		assertAssertsOn();
		GridlessADPFactoryTest tests = new GridlessADPFactoryTest();
		tests.parseArgs(args);
		tests.test_fromGrid_3BusSystem();
		tests.test_fromGrid_3BusSystemWith2OfEach();
	}

	private void test_fromGrid_3BusSystemWith2OfEach()
	{
		startTest();
		
		Grid grid = makeGrid();
		
		DistributedSource dg = new DistributedSource();
		dg.setName("DG2");
		ControllableDemand storage = new ControllableDemand();
		storage.setName("S2");
		Load load = new Load();
		load.setName("L2");

		setupUnits(load, storage, dg);
		
		grid.getBus("2").addChild(dg);
		grid.getBus("3").addChild(storage);
		grid.getBus("4").addChild(load);
		
		GridlessADP adp = GridlessADPFactory.fromGrid(grid, BASE_VOLTAGE, BASE_POWER, /*failIfNotConvergent=*/true);

		assertTrue(adp.getGridlessData().dgCount == 2, "fromGrid(grid).dgCount incorrect, is "+adp.getGridlessData().dgCount+", should be 2");
		assertTrue(adp.getGridlessData().loadCount == 2, "fromGrid(grid).loadCount incorrect, is "+adp.getGridlessData().loadCount+", should be 2");
		assertTrue(adp.getGridlessData().storageCount == 2, "fromGrid(grid).storageCount incorrect, is "+adp.getGridlessData().storageCount+", should be 2");
		
		endTest("test_fromGrid_3BusSystemWith2OfEach()");
	}

	private void test_fromGrid_3BusSystem()
	{
		startTest();
		
		Grid grid = makeGrid();
		
		GridlessADP adp = GridlessADPFactory.fromGrid(grid, BASE_VOLTAGE, BASE_POWER, /*failIfNotConvergent=*/true);
		
		assertTrue(adp != null, "fromGrid(grid) returned null");
		assertTrue(adp.getGridlessData().sensitivities != null, "fromGrid(grid).sensitivities is null");
		assertTrue(adp.getGridlessData().v_0 != null, "fromGrid(grid).v_0 is null");
		assertTrue(adp.getX0() != null, "fromGrid(grid).x_0 is null");
		
		assertTrue(adp.getGridlessData().dgCount == 1, "fromGrid(grid).dgCount incorrect, is "+adp.getGridlessData().dgCount+", should be 1");
		assertTrue(adp.getGridlessData().loadCount == 1, "fromGrid(grid).loadCount incorrect, is "+adp.getGridlessData().loadCount+", should be 1");
		assertTrue(adp.getGridlessData().storageCount == 1, "fromGrid(grid).storageCount incorrect, is "+adp.getGridlessData().storageCount+", should be 1");
		
		// Check storage parameters:
		assertTrue(adp.getGridlessData().storageChargeRate[0] == S_P_0, "Check storage parameters: Incorrect storage charge rate of "+adp.getGridlessData().storageChargeRate[0]+", should be "+S_P_0);
		assertTrue(adp.getGridlessData().storageDischargeRate[0] == S_MAX_DISCHARGE, "Check storage parameters: Incorrect storage max discharge rate of "+adp.getGridlessData().storageDischargeRate[0]+", should be "+S_MAX_DISCHARGE);
		assertTrue(adp.getGridlessData().storageChargeRate[0] == S_MAX_CHARGE, "Check storage parameters: Incorrect storage max charge rate of "+adp.getGridlessData().storageChargeRate[0]+", should be "+S_MAX_CHARGE);
		
		// Check correct x_0:
		assertTrue(adp.getX0().getDimension() == 1, "Check correct x_0: x_0 has incorrect dimension of "+adp.getX0().getDimension()+", should be 1");
		
		double soc_0 = adp.getX0().getEntry(0);
		double dgP_0 = adp.getGridlessData().power_0.getEntry(0);
		double storageP_0 = adp.getGridlessData().power_0.getEntry(1);
		double loadP_0 = adp.getGridlessData().power_0.getEntry(2);
		double dgQ_0 = adp.getGridlessData().power_0.getEntry(3);
		double storageQ_0 = adp.getGridlessData().power_0.getEntry(4);
		double loadQ_0 = adp.getGridlessData().power_0.getEntry(5);
		
		assertTrue(soc_0 == SOC_0, "Check correct x_0: Incorrect SOC_0 of "+soc_0+", should be "+SOC_0);
		assertTrue(dgP_0 == DG_P_0, "Check correct x_0: Incorrect DG P_0 of "+dgP_0+", should be "+DG_P_0);
		assertTrue(dgQ_0 == 0, "Check correct x_0: Incorrect DG Q_0 of "+dgQ_0+", should be "+0);
		assertTrue(storageP_0 == S_P_0, "Check correct x_0: Incorrect storage P_0 of "+storageP_0+", should be "+S_P_0);
		assertTrue(storageQ_0 == 0, "Check correct x_0: Incorrect storage Q_0 of "+storageQ_0+", should be "+0);
		assertTrue(loadP_0 == L_P_0, "Check correct x_0: Incorrect load P_0 of "+loadP_0+", should be "+L_P_0);
		assertTrue(loadQ_0 == L_Q_0, "Check correct x_0: Incorrect load Q_0 of "+loadQ_0+", should be "+L_Q_0);
		
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
		
		// Test correct sensitivity matrix:
		double SENSITIVITIES[][] = new double[][] {
				{0.0615049668,	0.0303326346,	0.0302831268,	-0.1254917978,	-0.062214813,	-0.0622315891},
				{0.0311522688,	0.0311585708,	0.0311749326,	-0.0618049959,	-0.0618050356,	-0.0617994913},
				{0.0318137905,	0.0318410568,	0.0653527807,	-0.0614662572,	-0.0614664292,	-0.127441907},
				{0.1216990136,	0.063452088,	0.0685124728,	0.0620673714,	0.0314809161,	0.0331956628},
				{0.0622655962,	0.064238859,	0.0693619897,	0.0318837095,	0.0318712621,	0.0336072708},
				{0.0645928275,	0.0666398427,	0.1386513019,	0.0330753911,	0.0330624785,	0.0677990438},

		};

		int nCols = adp.getGridlessData().sensitivities.getColumnDimension();
		int nRows = adp.getGridlessData().sensitivities.getRowDimension();
		assertTrue(nCols == 6, "Test sensitivity matrix order: Wrong number of columns "+nCols+", should be 6");
		assertTrue(nRows == 6, "Test sensitivity matrix order: Wrong number of rows "+nCols+", should be 6");
		
		if(verbose)
		{
			System.out.println("Converted Sensitivities:");
			for(int i = 0; i < nRows; ++i)
			{
				System.out.print("{");
				for(int j = 0; j < nCols; ++j)
				{
					double value = adp.getGridlessData().sensitivities.getEntry(i, j);
					System.out.print(value);
					System.out.print(",");
				}
				System.out.println("},");
			}
		}
			
		for(int i = 0; i < nRows; ++i)
		{
			for(int j = 0; j < nCols; ++j)
			{
				double v_ij = adp.getGridlessData().sensitivities.getEntry(i, j);
				assertTrue(equals(v_ij, SENSITIVITIES[i][j], 1e-6), "Test correct sensitivity matrix: Incorrect value of "+v_ij+" at row/col="+i+"/"+j+", should be "+SENSITIVITIES[i][j]);
			}
		}
		
		// Check correct v_0:
		assertTrue(adp.getGridlessData().v_0.getDimension() == 6, "Check correct v_0: v_0 has incorrect dimension of "+adp.getGridlessData().v_0.getDimension()+", should be 6");

		final double V_ABS_BUS2 = 0.9814111551;
		final double V_ARG_BUS2 = 0.00321262;
		final double V_ABS_BUS3 = 0.9692031462;
		final double V_ARG_BUS3 = -0.0030953313;
		final double V_ABS_BUS4 = 0.9355200999;
		final double V_ARG_BUS4 = -0.0080583551;
		
		final double v_arg_bus2 = adp.getGridlessData().v_0.getEntry(0);
		final double v_arg_bus3 = adp.getGridlessData().v_0.getEntry(1);
		final double v_arg_bus4 = adp.getGridlessData().v_0.getEntry(2);
		final double v_abs_bus2 = adp.getGridlessData().v_0.getEntry(3);
		final double v_abs_bus3 = adp.getGridlessData().v_0.getEntry(4);
		final double v_abs_bus4 = adp.getGridlessData().v_0.getEntry(5);
		
		assertTrue(equals(v_arg_bus2, V_ARG_BUS2, 1e-6), "Check correct v_0: Bus 2 voltage has incorerct arg value of "+v_arg_bus2+", should be "+V_ARG_BUS2);
		assertTrue(equals(v_arg_bus3, V_ARG_BUS3, 1e-6), "Check correct v_0: Bus 3 voltage has incorerct arg value of "+v_arg_bus3+", should be "+V_ARG_BUS3);
		assertTrue(equals(v_arg_bus4, V_ARG_BUS4, 1e-6), "Check correct v_0: Bus 4 voltage has incorerct arg value of "+v_arg_bus4+", should be "+V_ARG_BUS4);
		assertTrue(equals(v_abs_bus2, V_ABS_BUS2, 1e-6), "Check correct v_0: Bus 2 voltage has incorerct abs value of "+v_abs_bus2+", should be "+V_ABS_BUS2);
		assertTrue(equals(v_abs_bus3, V_ABS_BUS3, 1e-6), "Check correct v_0: Bus 3 voltage has incorerct abs value of "+v_abs_bus3+", should be "+V_ABS_BUS3);
		assertTrue(equals(v_abs_bus4, V_ABS_BUS4, 1e-6), "Check correct v_0: Bus 4 voltage has incorerct abs value of "+v_abs_bus4+", should be "+V_ABS_BUS4);
		
		endTest("test_fromGrid_3BusSystem()");
	}

	private Grid makeGrid()
	{
		// Grid with 3 busses in this order: Storage, Load, DG:
		Grid grid = Grid.grid().
				Bus("1").
				SlackSource("Slack", BASE_VOLTAGE, 0, 0, 0).
				Line("1-2", 1, 0.06*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
					Bus("2").
						ControllableDemand("S").
						Line("2-3", 1, 0.06*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
							Bus("3").
								Load("L").
							terminate().
						terminate().
						Line("2-4", 1, 0.06*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
						Bus("4").
							DistributedSource("DG").
						terminate().
					terminate().
					terminate().
				terminate().
			terminate().
		grid();
		
		// x_0:
		Load load = grid.getLoad("L");
		ControllableDemand storage = grid.getControllableDemand("S");
		DistributedSource dg = grid.getDistributedSource("DG");
		setupUnits(load, storage, dg);
		
		return grid;
	}

	private void setupUnits(Load load, ControllableDemand storage,
			DistributedSource dg)
	{
		load.setLoad(-L_P_0*BASE_POWER, -L_Q_0*BASE_POWER);
		
		storage.setMaxCapacity(S_MAX*BASE_POWER);
		storage.setCapacity(SOC_0*BASE_POWER);
		storage.setMaxChargeRate(-S_MAX_CHARGE*BASE_POWER);
		storage.setMaxDischargeRate(S_MAX_DISCHARGE*BASE_POWER);
		storage.setChargeRate(-S_P_0*BASE_POWER, 0);
		
		dg.setPmax(DG_MAX*BASE_POWER);
		dg.setPowerOutput(DG_P_0*BASE_POWER, 0);
	}
}
