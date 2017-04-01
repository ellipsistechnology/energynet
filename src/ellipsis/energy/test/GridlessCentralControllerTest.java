package ellipsis.energy.test;

import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeGridlessADPWith2DG2S;
import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeGridlessADPWith6DG2S;
import static ellipsis.util.VectorHelper.printVector;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Timer;

import ellipsis.energy.util.DiscreteGridlessADP;
import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.test.TestHelper;
import ellipsis.util.Pair;

public class GridlessCentralControllerTest extends GridTestHelper
{
	public static final RealVector S_max = new ArrayRealVector(new double[]{
			 // 6 DGs                   2 Ss    8 Loads                                                 6 DGs                   2 Ss    8 Loads
				0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0,-0.6,-0.255,-0.345,-0.192,-1.7325,-1.2645,-0.255,-0.255,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.435,-0.1875,-0.198,-0.129,-0.99,-0.693,-0.2265,-0.12});

	public static void main(String[] args) throws IOException
	{
		assertAssertsOn();
		
		GridlessCentralControllerTest tests = new GridlessCentralControllerTest();
		tests.parseArgs(args);
		
		tests.test_maxForecastPower();
		tests.test_ebeLambdaDeltaSProduct();
		
		tests.test_selectSubsets(1, 5);
		tests.test_selectSubsets(4, 10);
		tests.test_fillInSensitivitiesAndCounts(1, 5);
		tests.test_fillInSensitivitiesAndCounts(4, 10);
		tests.test_fillInForecastAndInitialState(1, 5);
		tests.test_fillInForecastAndInitialState(4, 10);
		tests.test_makeGridlessADP(1, 5);
		tests.test_makeGridlessADP(4, 10);
		
		tests.test_decompose(1, 5, 4);
		tests.test_decompose(4, 10, 4);
		tests.test_decompose(4, 10, 8);
		
		tests.test_epsilon(4, 10, 4);
		tests.test_epsilon(4, 10, 8);
		
		tests.test_optimise(4, 10, 4);
		tests.test_optimise(4, 10, 8);
	
		tests.test_externalVoltageChanges(4, 10, 4);
		tests.test_externalVoltageChanges(4, 10, 8);
		
		tests.test_iteration(4, 10, 4);
		tests.test_iteration(4, 10, 8);

		tests.test_loop(4, 10, 4, 1, false); // No external voltage updates.
		tests.test_loop(4, 14, 4, 1, true);  // Full network case (external voltages should be zero since nothing is external).
		tests.test_loop(4, 10, 4, 2, false); // No external voltage updates and compare two optimisations.
		tests.test_loop(4, 10, 4, 1, true);  // Small subsets.
		tests.test_loop(4, 10, 8, 1, true);  // Large subsets.
		
		System.out.println("\nPress enter to continue...");
		System.in.read();
	}

	private void test_fillInSensitivitiesAndCounts(int C_max, int U_max)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.C_max = C_max;
		cc.U_max = U_max;
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector max = cc.maxForecastPowerChange();
		RealMatrix LambdaProduct = cc.ebeLambdaDeltaSProduct(max);
		
		List<List<Integer>> Bs = cc.selectSubsets(LambdaProduct);
		for (List<Integer> B : Bs)
		{
			DiscreteGridlessADP adp_B = new DiscreteGridlessADP();
			
			cc.fillInSensitivitiesAndCounts(adp_B, B);
			
			// Check sensitivities have been filled in:
			// (Testing if they are correct is done in test_makeGridlessADP)
			int rowCount = adp_B.getGridlessData().sensitivities.getRowDimension();
			int colCount = adp_B.getGridlessData().sensitivities.getColumnDimension();
			for(int i = 0; i < rowCount; ++i)
			{
				for(int j = 0; j < colCount; ++j)
				{
					assertTrue(adp_B.getGridlessData().sensitivities.getEntry(i, j) != 0, "test_fillInSensitivitiesAndCounts(): Zero sensitivity found in "+adp_B);
				}
			}
			
			// Check correct counts:
			int dgCount = countInRange(B, 0, cc.gridlessADP().getGridlessData().dgCount);
			int storageCount = countInRange(B, cc.gridlessADP().getGridlessData().dgCount, cc.gridlessADP().getGridlessData().dgCount+cc.gridlessADP().getGridlessData().storageCount);
			int loadCount = countInRange(B, cc.gridlessADP().getGridlessData().dgCount+cc.gridlessADP().getGridlessData().storageCount, cc.gridlessADP().getGridlessData().dgCount+cc.gridlessADP().getGridlessData().storageCount+cc.gridlessADP().getGridlessData().loadCount);

			assertTrue(adp_B.getGridlessData().dgCount == dgCount, "test_fillInSensitivitiesAndCounts(): DG count was "+adp_B.getGridlessData().dgCount+", should be "+dgCount);
			assertTrue(adp_B.getGridlessData().storageCount == storageCount, "test_fillInSensitivitiesAndCounts(): "+adp_B.getGridlessData().storageCount+", should be "+storageCount);
			assertTrue(adp_B.getGridlessData().loadCount == loadCount, "test_fillInSensitivitiesAndCounts(): "+adp_B.getGridlessData().loadCount+", should be "+loadCount);
		}
		
		endTest("test_fillInSensitivitiesAndCounts("+C_max+", "+U_max+")");
	}

	private void test_fillInForecastAndInitialState(int C_max, int U_max)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.C_max = C_max;
		cc.U_max = U_max;
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector max = cc.maxForecastPowerChange();
		RealMatrix LambdaProduct = cc.ebeLambdaDeltaSProduct(max);
		
		List<List<Integer>> Bs = cc.selectSubsets(LambdaProduct);
		for (List<Integer> B : Bs)
		{
			DiscreteGridlessADP adp_B = new DiscreteGridlessADP();
			
			cc.fillInSensitivitiesAndCounts(adp_B, B);
			cc.fillInForecastAndInitialState(adp_B, B);
			
			// Check that x_0, v_0, means and SDs are filled in:
			// (Testing if they are correct is done in test_makeGridlessADP)
			assertTrue(adp_B.getX0() != null, "test_fillInForecastAndInitialState(): x_0 was null");
			assertTrue(adp_B.getGridlessData().v_0 != null, "test_fillInForecastAndInitialState(): v_0 was null");
			assertTrue(adp_B.getGridlessData().meanP_DG != null, "test_fillInForecastAndInitialState(): meanP_DG was null");
			assertTrue(adp_B.getGridlessData().meanP_L != null, "test_fillInForecastAndInitialState(): meanP_L was null");
			assertTrue(adp_B.getGridlessData().meanQ_L != null, "test_fillInForecastAndInitialState(): meanQ_L was null");
			assertTrue(adp_B.getGridlessData().sdP_DG != null, "test_fillInForecastAndInitialState(): sdP_DG was null");
			assertTrue(adp_B.getGridlessData().sdP_L != null, "test_fillInForecastAndInitialState(): sdP_L was null");
			assertTrue(adp_B.getGridlessData().sdQ_L != null, "test_fillInForecastAndInitialState(): sdQ_L was null");
			
			int unitCount = adp_B.getGridlessData().dgCount+adp_B.getGridlessData().storageCount+adp_B.getGridlessData().loadCount;
			int x_0Length = adp_B.getGridlessData().storageCount;
			int v_0Length = unitCount*2;
			assertTrue(adp_B.getX0().getDimension() == x_0Length, "test_fillInForecastAndInitialState(): x_0 length was "+adp_B.getX0().getDimension()+", should be "+x_0Length);
			assertTrue(adp_B.getGridlessData().v_0.getDimension() == v_0Length, "test_fillInForecastAndInitialState(): v_0 length was "+adp_B.getGridlessData().v_0.getDimension()+", should be "+v_0Length);
			
			assertTrue(adp_B.getGridlessData().meanP_DG.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but meanP_DG length was "+adp_B.getGridlessData().meanP_DG.length);
			assertTrue(adp_B.getGridlessData().meanP_L.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but meanP_L length was "+adp_B.getGridlessData().meanP_L.length);
			assertTrue(adp_B.getGridlessData().meanQ_L.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but meanQ_L length was "+adp_B.getGridlessData().meanQ_L.length);
			assertTrue(adp_B.getGridlessData().sdP_DG.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but sdP_DG length was "+adp_B.getGridlessData().sdP_DG.length);
			assertTrue(adp_B.getGridlessData().sdP_L.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but sdP_L length was "+adp_B.getGridlessData().sdP_L.length);
			assertTrue(adp_B.getGridlessData().sdQ_L.length == cc.localControllerTemplate.getT(), "test_fillInForecastAndInitialState(): T = "+cc.localControllerTemplate.getT()+", but sdQ_L length was "+adp_B.getGridlessData().sdQ_L.length);
			
			for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
			{
				assertTrue(adp_B.getGridlessData().meanP_DG[t].length == adp_B.getGridlessData().dgCount, "test_fillInForecastAndInitialState(): dgCount = "+adp_B.getGridlessData().dgCount+", but meanP_DG["+t+"] length was "+adp_B.getGridlessData().meanP_DG[t].length);
				assertTrue(adp_B.getGridlessData().meanP_L[t].length == adp_B.getGridlessData().loadCount, "test_fillInForecastAndInitialState(): loadCount = "+adp_B.getGridlessData().loadCount+", but meanP_L["+t+"] length was "+adp_B.getGridlessData().meanP_L[t].length);
				assertTrue(adp_B.getGridlessData().meanQ_L[t].length == adp_B.getGridlessData().loadCount, "test_fillInForecastAndInitialState(): loadCount = "+adp_B.getGridlessData().loadCount+", but meanQ_L["+t+"] length was "+adp_B.getGridlessData().meanQ_L[t].length);
				assertTrue(adp_B.getGridlessData().sdP_DG[t].length == adp_B.getGridlessData().dgCount, "test_fillInForecastAndInitialState(): dgCount = "+adp_B.getGridlessData().dgCount+", but sdP_DG["+t+"] length was "+adp_B.getGridlessData().sdP_DG[t].length);
				assertTrue(adp_B.getGridlessData().sdP_L[t].length == adp_B.getGridlessData().loadCount, "test_fillInForecastAndInitialState(): loadCount = "+adp_B.getGridlessData().loadCount+", but sdP_L["+t+"] length was "+adp_B.getGridlessData().sdP_L[t].length);
				assertTrue(adp_B.getGridlessData().sdQ_L[t].length == adp_B.getGridlessData().loadCount, "test_fillInForecastAndInitialState(): loadCount = "+adp_B.getGridlessData().loadCount+", but sdQ_L["+t+"] length was "+adp_B.getGridlessData().sdQ_L[t].length);
				
				for(int i = 0; i < adp_B.getGridlessData().dgCount; ++i)
				{
					if(t < 9)
					{
						assertTrue(adp_B.getGridlessData().meanP_DG[t][i] != 0, "test_fillInForecastAndInitialState(): meanP_DG["+t+"]["+i+"] = 0");
						assertTrue(adp_B.getGridlessData().sdP_DG[t][i] != 0, "test_fillInForecastAndInitialState(): sdP_DG["+t+"]["+i+"] = 0");
					}
				}
				
				for(int i = 0; i < adp_B.getGridlessData().loadCount; ++i)
				{
					assertTrue(adp_B.getGridlessData().meanP_L[t][i] != 0, "test_fillInForecastAndInitialState(): meanP_L["+t+"]["+i+"] = 0");
					assertTrue(adp_B.getGridlessData().meanQ_L[t][i] != 0, "test_fillInForecastAndInitialState(): meanQ_L["+t+"]["+i+"] = 0");
					assertTrue(adp_B.getGridlessData().sdP_L[t][i] != 0, "test_fillInForecastAndInitialState(): sdP_L["+t+"]["+i+"] = 0");
					assertTrue(adp_B.getGridlessData().sdQ_L[t][i] != 0, "test_fillInForecastAndInitialState(): sdQ_L["+t+"]["+i+"] = 0");
				}
			}
		}
		
		endTest("test_fillInForecastAndInitialState("+C_max+", "+U_max+")");
	}

	private void test_maxForecastPower()
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector max = cc.maxForecastPower();
		int dimension = max.getDimension();
		assertTrue(S_max.getDimension() == dimension, "test_maxForecastPower(): Incorrect dimension for result of maxForecastPower(): "+dimension);
		for(int i = 0; i < dimension; ++i)
			assertTrue(equals(S_max.getEntry(i), max.getEntry(i), 1e-6), "test_maxForecastPower(): Incorrect max value at index "+i+" of "+max.getEntry(i)+", should be "+S_max.getEntry(i));
		
		endTest("test_maxForecastPower()");
	}

	private void test_ebeLambdaDeltaSProduct()
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector DeltaS_max = S_max.subtract(cc.S_0());
		
		RealVector max = cc.maxForecastPowerChange();
		RealMatrix sumElements = cc.ebeLambdaDeltaSProduct(max);
		int columnDimension = sumElements.getColumnDimension();
		int rowDimension = sumElements.getRowDimension();
		assertTrue(columnDimension == rowDimension, "test_ebeLambdaDeltaSProduct(): Column and row dimensions do not match");
		assertTrue(columnDimension == cc.gridlessADP().getGridlessData().sensitivities.getColumnDimension(), "test_ebeLambdaDeltaSProduct(): Incorrect column dimension of "+columnDimension+", should be "+cc.gridlessADP().getGridlessData().sensitivities.getColumnDimension());
		assertTrue(rowDimension == cc.gridlessADP().getGridlessData().sensitivities.getRowDimension(), "test_ebeLambdaDeltaSProduct(): Incorrect row dimension of "+rowDimension+", should be "+cc.gridlessADP().getGridlessData().sensitivities.getRowDimension());
		for(int i = 0; i < rowDimension; ++i)
		{
			for(int j = 0; j < columnDimension; ++j)
			{
				double Lambda_ij = sumElements.getEntry(i, j);
				double DeltaS = DeltaS_max.getEntry(j);
				double correctValue = cc.gridlessADP().getGridlessData().sensitivities.getEntry(i, j)*DeltaS;
				assertTrue(equals(Lambda_ij, correctValue, 1e-6), "test_ebeLambdaDeltaSProduct(): Incorrect value in sum component ("+i+","+j+"), was "+Lambda_ij+", should be "+correctValue);
			}
		}
		
		endTest("test_ebeLambdaDeltaSProduct()");
	}

	private void test_selectSubsets(int C_max, int U_max)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.C_max = C_max;
		cc.U_max = U_max;
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector max = cc.maxForecastPowerChange();
		RealMatrix LambdaProduct = cc.ebeLambdaDeltaSProduct(max);
		
		List<List<Integer>> Bs = cc.selectSubsets(LambdaProduct);
		assertTrue(Bs.size() == 8, "test_selectSubsets(): Incorrect number of subsets given, was "+Bs.size()+", should have been 8 (6 DGs and 2 storage)");
		Iterator<List<Integer>> iterator = Bs.iterator();
		int b = 0;
		
		// DG:
		for(int i = 0; i < 6; ++i)
		{
			List<Integer> B = iterator.next();
			
			// Check pilot:
			assertTrue(cc.pilots.get(b) == i, "test_selectSubsets(): Pilot for subset "+b+" was "+cc.pilots.get(b)+", should be "+i);
			assertTrue(B.contains(i), "test_selectSubsets(): Subset for DG"+i+" did not contain pilot "+i);
			
			checkSubset(cc, B, "DG"+i);
			
			++b;
		}
		
		// Storage:
		for(int i = 0; i < 2; ++i)
		{
			List<Integer> B = iterator.next();
			
			// Check pilot:
			assertTrue(cc.pilots.get(b) == i+6, "test_selectSubsets(): Pilot for subset "+b+" was "+cc.pilots.get(b)+", should be "+(i+6));
			assertTrue(B.contains(i+6), "test_selectSubsets(): Subset for DG"+i+" did not contain pilot "+(i+6));
			
			checkSubset(cc, B, "S"+i);

			++b;
		}
		
		endTest("test_selectSubsets("+C_max+", "+U_max+")");
	}

	public void checkSubset(GridlessCentralCoordinator cc, List<Integer> B, String subsetName)
	{
		int controllableCount = countControllables(B, cc.gridlessADP());
		assertTrue(controllableCount <= cc.C_max, "test_selectSubsets(): "+subsetName+": Too many controllable units: "+controllableCount);
		assertTrue(B.size() <= cc.U_max, "test_selectSubsets(): "+subsetName+": Too many units: "+B.size());
		assertTrue(controllableCount == cc.C_max && B.size() == cc.U_max, "test_selectSubsets(): "+subsetName+": Wrong size; unit count = "+B.size()+", controllable units = "+controllableCount);
	}

	private int countControllables(List<Integer> B, GridlessADP adp)
	{
		return countInRange(B, 0, adp.getGridlessData().dgCount+adp.getGridlessData().storageCount);
	}
	
	private int countInRange(List<Integer> B, int low, int high)
	{
		int count = 0;
		for (Integer b : B)
		{
			if(low <= b && b < high)
				++count;
		}
		
		return count;
	}

	private void test_makeGridlessADP(int C_max, int U_max)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.C_max = C_max;
		cc.U_max = U_max;
		cc.localControllerTemplate = makeGridlessADPWith6DG2S(false);
		
		RealVector max = cc.maxForecastPowerChange();
		RealMatrix LambdaProduct = cc.ebeLambdaDeltaSProduct(max);
		
		List<List<Integer>> Bs = cc.selectSubsets(LambdaProduct);

		cc.makeGridlessADP(Bs);
		printDecomposition(cc);
		assertTrue(cc.decomposition.length == Bs.size(), "test_makeGridlessADP(): Incorrect number of GridlessADPs, was "+cc.decomposition.length+", should be "+Bs.size());
		
		int dgCount = cc.gridlessADP().getGridlessData().dgCount;
		int storageCount = cc.gridlessADP().getGridlessData().storageCount;
		int loadCount = cc.gridlessADP().getGridlessData().loadCount;
		int unitCount = cc.gridlessADP().getGridlessData().dgCount+cc.gridlessADP().getGridlessData().storageCount+cc.gridlessADP().getGridlessData().loadCount;
		
		int b = 0;
		for (List<Integer> B : Bs)
		{
			int dgCount_B = cc.decomposition(b).getGridlessData().dgCount;
			int storageCount_B = cc.decomposition(b).getGridlessData().storageCount;
			int loadCount_B = cc.decomposition(b).getGridlessData().loadCount;
			
			int unitCount_B = dgCount_B+storageCount_B+loadCount_B;
			assertTrue(B.size() == unitCount_B, "test_makeGridlessADP(): Unit count was "+unitCount_B+", should be "+B.size());
			
			// Check counts and general parameters:
			{
				int dgCount_correct = countInRange(B, 0, dgCount);
				int storageCount_correct = countInRange(B, dgCount, dgCount+storageCount);
				int loadCount_correct = countInRange(B, dgCount+storageCount, dgCount+storageCount+loadCount);
				
				assertTrue(dgCount_B == dgCount_correct, "test_makeGridlessADP(): Incorrect number of DG for decomposition "+b+", was "+dgCount_B+", should be "+dgCount_correct);
				assertTrue(storageCount_B == storageCount_correct, "test_makeGridlessADP(): Incorrect number of storage for decomposition "+b+", was "+storageCount_B+", should be "+storageCount_correct);
				assertTrue(loadCount_B == loadCount_correct, "test_makeGridlessADP(): Incorrect number of loads for decomposition "+b+", was "+loadCount_B+", should be "+loadCount_correct);
				
				assertTrue(cc.localControllerTemplate.getT() == cc.decomposition(b).getT(), "test_makeGridlessADP(): Decomposed GridlessADP.T is "+cc.decomposition(b).getT()+" for decomposition "+b+", should be "+cc.localControllerTemplate.getT());
				assertTrue(cc.gridlessADP().getNr() == cc.decomposition(b).getNr(), "test_makeGridlessADP(): Decomposed GridlessADP.N_r is "+cc.decomposition(b).getNr()+" for decomposition "+b+", should be "+cc.gridlessADP().getNr());
				assertTrue(cc.gridlessADP().getGamma() == cc.decomposition(b).getGamma(), "test_makeGridlessADP(): Decomposed GridlessADP.gamma is "+cc.decomposition(b).getGamma()+" for decomposition "+b+", should be "+cc.gridlessADP().getGamma());
			}
			
			// Check storage parameters:
			{
				Iterator<Integer> iterator = B.iterator();
				for(int i = 0; i < unitCount_B; ++i) // Need to iterate over all units to get correct original index.
				{
					int i_original = iterator.next();
					if(dgCount_B <= i && i < storageCount_B) // Is it storage?
					{
						assertTrue(cc.decomposition(b).getGridlessData().storageChargeRate[i_original-dgCount] == cc.gridlessADP().getGridlessData().storageChargeRate[i-dgCount_B], "test_makeGridlessADP(): Storage charge rate was "+cc.decomposition(b).getGridlessData().storageChargeRate[i-dgCount_B]+", should be "+cc.gridlessADP().getGridlessData().storageChargeRate[i_original-dgCount]);
						assertTrue(cc.decomposition(b).getGridlessData().storageDischargeRate[i_original-dgCount] == cc.gridlessADP().getGridlessData().storageDischargeRate[i-dgCount_B], "test_makeGridlessADP(): Storage discharge rate was "+cc.decomposition(b).getGridlessData().storageDischargeRate[i-dgCount_B]+", should be "+cc.gridlessADP().getGridlessData().storageDischargeRate[i_original-dgCount]);
						assertTrue(cc.decomposition(b).getGridlessData().storageMaxCapacity[i_original-dgCount] == cc.gridlessADP().getGridlessData().storageMaxCapacity[i-dgCount_B], "test_makeGridlessADP(): Storage max capacity was "+cc.decomposition(b).getGridlessData().storageMaxCapacity[i-dgCount_B]+", should be "+cc.gridlessADP().getGridlessData().storageMaxCapacity[i_original-dgCount]);
					}
				}
			}
			
			// Check forecasts:
			{
				int loadOffset_B = dgCount_B+storageCount_B;
				int loadOffset = dgCount+storageCount;
				for(int t = 0; t < cc.decomposition(b).getT(); ++t)
				{
					Iterator<Integer> iterator = B.iterator();
					for(int i = 0; i < unitCount_B; ++i)
					{
						int i_original = iterator.next();
						
						if(i < dgCount_B)
						{
							assertTrue(cc.decomposition(b).getGridlessData().meanP_DG[t][i] == cc.gridlessADP().getGridlessData().meanP_DG[t][i_original], "test_makeGridlessADP(): Mean of P_DG was "+cc.decomposition(b).getGridlessData().meanP_DG[t][i]+", should be"+cc.gridlessADP().getGridlessData().meanP_DG[t][i_original]);
							assertTrue(cc.decomposition(b).getGridlessData().sdP_DG[t][i] == cc.gridlessADP().getGridlessData().sdP_DG[t][i_original], "test_makeGridlessADP(): SD of P_DG was "+cc.decomposition(b).getGridlessData().sdP_DG[t][i]+", should be"+cc.gridlessADP().getGridlessData().sdP_DG[t][i_original]);
						}
						else if(i >= loadOffset_B)
						{
							assertTrue(cc.decomposition(b).getGridlessData().meanP_L[t][i-loadOffset_B] == cc.gridlessADP().getGridlessData().meanP_L[t][i_original-loadOffset], "test_makeGridlessADP(): Mean of P_L was "+cc.decomposition(b).getGridlessData().meanP_L[t][i-loadOffset_B]+", should be"+cc.gridlessADP().getGridlessData().meanP_L[t][i_original-loadOffset]);
							assertTrue(cc.decomposition(b).getGridlessData().sdP_L[t][i-loadOffset_B] == cc.gridlessADP().getGridlessData().sdP_L[t][i_original-loadOffset], "test_makeGridlessADP(): SD of P_L was "+cc.decomposition(b).getGridlessData().sdP_L[t][i-loadOffset_B]+", should be"+cc.gridlessADP().getGridlessData().sdP_L[t][i_original-loadOffset]);
							assertTrue(cc.decomposition(b).getGridlessData().meanQ_L[t][i-loadOffset_B] == cc.gridlessADP().getGridlessData().meanQ_L[t][i_original-loadOffset], "test_makeGridlessADP(): Mean of Q_L was "+cc.decomposition(b).getGridlessData().meanQ_L[t][i-loadOffset_B]+", should be"+cc.gridlessADP().getGridlessData().meanQ_L[t][i_original-loadOffset]);
							assertTrue(cc.decomposition(b).getGridlessData().sdQ_L[t][i-loadOffset_B] == cc.gridlessADP().getGridlessData().sdQ_L[t][i_original-loadOffset], "test_makeGridlessADP(): SD of Q_L was "+cc.decomposition(b).getGridlessData().sdQ_L[t][i-loadOffset_B]+", should be"+cc.gridlessADP().getGridlessData().sdQ_L[t][i_original-loadOffset]);
						}
					}
				}
			}
			
			// Check sensitivities:
			{
				assertTrue(cc.decomposition(b).getGridlessData().sensitivities.getColumnDimension() == 2*unitCount_B, "test_makeGridlessADP(): Incorrect column dimension of "+cc.decomposition(b).getGridlessData().sensitivities.getColumnDimension()+", should be "+2*unitCount_B);
				assertTrue(cc.decomposition(b).getGridlessData().sensitivities.getRowDimension() == 2*unitCount_B, "test_makeGridlessADP(): Incorrect row dimension of "+cc.decomposition(b).getGridlessData().sensitivities.getRowDimension()+", should be "+2*unitCount_B);
				Iterator<Integer> iterator_i = B.iterator();
				for(int i = 0; i < unitCount_B; ++i)
				{
					int i_original = iterator_i.next();
					Iterator<Integer> iterator_j = B.iterator();
					for(int j = 0; j < unitCount_B; ++j)
					{
						int j_original = iterator_j.next();
						
						// Top-left delta/P:
						double value_delta_P = cc.decomposition(b).getGridlessData().sensitivities.getEntry(i, j);
						double correctValue_delta_P = cc.gridlessADP().getGridlessData().sensitivities.getEntry(i_original, j_original);
						assertTrue(value_delta_P == correctValue_delta_P, "test_makeGridlessADP(): Wrong sensitivity at ("+i+","+j+"), was "+value_delta_P+", should be "+correctValue_delta_P+" at original coordinate ("+i_original+","+j_original+")");

						// Top-right delta/Q:
						double value_delta_Q = cc.decomposition(b).getGridlessData().sensitivities.getEntry(i, j+unitCount_B);
						double correctValue_delta_Q = cc.gridlessADP().getGridlessData().sensitivities.getEntry(i_original, j_original+unitCount);
						assertTrue(value_delta_Q == correctValue_delta_Q, "test_makeGridlessADP(): Wrong sensitivity at ("+i+","+j+unitCount_B+"), was "+value_delta_Q+", should be "+correctValue_delta_Q+" at original coordinate ("+i_original+","+j_original+unitCount+")");

						// Bottom-left abs/P:
						double value_abs_P = cc.decomposition(b).getGridlessData().sensitivities.getEntry(i+unitCount_B, j);
						double correctValue_abs_P = cc.gridlessADP().getGridlessData().sensitivities.getEntry(i_original+unitCount, j_original);
						assertTrue(value_abs_P == correctValue_abs_P, "test_makeGridlessADP(): Wrong sensitivity at ("+i+unitCount_B+","+j+"), was "+value_abs_P+", should be "+correctValue_abs_P+" at original coordinate ("+i_original+unitCount+","+j_original+")");

						// Bottom-right abs/Q:
						double value_abs_Q = cc.decomposition(b).getGridlessData().sensitivities.getEntry(i+unitCount_B, j+unitCount_B);
						double correctValue_abs_Q = cc.gridlessADP().getGridlessData().sensitivities.getEntry(i_original+unitCount, j_original+unitCount);
						assertTrue(value_abs_Q == correctValue_abs_Q, "test_makeGridlessADP(): Wrong sensitivity at ("+i+unitCount_B+","+j+unitCount_B+"), was "+value_abs_Q+", should be "+correctValue_abs_Q+" at original coordinate ("+i_original+unitCount+","+j_original+unitCount+")");
					}
				}
			}
			
			// Initial voltages:
			{
				assertTrue(cc.decomposition(b).getGridlessData().v_0.getDimension() == unitCount_B*2, "test_makeGridlessADP(): v_0 has dimension "+cc.decomposition(b).getGridlessData().v_0.getDimension()+", but unit count was "+unitCount_B);
				Iterator<Integer> iterator_i = B.iterator();
				for(int i = 0; i < unitCount_B; ++i)
				{
					int i_original = iterator_i.next();
					assertTrue(cc.decomposition(b).getGridlessData().v_0.getEntry(i) == cc.gridlessADP().getGridlessData().v_0.getEntry(i_original), "test_makeGridlessADP(): v_0["+i+"]_B = "+cc.decomposition(b).getGridlessData().v_0.getEntry(i)+" != "+cc.gridlessADP().getGridlessData().v_0.getEntry(i_original)+" = v_0["+i_original+"]");
				}
			}
			
			// External voltages:
			{
				assertTrue(cc.decomposition(b).getGridlessData().deltaV_external.length == cc.localControllerTemplate.getT()+1, "test_makeGridlessADP(): deltaV_external has length "+cc.decomposition(b).getGridlessData().deltaV_external.length+", but T = "+cc.localControllerTemplate.getT());
				for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
				{
					assertTrue(cc.decomposition(b).getGridlessData().deltaV_external[t].getDimension() == unitCount_B*2, "test_makeGridlessADP(): deltaV_external["+t+"] has dimension "+cc.decomposition(b).getGridlessData().deltaV_external[t].getDimension()+", but unit count was "+unitCount_B);
					for(int i = 0; i < unitCount_B; ++i)
					{
						assertTrue(cc.decomposition(b).getGridlessData().deltaV_external[t].getEntry(i) == 0, "test_makeGridlessADP(): deltaV_external["+t+"] = "+cc.decomposition(b).getGridlessData().deltaV_external[t].getEntry(i)+", should be 0");
					}
				}
			}
			
			// Initial state:
			{
				assertTrue(cc.decomposition(b).getX0().getDimension() == storageCount_B, "test_makeGridlessADP(): x_0 has dimension "+cc.decomposition(b).getX0().getDimension()+", should be "+(2*unitCount_B+storageCount_B));
				Iterator<Integer> iterator_i = B.iterator();
				for(int i = 0; i < dgCount_B+storageCount_B; ++i)
				{
					int i_original = iterator_i.next();
					if(i < dgCount_B) // Only check SOC elements in x_0; unit order is DG, Storage, Loads.
						continue;
					
					double x_b_0_i = cc.decomposition(b).getX0().getEntry(i-dgCount_B);
					double x_0_i = cc.gridlessADP().getX0().getEntry(i_original-dgCount);
					assertTrue(x_b_0_i == x_0_i, "test_makeGridlessADP(): x_0["+i+"]_B = "+x_b_0_i+" != "+x_0_i+" = x_0["+i_original+"]");
//					assertTrue(cc.decomposition(b).x_0.getEntry(storageCount_B+unitCount_B+i) == cc.adp().x_0.getEntry(storageCount+unitCount+i_original), "test_makeGridlessADP(): x_0["+(storageCount_B+unitCount_B+i)+"]_B = "+cc.decomposition(b).x_0.getEntry(storageCount_B+unitCount_B+i)+" != "+cc.adp().x_0.getEntry(storageCount+unitCount+i_original)+" = x_0["+storageCount+unitCount+i_original+"]");
//					if(dgCount_B <= i && i < dgCount_B+storageCount_B)
//					{
//						int s = i - dgCount_B;
//						int s_original = i_original - dgCount;
//						assertTrue(cc.decomposition(b).x_0.getEntry(s) == cc.adp().x_0.getEntry(s_original), "test_makeGridlessADP(): Incorrect SOC state x_0["+s+"]_B = "+cc.decomposition(b).x_0.getEntry(s)+", should be x_0["+s_original+"] = "+cc.adp().x_0.getEntry(s_original));
//					}
				}
			}
			
			// Initial power:
			{
				assertTrue(cc.decomposition(b).getGridlessData().power_0.getDimension() == 2*unitCount_B, "test_makeGridlessADP(): power_0 has dimension "+cc.decomposition(b).getGridlessData().power_0.getDimension()+", should be "+(2*unitCount_B));
				Iterator<Integer> iterator_i = B.iterator();
				for(int i = 0; i < unitCount_B; ++i)
				{
					int i_original = iterator_i.next();
					double power_b_0_i = cc.decomposition(b).getGridlessData().power_0.getEntry(i);
					double power_0_i = cc.gridlessADP().getGridlessData().power_0.getEntry(i_original);
					assertTrue(power_b_0_i == power_0_i, "test_makeGridlessADP(): power_0["+i+"]_B = "+power_b_0_i+" != "+power_0_i+" = power_0["+i_original+"]");
					assertTrue(cc.decomposition(b).getGridlessData().power_0.getEntry(unitCount_B+i) == cc.gridlessADP().getGridlessData().power_0.getEntry(unitCount+i_original), "test_makeGridlessADP(): power_0["+(unitCount_B+i)+"]_B = "+cc.decomposition(b).getGridlessData().power_0.getEntry(unitCount_B+i)+" != "+cc.gridlessADP().getGridlessData().power_0.getEntry(unitCount+i_original)+" = power_0["+(unitCount+i_original)+"]");
//					if(dgCount_B <= i && i < dgCount_B+storageCount_B)
//					{
//						int s = i - dgCount_B;
//						int s_original = i_original - dgCount;
//						assertTrue(cc.decomposition(b).x_0.getEntry(s) == cc.adp().x_0.getEntry(s_original), "test_makeGridlessADP(): Incorrect SOC state x_0["+s+"]_B = "+cc.decomposition(b).x_0.getEntry(s)+", should be x_0["+s_original+"] = "+cc.adp().x_0.getEntry(s_original));
//					}
				}
			}
			
			// Next subset:
			++b;
		} // End subset iterations.
		
		endTest("test_makeGridlessADP("+C_max+", "+U_max+")");
	}

	private void test_decompose(int C_max, int U_max, int controllableCount)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		
		cc.decompose();
		
		assertTrue(cc.decomposition != null, "test_decompose(): Decomposition was null");
		
		printDecomposition(cc);
		
		assertTrue(cc.decomposition.length == controllableCount, "test_decompose(): Wrong number of subsets: "+cc.decomposition.length+", should be "+controllableCount);
		for(int b = 0; b < cc.decomposition.length; ++b)
		{
			int controllableCount_B = cc.decomposition(b).getGridlessData().dgCount+cc.decomposition(b).getGridlessData().storageCount;
			int unitCount_B = controllableCount_B + cc.decomposition(b).getGridlessData().loadCount;
			assertTrue(controllableCount_B == C_max || unitCount_B == U_max, "test_decompose(): Wrong number of controllable units: was "+controllableCount_B+", should be "+C_max);
		}
		
		// TODO test specific subsets
		
		endTest("test_decompose("+C_max+", "+U_max+", "+controllableCount+")");
	}

	private void test_epsilon(int C_max, int U_max, int controllableCount)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		
		cc.decompose();
		
		assertTrue(cc.decomposition != null, "test_decompose(): Decomposition was null");
		
		assertTrue(cc.decomposition.length == controllableCount, "test_decompose(): Wrong number of subsets: "+cc.decomposition.length+", should be "+controllableCount);
		double[] epsilon = new double[cc.decomposition.length];
		for(int b = 0; b < cc.decomposition.length; ++b)
		{
			epsilon[b] = cc.epsilon(b);
			
			assertTrue(epsilon[b] >= 0, "test_epsilon(): epsilon = "+epsilon[b]+" < 0");
			
			if(verbose)
				System.out.println("epsilon["+b+"] = "+epsilon[b]);
		}

		// Test that epsilon increases with removal from B:
		double[] epsilon_increased = new double[cc.decomposition.length];
		RealVector DeltaS_t_max = cc.maxForecastPowerChange();
		RealMatrix LambdaProduct = cc.ebeLambdaDeltaSProduct(DeltaS_t_max);
		List<List<Integer>> Bs = cc.selectSubsets(LambdaProduct); // new Bs
		for(int b = 0; b < cc.decomposition.length; ++b) // remove one from each B
		{
			List<Integer> B = Bs.get(b);
			List<Integer> listB = new ArrayList<>(B);
			listB.remove(listB.size()-1);
			List<Integer> newB = new ArrayList<>(listB); // new B short last element
			Bs.set(b, newB);
		}
		cc.makeGridlessADP(Bs);
		for(int b = 0; b < cc.decomposition.length; ++b)
		{
			epsilon_increased[b] = cc.epsilon(b);
			assertTrue(epsilon[b] < epsilon_increased[b], "test_epsilon(): epsilon did not increase with removed node for decomposition "+b+": epsilon = "+epsilon[b]+", epsilon_increased = "+epsilon_increased[b]);
			
			if(verbose)
				System.out.println("epsilon_increased["+b+"] = "+epsilon_increased[b]);
		}
		
		endTest("test_epsilon("+C_max+", "+U_max+", "+controllableCount+")");
	}

	public void printDecomposition(GridlessCentralCoordinator controller)
	{
		if(verbose)
		{
			for(int i = 0; i < controller.decomposition.length; ++i)
			{
				System.out.println("Decomposition "+i+":");
				System.out.println("\tPilot = "+controller.pilots.get(i));
				System.out.println("\tDG Count = "+controller.decomposition(i).getGridlessData().dgCount);
				System.out.println("\tStorage Count = "+controller.decomposition(i).getGridlessData().storageCount);
				System.out.println("\tLoad Count = "+controller.decomposition(i).getGridlessData().loadCount);
			}
		}
	}

	private void test_optimise(int C_max, int U_max, int controllableCount)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		
		cc.decompose();
		
		double[][] bandwidth = cc.makeBandwidthArray(cc.gridlessADP().getGridlessData().storageCount, 0, 0, 0.5, 0, 0, 0);
		cc.optimise(0.5);
		
		if(verbose)
		{
			for(int i = 0; i < cc.schedule.length; ++i)
			{
				System.out.println("\nSchedule for ["+i+"]:");
				printSchedule(cc.schedule[i], cc.ctgs[i]);
			}
		}
		
		// Compare with deterministic calculation of subset:
		for(int i = 0; i < cc.schedule.length; ++i)
		{
			List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)cc.decomposition(i));
			for(int t = 0; t < schedule.size(); ++t)
			{
				Pair<RealVector, Double> pair = schedule.get(t);
				Double ctg = pair.getValue();
				assertTrue(equals(ctg, cc.ctgs[i][t], t==cc.localControllerTemplate.getT()?1:0.5), "Incorrect CTG at time "+t+" for subset "+i+", DP => "+ctg+", ADP => "+cc.ctgs[i][t]);
			}
			if(verbose)
			{
				printScheduleComparisonAndVoltages((DiscreteGridlessADP)cc.decomposition(i), schedule, cc.ctgs[i], cc.schedule[i]);
			}
		}
		
		// Combine subset controls and calculate full network CTG:
		RealVector[] combinedControl = cc.controlSchedule();
		double[] trueCtgForApproxSchedule = cc.gridlessADP().ctgs(combinedControl);
		assertTrue(trueCtgForApproxSchedule.length == cc.localControllerTemplate.getT()+1, "trueCtgForApproxSchedule.length = "+trueCtgForApproxSchedule.length+", but T = "+cc.localControllerTemplate.getT());
		
		// Non-divided problem (i.e. full network):
		if(verbose)
			System.out.println("\nFull network schedule:");
		Timer timer = Timer.startNewTimer();
//		double[] bandwidth = GridlessCentralController.makeBandwidthArray(bw, cc.adp().storageCount);
		cc.localControllerTemplate.train(bandwidth);
		double[] approxCtgs = new double[cc.localControllerTemplate.getT()+1];
		RealVector[] approxSchedule = cc.localControllerTemplate.schedule(approxCtgs);
		if(verbose)
		{
			System.out.println("Full network ADP ran in "+timer.stop());
			printSchedule(approxSchedule, approxCtgs);
		}
		
		// Deterministic comparison for full network:
		{
			if(verbose)
				System.out.println("\nDeterministic comparison with full network:");
			List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)cc.localControllerTemplate);
			for(int t = 0; t < schedule.size(); ++t)
			{
				Pair<RealVector, Double> pair = schedule.get(t);
				Double ctg = pair.getValue();
				assertTrue(equals(ctg, approxCtgs[t], t==cc.localControllerTemplate.getT()?1:0.5), "Incorrect CTG at time "+t+" for full network, DP CTG => "+ctg+", ADP CTG => "+approxCtgs[t]);
				assertTrue(equals(ctg, trueCtgForApproxSchedule[t], 0.5), "Incorrect CTG at time "+t+" for combined controls, DP CTG => "+ctg+", combined distributed ADP CTG => "+trueCtgForApproxSchedule[t]);
			}
			if(verbose)
			{
				printScheduleComparisonAndVoltages((DiscreteGridlessADP)cc.gridlessADP(), schedule, approxCtgs, approxSchedule);
				System.out.println("\nDeterministic comparison with combined control schedule:");
				printScheduleComparisonAndVoltages((DiscreteGridlessADP)cc.gridlessADP(), schedule, trueCtgForApproxSchedule, combinedControl);
			}
		}

		endTest("test_optimise("+C_max+","+U_max+","+controllableCount+")");
	}

	private void printSchedule(RealVector[] schedule, double[] ctgs)
	{
		for(int t = 0; t < schedule.length; ++t)
		{
			System.out.print(t);
			System.out.print(",");
			System.out.print(schedule[t]);
			System.out.print(",");
			System.out.println(ctgs[t]);
		}
	}

	private void test_externalVoltageChanges(int C_max, int U_max, int controllableCount)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		
		cc.decompose();
		
		cc.optimise(0.5);
		
		RealVector[] u = cc.controlSchedule();
		RealVector[][] externalVoltageChanges = cc.externalVoltageChanges(u); // [subset][time]
		
		// Check dimensions:
		for(int i = 0; i < externalVoltageChanges.length; ++i)
		{
			assertTrue(externalVoltageChanges[i] != null, "test_externalVoltageChanges(): externalVoltageChanges["+i+"] was null");
			assertTrue(externalVoltageChanges[i].length == cc.localControllerTemplate.getT()+1, "test_externalVoltageChanges(): Incorrect dimension for externalVoltageChange["+i+"], was "+externalVoltageChanges[i].length+", should be T+1 = "+(cc.localControllerTemplate.getT()+1));
			for(int t = 0; t < externalVoltageChanges[i].length; ++t)
			{
				assertTrue(externalVoltageChanges[i][t] != null, "test_externalVoltageChanges(): externalVoltageChanges["+i+"]["+t+"] was null");
				assertTrue(externalVoltageChanges[i][t].getDimension() == cc.decomposition(i).voltageDimension(), "test_externalVoltageChanges(): Incorrect dimension for external voltages; "+externalVoltageChanges[i][t].getDimension()+", should be "+cc.decomposition(i).voltageDimension());
			}
		}
		
		// Add each element to the approximate voltage changes 
		// and compare with the true voltage:
//		int i = 0;
//		for (List<Integer> B : cc.Bs)
		for(int i = 0; i < cc.decomposition.length; ++i)
		{
			RealVector x_B_t = cc.decomposition(i).getX0();
			RealVector x_t = cc.gridlessADP().getX0();
			for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
			{
				// Control and noise terms:
				RealVector u_B_t = new ArrayRealVector(cc.decomposition(i).controlDimension());
				{
					int j = 0;
					for (Integer b : cc.decomposition(i).getBusIndeces())
					{
						if(b >= cc.gridlessADP().getGridlessData().dgCount+cc.gridlessADP().getGridlessData().storageCount)
							break;
						double P = u[t].getEntry(b);
						u_B_t.setEntry(j, P);
						++j;
					}
				}
				RealVector wZero_B_t = cc.decomposition(i).wZero(t);
				RealVector wZero = cc.gridlessADP().wZero(t);
				
				// Get subset approximate v:
				RealVector approxV = cc.decomposition(i).voltagesFromControlAndNoise(t, u_B_t, wZero_B_t);
				
				// Calculate what should be true v (approx + error due to external changes):
				RealVector v = approxV.add(externalVoltageChanges[i][t]);
				
				// Get the true voltage from the full network (shrinking to subset vector):
				RealVector trueV = new ArrayRealVector(v.getDimension());
				{
					RealVector trueV_full = cc.gridlessADP().voltagesFromControlAndNoise(t, u[t], wZero);
					int j = 0;
					int absOffset = v.getDimension()/2;
					int absOffset_full = trueV_full.getDimension()/2;
					for (Integer b : cc.decomposition(i).getBusIndeces())
					{
						double arg = trueV_full.getEntry(b);
						double abs = trueV_full.getEntry(b+absOffset_full);
						trueV.setEntry(j, arg);
						trueV.setEntry(j+absOffset, abs);
						++j;
					}
				}
				
				// Find difference between calculated and true voltages (should be zero):
				RealVector diff = trueV.subtract(v);
				assertTrue(equals(diff.getNorm(), 0.0, 1e-6), "test_externalVoltageChanges(): Approximate v + external v did not equal true v for subset "+i+" at time "+t+";\n" +
						"\tapprox v =,"+printVector(approxV)+"\n" +
						"\texternal v =,"+printVector(externalVoltageChanges[i][t])+"\n" +
						"\tapprox v + external v =,"+printVector(v)+"\n" +
						"\ttrue v =,"+printVector(trueV));
				
				// Next state for subset:
				x_B_t = cc.decomposition(i).f(t, x_B_t, u_B_t, wZero_B_t);
				
				// Next state:
				x_t = cc.gridlessADP().f(t, x_t, u[t], wZero);
			}
			
			++i;
		}
		
		endTest("test_externalVoltageChanges("+C_max+","+U_max+","+controllableCount+")");
	}

	private void test_iteration(int C_max, int U_max, int controllableCount)
	{
		startTest();
		
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator();
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		
		cc.decompose();
		
		RealVector oldX_0 = new ArrayRealVector(cc.gridlessADP().getX0());
		RealVector oldV_0 = new ArrayRealVector(cc.gridlessADP().getGridlessData().v_0);
		
		cc.iteration(0.5);
		
		RealVector[] u = cc.controlSchedule();
		
		// Test global power_0 is updated (P_DGs and P_Ss from control vector):
		for(int i = 0; i < u[0].getDimension(); ++i)
		{
			assertTrue(cc.gridlessADP().getGridlessData().power_0.getEntry(i) == u[0].getEntry(i), "test_iteration(): Incorrect global x_0["+i+"] = "+cc.gridlessADP().getGridlessData().power_0.getEntry(i)+", should be "+u[0].getEntry(i));
		}

		// Test global v_0 is updated:
		// FIXME implement test global v_0 is updated
		int unimplemented;
//		RealVector newX_0 = cc.adp().x_0; // Reset to old v_0 and x_0 for calculating voltage.
//		RealVector newV_0 = cc.adp().v_0;
//		cc.adp().x_0 = oldX_0;
//		cc.adp().v_0 = oldV_0;
//		RealVector v_0 = cc.adp().voltages(0, newX_0);
//		cc.adp().x_0 = newX_0;
//		cc.adp().v_0 = newV_0;
//		for(int i = 0; i < v_0.getDimension(); ++i)
//		{
//			assertTrue(v_0.getEntry(i) == cc.adp().v_0.getEntry(i), "test_iteration(): Incorrect global v_0["+i+"] = "+cc.adp().v_0.getEntry(i)+", should be "+v_0.getEntry(i));
//		}

		int QOffset_full = cc.gridlessADP().unitCount();
		int absOffset_full = cc.gridlessADP().unitCount();
		for(int i = 0; i < cc.decomposition.length; ++i)
		{
			List<Integer> B = cc.decomposition[i].getBusIndeces();
			int j = 0;
			int QOffset_B = cc.decomposition(i).unitCount();
			int absOffset_B = cc.decomposition(i).unitCount();
			for (Integer b : B)
			{
				// Test local subset power_0 is updated:
				double P_b = cc.gridlessADP().getGridlessData().power_0.getEntry(b);
				double P_B_j = cc.decomposition(i).getGridlessData().power_0.getEntry(j);
				double Q_b = cc.gridlessADP().getGridlessData().power_0.getEntry(b+QOffset_full);
				double Q_B_j = cc.decomposition(i).getGridlessData().power_0.getEntry(j+QOffset_B);
				assertTrue(P_b==P_B_j, "test_iteration(): decomposition["+i+"].power_0["+(j)+"] = "+P_B_j+", should be global power_0["+(b)+"] = "+P_b);
				assertTrue(Q_b==Q_B_j, "test_iteration(): decomposition["+i+"].power_0["+(+j+QOffset_B)+"] = "+Q_B_j+", should be global power_0["+(b+QOffset_full)+"] = "+Q_b);
				
				// Test local subset v_0 is updated:
				double arg_b = cc.gridlessADP().getGridlessData().v_0.getEntry(b);
				double arg_B_j = cc.decomposition(i).getGridlessData().v_0.getEntry(j);
				double abs_b = cc.gridlessADP().getGridlessData().v_0.getEntry(b+QOffset_full);
				double abs_B_j = cc.decomposition(i).getGridlessData().v_0.getEntry(j+QOffset_B);
				assertTrue(arg_b==arg_B_j, "test_iteration(): decomposition["+i+"].v_0["+(j)+"] = "+arg_B_j+", should be global v_0["+(b)+"] = "+arg_b);
				assertTrue(abs_b==abs_B_j, "test_iteration(): decomposition["+i+"].v_0["+(j+absOffset_B)+"] = "+abs_B_j+", should be global v_0["+(b+absOffset_full)+"] = "+abs_b);
				
				// Test external voltages are updated:
				RealVector[][] externalVoltageChanges = cc.externalVoltageChanges(u);
				for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
				{
					assertTrue(cc.decomposition(i).getGridlessData().deltaV_external[t].equals(externalVoltageChanges[i][t]), "test_iteration(): decomposition["+i+"].deltaV_external["+t+"] = "+cc.decomposition(i).getGridlessData().deltaV_external[t]+", should be "+externalVoltageChanges[i][t]);
				}
				
				++j;
			}
		}
		
		endTest("test_iteration()");
	}

	private void test_loop(
			int C_max, int U_max, int controllableCount, 
			int executions, 
			final boolean withExternalVoltageUpdates)
	{
		startTest();
		
		int iterations = 5;
		double alpha = 1;//0.2;
		
		// TODO compare the results of different alphas on the following test
		
		// Create and setup central controller:
		// Disable the external voltage updates by forcing them all to zero:
		GridlessCentralCoordinator cc = new GridlessCentralCoordinator(){
			@Override
			public RealVector[][] externalVoltageChanges(RealVector[] u)
			{
				RealVector[][] externalVoltageChanges = super.externalVoltageChanges(u);
				if(!withExternalVoltageUpdates)
				{
					for (int i = 0; i < externalVoltageChanges.length; i++)
					{
						for(int t = 0; t < externalVoltageChanges[i].length; ++t)
						{
							externalVoltageChanges[i][t].set(0);
						}
					}
				}
				return externalVoltageChanges;
			}
		};
		
		cc.localControllerTemplate = controllableCount == 4 ? makeGridlessADPWith2DG2S() : makeGridlessADPWith6DG2S(false);
		cc.C_max = C_max;
		cc.U_max = U_max;
		cc.alpha = alpha;

		cc.decompose();
		
		// Callback for logging and assertions after each loop:
		// Note: This must be after decompose() is called since its initialisation
		//       performs a dp() on each subset.
		test_loop_LoopFinishedCallback callback = new test_loop_LoopFinishedCallback(cc);
		callback.withExternalVoltageUpdates = withExternalVoltageUpdates;
		cc.addLoopFinishedCallback(callback);
		
		// Optimisation loop:
		RealVector x_0 = new ArrayRealVector(cc.gridlessADP().getX0());
		RealVector[] schedule = null;
		double[][] ctgs = null;
		for(int i = 0; i < executions; ++i)
		{
			// Reset initial state in case this is not the first execution:
			// FIXME test_loop() needs to use updatePowerAndVoltage(u_0)
			int unimplemented;
//			cc.updatePowerAndVoltage(u_0);
			
			// Loop:
			cc.loop(iterations, 0.5);
			
			// Check that the schedule is identical to previous executions:
			RealVector[] newSchedule = cc.controlSchedule();
			if(schedule != null)
			{
//				for (int t = 0; t < newSchedule.length; t++)
//				{
//					for(int j = 0; j < schedule[t].getDimension(); ++j)
//					{
//						assertTrue(schedule[t].getEntry(j) == newSchedule[t].getEntry(j), "test_loop(): Shedule found in execution "+i+" did not match previous schedule for time "+t+" ("+schedule[t].getEntry(j)+" != "+newSchedule[t].getEntry(j)+");\n"
//								+ "previous schedule=,"+printVector(schedule[t])+"\n"
//								+ "new schedule=,"+printVector(newSchedule[t]));
//					}
//				}
				
				for(int j = 0; j < cc.decomposition.length; ++j)
				{
					for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
					{
						assertTrue(equals(cc.ctgs[i][t], ctgs[i][t], 1e-3), "test_loop(): CTG found in execution "+i+" did not match previous CTG for time "+t+" ("+ctgs[i][t]+" != "+cc.ctgs[i][t]+")");
					}
				}
			}
			schedule = newSchedule;
			ctgs = cc.ctgs;
		}
		
		if(verbose)
		{
			System.out.println(callback.externalVoltageData.toString());
			System.out.println(callback.controlData.toString());
			System.out.println(callback.ctgData.toString());
			System.out.println(callback.voltageData.toString());
		}
		
		endTest("test_loop(C_max="+C_max+",U_max="+U_max+",controllableCount="+controllableCount+",executions="+executions+",withExternalVoltageUpdates="+withExternalVoltageUpdates+")");
	}
	
	private class test_loop_LoopFinishedCallback implements Runnable
	{
		public boolean withExternalVoltageUpdates;
		private GridlessCentralCoordinator cc;
		private int iteration = 0;
		public StringBuffer externalVoltageData = new StringBuffer("\ndeltaV_external:\n");
		public StringBuffer controlData = new StringBuffer("\nControls:\n");
		public StringBuffer ctgData = new StringBuffer("CTGs\n");
		public StringBuffer voltageData = new StringBuffer("\n");
		
		List<List<Pair<RealVector, Double>>> previousDPResults; // One for each subset

		public test_loop_LoopFinishedCallback(GridlessCentralCoordinator cc)
		{
			this.cc = cc;
			printX_0(cc);
		}

		public void printX_0(GridlessCentralCoordinator cc)
		{
			if(verbose)
			{
				System.out.println("Iteration "+iteration);
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					System.out.print("DG,S,L:,"+cc.decomposition(i).getGridlessData().dgCount+","+cc.decomposition(i).getGridlessData().storageCount+","+cc.decomposition(i).getGridlessData().loadCount+",");
					System.out.print("x_0:,"+printVector(cc.decomposition(i).getX0()));
					System.out.println("v_0:,"+printVector(cc.decomposition(i).getGridlessData().v_0));
				}
			}
		}

		/**
		 * Performs DP for each decomposition.
		 * @param cc
		 * @return The result of {@link GridTestHelper#dp(DiscreteGridlessADP)} for each decomposition.
		 */
		public List<List<Pair<RealVector, Double>>> dpComparison(GridlessCentralCoordinator cc)
		{
			List<List<Pair<RealVector, Double>>> dpResults = new ArrayList<>();
			for(int i = 0; i < cc.decomposition.length; ++i)
			{
				List<Pair<RealVector, Double>> dpResult = dp((DiscreteGridlessADP)cc.decomposition(i));
				dpResults.add(dpResult);
			}
			return dpResults;
		}

		@Override
		public void run()
		{
			// Test that the DP of each subset with the updated x_0 does not change:
			// (This only applies if external voltage changes are disabled.)
			if(!withExternalVoltageUpdates)
			{
				printX_0(cc);
				
				List<List<Pair<RealVector, Double>>> dpResults = dpComparison(cc);
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					List<Pair<RealVector, Double>> newResult = dpResults.get(i);
					List<Pair<RealVector, Double>> oldResult = null;
					if(iteration > 0)
					{
						oldResult = previousDPResults.get(i);
						assertTrue(newResult.size() == oldResult.size(), "End loop "+iteration+": Size mismatch for DO comparison.");
					}

					for (int t = 0; t < newResult.size(); t++)
					{
						Pair<RealVector, Double> newCTG = newResult.get(t);
						
						if(iteration > 0)
						{
							// Compare old and new DP results:
							Pair<RealVector, Double> oldCTG = oldResult.get(t);
							if(t < cc.decomposition(i).getT()) // No u_T.
								assertTrue(oldCTG.getKey().equals(newCTG.getKey()), "End loop "+iteration+": Control from DP for subset "+i+" at time "+t+" was different to previous iteration; new = "+printVector(newCTG.getKey())+", old = "+printVector(oldCTG.getKey()));
							if(t > 0) // Cost at t=0 is dependant on x_0 which changes.
							{
								assertTrue(TestHelper.equals(oldCTG.getValue(), newCTG.getValue(), 0.5), "End loop "+iteration+": CTG from DP for subset "+i+" at time "+t+" was different to previous iteration; new = "+newCTG.getValue()+", old = "+oldCTG.getValue()+
										"\n\told control =,"+printVector(oldCTG.getKey())+
										"\n\tnew control =,"+printVector(newCTG.getKey()));
							}
						}
						
						// Compare DP with ADP:
						assertTrue(TestHelper.equals(cc.ctgs[i][t], newCTG.getValue(), 0.5), "End loop "+iteration+": ADP CTG did not match DP CTG for subset "+i+" at time "+t+"; ADP CTG = "+cc.ctgs[i][t]+", DP CTG = "+newCTG.getValue());
					}
				}
				previousDPResults = dpResults;
			}
			
			// TODO Test for convergence:
			if(verbose)
			{
				// External voltages:
				externalVoltageData.append("iteration "+iteration+":");
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					externalVoltageData.append(",Subset "+i+":,");
					for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
					{
						externalVoltageData.append("time "+t+":,");
						externalVoltageData.append(",");
						externalVoltageData.append(printVector(cc.decomposition(i).getGridlessData().deltaV_external[t]));
						externalVoltageData.append(",,");
					}
				}
				externalVoltageData.append(",\n");
				
				// Scheduled controls:
				controlData.append("iteration "+iteration+":");
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					controlData.append(",Subset "+i+":,");
					for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
					{
						controlData.append("time "+t+":,");
						controlData.append(",");
						controlData.append(printVector(cc.schedule[i][t]));
						controlData.append(",,");
					}
				}
				controlData.append(",\n");
			}
			
			// TODO test for ctg reduction
			if(verbose)
			{
				ctgData.append("iteration "+iteration+":");
				ctgData.append(",Full network CTG:,");
				double[] ctg = ctg(cc);
				appendCTGs(ctgData, ctg);
				for(int i = 0; i < cc.decomposition.length; ++i)
				{
					ctgData.append(",Subset "+i+":,");
					double[] ctgs = cc.ctgs[i];
					appendCTGs(ctgData, ctgs);
				}
				ctgData.append(",\n");
			}
			
			// TODO test for voltage error improvements
			
			++iteration;
		}

		/**
		 * For debug only.
		 */
		@SuppressWarnings("unused")
		private void compareOldAndNewDP(
				DiscreteGridlessADP adp,
				List<Pair<RealVector, Double>> oldResult,
				RealVector oldx_0, 
				List<Pair<RealVector, Double>> newResult, 
				RealVector newx_0)
		{
			RealVector[] oldx = new RealVector[adp.T+1];
			oldx[0] = oldx_0;
			RealVector[] newx = new RealVector[adp.T+1];
			newx[0] = newx_0;
			
			for(int t = 0; t <= adp.T; ++t)
			{
				RealVector u_told = oldResult.get(t).getKey();
				RealVector u_tnew = newResult.get(t).getKey();
				
				if(t < adp.T)
				{
					RealVector wZero = adp.wZero(t);
					oldx[t+1] = adp.f(t, oldx[t], u_told, wZero);
					newx[t+1] = adp.f(t, newx[t], u_tnew, wZero);
				}
			}
			
			double ctgNew = 0;
			double ctgOld = 0;
			for(int t = adp.T; t >= 0; --t)
			{
				ctgOld += adp.g(t, oldx[t], null);
				ctgNew += adp.g(t, newx[t], null);
				RealVector oldu_t = oldResult.get(t).getKey();
				RealVector newu_t = newResult.get(t).getKey();
				System.out.println(t+","+printVector(oldx[t])+","+printVector(oldu_t)+","+ctgOld);
				System.out.println(t+","+printVector(newx[t])+","+printVector(newu_t)+","+ctgNew);
			}
		}

		public void appendCTGs(StringBuffer out, double[] ctgs)
		{
			for(int t = 0; t < cc.localControllerTemplate.getT(); ++t)
			{
				out.append("time "+t+":,");
				out.append(",");
				out.append(ctgs[t]);
				out.append(",");
			}
		}
	}
	
	private double[] ctg(GridlessCentralCoordinator cc)
	{
		List<Pair<RealVector, Double>> schedule = dp((DiscreteGridlessADP)cc.gridlessADP());
		double[] ctgs = new double[cc.localControllerTemplate.getT()];
		for (int t = 0; t < ctgs.length; t++)
		{
			ctgs[t] = schedule.get(t).getValue();
		}
		return ctgs;
	}
}
