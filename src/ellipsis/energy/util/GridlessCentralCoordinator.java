package ellipsis.energy.util;

import static java.lang.Math.max;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.MatrixHelper;
import ellipsis.util.VectorHelper;

public class GridlessCentralCoordinator extends CentralCoordinator
{
	// Parameters:
//	public GridlessADP adp;
	public int C_max; // Maximum number of controllable units allowed.
	public int U_max; // Maximum number of units.
	public int passTimeRate = 5; // FIXME Every passTimeRate iterations.
	public boolean ignoreExternal = false; // FIXME ignoreExternal
	
	// Output variables:
//	public GridlessADP[] decomposition;
//	public List<List<Integer>> Bs;
	public List<Integer> pilots; // Global index of pilot busses.
	public List<Integer> pilots_i; // Local index of pilot busses.
	
	private RealVector[] previousUndampenedSchedule;
	private RealVector[] previousDeltaS_hat;
	
	public GridlessADP gridlessADP()
	{
		return (GridlessADP) localControllerTemplate;
	}
	
	public CentralTemplateGridlessADP centralTemplate()
	{
		if(localControllerTemplate instanceof CentralTemplateGridlessADP)
			return (CentralTemplateGridlessADP) localControllerTemplate;
		return null;
	}
	
	public GridlessADP decomposition(int i)
	{
		if(decomposition[i] instanceof GridlessADP)
			return (GridlessADP) decomposition[i];
		else
			return null;
	}
	
	public void decompose()
	{
		// Define maximum forecast power shift:
		RealVector DeltaS_t_max = maxForecastPowerChange();
				
		// Create array of products of sensitivities and power shifts:
		RealMatrix LambdaProduct = ebeLambdaDeltaSProduct(DeltaS_t_max);
		
		// Select C_max most significant elements for each controllable unit:
		List<List<Integer>> Bs = selectSubsets(LambdaProduct); // These are DG and storage only.
		
		// Create local controllers:
		int transformerCount = centralTemplate() != null ? centralTemplate().getTransformerCount() : 0;
		decomposition = new LocalController[Bs.size()+transformerCount];
		probabilityOfUpdate = new double[Bs.size()+transformerCount];
		makeGridlessADP(Bs);
		makeTransformers();
	}

	public void setupTransformerLoadVoltages(GridlessTransformerRegulator transformer)
	{
		GridlessADP gridlessADP = gridlessADP();
		GridlessData data = gridlessADP.getGridlessData();
		
		// Get local sensitivities and power_0s (local should be load busses):
		List<Integer> B = transformer.getBusIndeces();
		int qOffset_B = B.size();
		double[][] sensitivitiesArray_B = new double[qOffset_B*2][qOffset_B*2];
		double[] pq0Array_B = new double[qOffset_B*2];
		int power_0Offset = data.power_0.getDimension()/2;
		int i_B = 0;
		for (Integer i : B)
		{
			int j_B = 0;
			for (Integer j : B)
			{
				sensitivitiesArray_B[i_B][j_B] = data.sensitivities.getEntry(i, j);
				sensitivitiesArray_B[i_B][j_B+qOffset_B] = data.sensitivities.getEntry(i, j+power_0Offset);
				sensitivitiesArray_B[i_B+qOffset_B][j_B] = data.sensitivities.getEntry(i+power_0Offset, j);
				sensitivitiesArray_B[i_B+qOffset_B][j_B+qOffset_B] = data.sensitivities.getEntry(i+power_0Offset, j+power_0Offset);
				++j_B;
			}
			pq0Array_B[i_B] = data.power_0.getEntry(i);
			pq0Array_B[i_B+qOffset_B] = data.power_0.getEntry(i+power_0Offset);
			++i_B;
		}
		RealMatrix sensitivities_B = new Array2DRowRealMatrix(sensitivitiesArray_B, false);
		transformer.setSensitivities(sensitivities_B);
		
		RealVector[] loadPower_B = new RealVector[gridlessADP.getT()];
		for(int t = 0; t < gridlessADP.getT(); ++t)
		{
			// Get load power changes:
			RealVector pq_B = 
							VectorHelper.vector(0.0)
					.append(new ArrayRealVector(data.meanP_L[t]))
					.append(VectorHelper.vector(0.0))
					.append(new ArrayRealVector(data.meanQ_L[t]));
			
			// Set voltages:
			loadPower_B[t] = pq_B;
		}
		transformer.setLoadPower(loadPower_B);
	}

	private void makeTransformers()
	{
		if(centralTemplate() != null)
		{
			int transformerCount = centralTemplate().getTransformerCount();
			GridlessADP gridlessADP = gridlessADP();
			GridlessData gridlessData = gridlessADP.getGridlessData();
			int transformerStartIndex = gridlessData.dgCount+gridlessData.storageCount;
			
			List<GridlessTransformerRegulator> transformers = centralTemplate().getTransformers();
			for(int x = 0; x < transformerCount; ++x)
			{
				int index = transformerStartIndex+x;
				decomposition[index] = transformers.get(x).copy();
				
				// Busses are loads:
				ArrayList<Integer> B = new ArrayList<Integer>();
				B.add(index);
				int loadPowerOffset = gridlessADP.loadPowerOffset();
				for(int i = 0; i < gridlessData.loadCount; ++i)
					B.add(loadPowerOffset+i);
				decomposition[index].setBusIndeces(B);
				
//				decomposition[index].setV0(localControllerTemplate.getV0().copy());
				fillInInitialState((LocalController) decomposition[index], B);
				pilots.add(index);
				pilots_i.add(0);
				
				setupTransformerLoadVoltages((GridlessTransformerRegulator) decomposition[index]);
			}
		}
	}

	public RealVector maxForecastPowerChange()
	{
		RealVector S_t_max = maxForecastPower();
		RealVector DeltaS_t_max = S_t_max.subtract(S_0());
		return DeltaS_t_max;
	}

	public RealVector S_0()
	{
//		int unitCount = adp.dgCount+adp.storageCount+adp.loadCount;
//		return adp.x_0.getSubVector(adp.storageCount, 2*unitCount);
		return gridlessADP().getGridlessData().power_0;
	}

	public void makeGridlessADP(List<List<Integer>> Bs)
	{
		int b = 0;
		for (List<Integer> B : Bs)
		{
			GridlessADP adp_B = (GridlessADP) gridlessADP().makeNew();
			
			adp_B.setCostFunction(gridlessADP().getCostFunction());
			
			fillInSensitivitiesAndCounts(adp_B, B);
			fillInForecastAndInitialState(adp_B, B);
			fillInInitialState(adp_B, B);
			fillInNames(adp_B, B);
			
			// Build GridlessADP:
			adp_B.setname(adp_B.getName() + "("+gridlessADP().getGridlessData().voltageNames.get(pilots.get(b))+")"); // We use votlageNames here because this represents all units.
			adp_B.getGridlessData().storageChargeRate = gridlessADP().getGridlessData().storageChargeRate;
			adp_B.getGridlessData().storageDischargeRate = gridlessADP().getGridlessData().storageDischargeRate;
			adp_B.getGridlessData().storageMaxCapacity = gridlessADP().getGridlessData().storageMaxCapacity;
			adp_B.setT(localControllerTemplate.getT());
			adp_B.setNr(gridlessADP().getNr());
			
			decomposition[b] = adp_B;
			probabilityOfUpdate[b] = 1; // default
			
			adp_B.setBusIndeces(B);
			
			++b;
		}
	}

	private void fillInNames(GridlessADP adp_B, List<Integer> sortedB)
	{
		// Note: In order of DG, Storage, Loads.
		int storageStartIndex = gridlessADP().getGridlessData().dgCount;
		int loadStartIndex = gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount;
		
		adp_B.getGridlessData().stateNames = new HashMap<>();
		adp_B.getGridlessData().controlNames = new HashMap<>();
		adp_B.getGridlessData().noiseNames = new HashMap<>();
		adp_B.getGridlessData().voltageNames = new HashMap<>();
		
		int i_new = 0;
		for (Integer i : sortedB)
		{
			boolean isDG = i < storageStartIndex;
			boolean isLoad = i >= loadStartIndex;
			boolean isStorage = i >= storageStartIndex && i < loadStartIndex;
			
			// State names (storage SOC):
			if(isStorage)
				adp_B.getGridlessData().stateNames.put(i_new, gridlessADP().getGridlessData().stateNames.get(i));
			
			// Control names (DG and Storage):
			if(isDG || isStorage)
				adp_B.getGridlessData().controlNames.put(i_new, gridlessADP().getGridlessData().controlNames.get(i));
			
			// Noise names (DG variation, loads):
			if(isDG || isLoad)
				adp_B.getGridlessData().noiseNames.put(i_new, gridlessADP().getGridlessData().noiseNames.get(i));
			
			// Voltage names (all units):
			adp_B.getGridlessData().voltageNames.put(i_new, gridlessADP().getGridlessData().voltageNames.get(i));
			
			++i_new;
		}
	}

	public void fillInForecastAndInitialState(GridlessADP adp_B, List<Integer> sortedB)
	{
		{
			int i_new = 0;
			
			int T = localControllerTemplate.getT();
			adp_B.setX0(new ArrayRealVector(adp_B.getGridlessData().storageCount));
			adp_B.getGridlessData().meanP_DG = new double[T][];
			adp_B.getGridlessData().meanP_L = new double[T][];
			adp_B.getGridlessData().meanQ_L = new double[T][];
			adp_B.getGridlessData().sdP_DG = new double[T][];
			adp_B.getGridlessData().sdP_L = new double[T][];
			adp_B.getGridlessData().sdQ_L = new double[T][];
			
			for(int t = 0; t < T; ++t)
			{
				adp_B.getGridlessData().meanP_DG[t] = new double[adp_B.getGridlessData().dgCount];
				adp_B.getGridlessData().meanP_L[t] = new double[adp_B.getGridlessData().loadCount];
				adp_B.getGridlessData().meanQ_L[t] = new double[adp_B.getGridlessData().loadCount];
				adp_B.getGridlessData().sdP_DG[t] = new double[adp_B.getGridlessData().dgCount];
				adp_B.getGridlessData().sdP_L[t] = new double[adp_B.getGridlessData().loadCount];
				adp_B.getGridlessData().sdQ_L[t] = new double[adp_B.getGridlessData().loadCount];
			}
			
			for (Integer i : sortedB)
			{
				for(int t = 0; t < T; ++t)
				{
					// Forecast:
					if(i < gridlessADP().getGridlessData().dgCount) // DG
					{
						adp_B.getGridlessData().meanP_DG[t][i_new] = gridlessADP().getGridlessData().meanP_DG[t][i];
						adp_B.getGridlessData().sdP_DG[t][i_new] = gridlessADP().getGridlessData().sdP_DG[t][i];
					}
					else if(i >= gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount) // Loads
					{
						adp_B.getGridlessData().meanP_L[t][i_new-adp_B.getGridlessData().dgCount-adp_B.getGridlessData().storageCount] = gridlessADP().getGridlessData().meanP_L[t][i-gridlessADP().getGridlessData().dgCount-gridlessADP().getGridlessData().storageCount];
						adp_B.getGridlessData().sdP_L[t][i_new-adp_B.getGridlessData().dgCount-adp_B.getGridlessData().storageCount] = gridlessADP().getGridlessData().sdP_L[t][i-gridlessADP().getGridlessData().dgCount-gridlessADP().getGridlessData().storageCount];
						adp_B.getGridlessData().meanQ_L[t][i_new-adp_B.getGridlessData().dgCount-adp_B.getGridlessData().storageCount] = gridlessADP().getGridlessData().meanQ_L[t][i-gridlessADP().getGridlessData().dgCount-gridlessADP().getGridlessData().storageCount];
						adp_B.getGridlessData().sdQ_L[t][i_new-adp_B.getGridlessData().dgCount-adp_B.getGridlessData().storageCount] = gridlessADP().getGridlessData().sdQ_L[t][i-gridlessADP().getGridlessData().dgCount-gridlessADP().getGridlessData().storageCount];
					}
				}
	
				if(i >= gridlessADP().getGridlessData().dgCount && i < (gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount)) // Storage?
				{
					adp_B.getX0().setEntry(i_new-adp_B.getGridlessData().dgCount, gridlessADP().getX0().getEntry(i-gridlessADP().getGridlessData().dgCount)); // SOC
				}
				
				++i_new;
			}
		}
	}

	private void fillInInitialState(LocalController adp_B, List<Integer> sortedB)
	{
		int unitCount = gridlessADP().unitCount();
		int unitCount_B = adp_B.unitCount();
		int i_new = 0;
		int T = localControllerTemplate.getT();

		adp_B.setV0(new ArrayRealVector(unitCount_B*2));
		adp_B.setPower0(new ArrayRealVector(unitCount_B*2));
		
		// Initialise deltaVExternal:
		{
			RealVector[] deltaVExt = new RealVector[T+1];
			for(int t = 0; t <= T; ++t)
			{
				deltaVExt[t] = new ArrayRealVector(2*unitCount_B);
			}
			adp_B.setDeltaVExternal(deltaVExt);
		}
		
		for (Integer i : sortedB)
		{
			for(int t = 0; t <= T; ++t)
			{
				// External delta v:
				double delta = gridlessADP().getGridlessData().deltaV_external[t].getEntry(i);
				double abs = gridlessADP().getGridlessData().deltaV_external[t].getEntry(unitCount+i);
				adp_B.setDeltaVExternal(t, i_new, delta);
				adp_B.setDeltaVExternal(t, unitCount_B+i_new, abs);
			}
			
			// Initial state:
			double P_0 = gridlessADP().getGridlessData().power_0.getEntry(i);
			double Q_0 = gridlessADP().getGridlessData().power_0.getEntry(unitCount+i);
			double delta = gridlessADP().getGridlessData().v_0.getEntry(i);
			double abs = gridlessADP().getGridlessData().v_0.getEntry(unitCount+i);
			
//				adp_B.x_0.setEntry(adp_B.storageCount+i_new, P);
//				adp_B.x_0.setEntry(adp_B.storageCount+unitCount_B+i_new, Q);
			adp_B.setPower0(i_new, P_0);
			adp_B.setPower0(unitCount_B+i_new, Q_0);
			adp_B.setV0(i_new, delta);
			adp_B.setV0(unitCount_B+i_new, abs);
			
			++i_new;
		}
	}

	public void fillInSensitivitiesAndCounts(GridlessADP adp_B, List<Integer> B)
	{
		int unitCount_B = B.size();
		int dgCount = gridlessADP().getGridlessData().dgCount;
		int storageCount = gridlessADP().getGridlessData().storageCount;
		int loadCount = gridlessADP().getGridlessData().loadCount;
		int unitCount = gridlessADP().unitCount();
		adp_B.getGridlessData().sensitivities = new Array2DRowRealMatrix(unitCount_B*2, unitCount_B*2);
//		adp_B.getGridlessData().slackAdmittances = new ArrayFieldVector<Complex>(new Complex[unitCount_B]);
		adp_B.getGridlessData().slackSensitivities = new Array2DRowRealMatrix(2, unitCount_B*2); // 2 rows for arg and abs of slack
		int dgCount_B = 0;
		int storageCount_B = 0;
		int loadCount_B = 0;
		int i_new = 0; // row in new sensitivity matrix
		RealMatrix sensitivities = gridlessADP().getGridlessData().sensitivities;
		RealMatrix slackSensitivities = gridlessADP().getGridlessData().slackSensitivities;
		for (Integer i : B)
		{
			// Copy relevant elements of the sensitivity matrix:
			int j_new = 0; // column in new sensitivity matrix
			for (Integer j : B)
			{
				adp_B.getGridlessData().sensitivities.setEntry(i_new,             j_new,             sensitivities.getEntry(i,           j));
				adp_B.getGridlessData().sensitivities.setEntry(i_new,             j_new+unitCount_B, sensitivities.getEntry(i,           j+unitCount));
				adp_B.getGridlessData().sensitivities.setEntry(i_new+unitCount_B, j_new,             sensitivities.getEntry(i+unitCount, j));
				adp_B.getGridlessData().sensitivities.setEntry(i_new+unitCount_B, j_new+unitCount_B, sensitivities.getEntry(i+unitCount, j+unitCount));
				
				++j_new;
			}
			
			// Copy relevant elements of the slack admittance vector:
//			adp_B.getGridlessData().slackAdmittances.setEntry(i_new, adp.slackAdmittances.getEntry(i));
			
			// Copy relevant elements of the slack sensitivity matrix:
			adp_B.getGridlessData().slackSensitivities.setEntry(0, i_new, slackSensitivities.getEntry(0, i));
			adp_B.getGridlessData().slackSensitivities.setEntry(0, unitCount_B+i_new, slackSensitivities.getEntry(0, unitCount+i));
			adp_B.getGridlessData().slackSensitivities.setEntry(1, i_new, slackSensitivities.getEntry(1, i));
			adp_B.getGridlessData().slackSensitivities.setEntry(1, unitCount_B+i_new, slackSensitivities.getEntry(1, unitCount+i));
			
			// Count DG, storage and loads:
			if(i < dgCount)
				++dgCount_B;
			else if(i < dgCount+storageCount)
				++storageCount_B;
			else if(i < dgCount+storageCount+loadCount)
				++loadCount_B;
			
			++i_new;
		}
		
		adp_B.getGridlessData().slackAdmittance = gridlessADP().getGridlessData().slackAdmittance;
		
		adp_B.getGridlessData().dgCount = dgCount_B;
		adp_B.getGridlessData().storageCount = storageCount_B;
		adp_B.getGridlessData().loadCount = loadCount_B;
	}

	public List<List<Integer>> selectSubsets(RealMatrix LambdaProduct)
	{
//for(int i = 0; i < LambdaProduct.getRowDimension(); ++i)
//{
//	for(int j = 0; j < LambdaProduct.getColumnDimension(); ++j)
//	{
//		System.out.print(LambdaProduct.getEntry(i, j));
//		System.out.print(',');
//	}
//	System.out.println();
//}
		int dgCount = gridlessADP().getGridlessData().dgCount;
		int storageCount = gridlessADP().getGridlessData().storageCount;
		int unitCount = gridlessADP().unitCount();
		if(LambdaProduct.getColumnDimension() != 2*unitCount)
			throw new RuntimeException("LambdaProduct dimension was not 2*unitCount: "+LambdaProduct.getColumnDimension()+" != 2*"+unitCount);
		
		List<List<Integer>> Bs = new ArrayList<>();
		pilots = new ArrayList<>();
		pilots_i = new ArrayList<>();
		
		for(int dg = 0; dg < dgCount; ++dg)
		{
			SortedSet<Integer> Bset = createSubset(LambdaProduct, dg);
			List<Integer> B = new ArrayList<>(Bset);
			Bs.add(B);
			pilots.add(dg);
			pilots_i.add(B.indexOf(dg));
		}

		for(int s = 0; s < storageCount; ++s)
		{
			int pilot = dgCount+s;
			SortedSet<Integer> Bset = createSubset(LambdaProduct, pilot);
			List<Integer> B = new ArrayList<>(Bset);
			Bs.add(B);
			pilots.add(pilot);
			pilots_i.add(B.indexOf(pilot));
		}
		
		return Bs;
	}

	private SortedSet<Integer> createSubset(RealMatrix LambdaProduct, int pilot)
	{
		int unitCount = gridlessADP().unitCount();
		SortedSet<Integer> B = new TreeSet<>();
		Set<Integer> ignored = new HashSet<>();
		B.add(pilot);
		int controllableCount = 1; // Already have the controllable bus at index offset.
		double[] argRow = LambdaProduct.getRow(pilot);
		double[] absRow = LambdaProduct.getRow(unitCount+pilot);
		while(B.size() < U_max && B.size()+ignored.size() < unitCount)
		{
			double max = -1;
			int max_i = -1;
			for(int i = 0; i < unitCount; ++i)
			{
				if(!B.contains(i) && !ignored.contains(i))
				{
					double abs = absRow[i]+absRow[unitCount+i]; // P and Q terms
					double arg = argRow[i]+argRow[unitCount+i];
					double value = max(Math.abs(abs), Math.abs(arg));
					if(value > max)
					{
						max = value;
						max_i = i;
					}
				}
			}
			if(max_i == -1)
				throw new RuntimeException("Error searching for most significant elements of LambdaProduct; is C_max ("+C_max+") > unitCount ("+unitCount+")?");
			
			// Add to B if we won't exceed number of controllable units:
			boolean isControllable = max_i < gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount; // Must be DG or storage and therefore controllable.
			if(controllableCount < C_max || !isControllable)
			{
				B.add(max_i);
				if(isControllable)
					++controllableCount;
			}
			else
			{
				ignored.add(max_i);
			}
		}
		return B;
	}

	public RealMatrix ebeLambdaDeltaSProduct(RealVector DeltaS_t_max)
	{
		RealMatrix LambdaProduct = new Array2DRowRealMatrix(gridlessADP().getGridlessData().sensitivities.getData());
		int dimension = DeltaS_t_max.getDimension();
		RealMatrix DeltaS_t_max_diagonal = new Array2DRowRealMatrix(dimension, dimension);
		for(int i = 0; i < dimension; ++i)
		{
			DeltaS_t_max_diagonal.setEntry(i, i, DeltaS_t_max.getEntry(i));
		}
		
		LambdaProduct = LambdaProduct.multiply(DeltaS_t_max_diagonal);
		
		return LambdaProduct;
	}

	public RealVector maxForecastPower()
	{
		int unitCount = gridlessADP().unitCount();
		RealVector DeltaS_t_max = new ArrayRealVector(unitCount*2);
		double maxNorm = 0;
		for(int t = 0; t < localControllerTemplate.getT(); ++t)
		{
			RealVector DeltaS_t = new ArrayRealVector(unitCount*2);
			
			// DG forecast:
			for(int i = 0; i < gridlessADP().getGridlessData().dgCount; ++i)
			{
				DeltaS_t.setEntry(i, gridlessADP().getGridlessData().meanP_DG[t][i]);
				DeltaS_t.setEntry(unitCount+i, 0);
			}

			// Storage
			for(int i = 0; i < gridlessADP().getGridlessData().storageCount; ++i)
			{
				DeltaS_t.setEntry(gridlessADP().getGridlessData().dgCount+i, 0);
				DeltaS_t.setEntry(unitCount+gridlessADP().getGridlessData().dgCount+i, 0);
			}
			
			// Loads:
			int loadOffset = gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount;
			for(int i = 0; i < gridlessADP().getGridlessData().loadCount; ++i)
			{
				DeltaS_t.setEntry(loadOffset + i, gridlessADP().getGridlessData().meanP_L[t][i]);
				DeltaS_t.setEntry(loadOffset + unitCount + i, gridlessADP().getGridlessData().meanQ_L[t][i]);
			}
			
			// Find max norm:
			double norm = DeltaS_t.getNorm();
			if(norm > maxNorm)
			{
				maxNorm = norm;
				DeltaS_t_max = DeltaS_t;
			}
		}

		return DeltaS_t_max;
	}
	
	/**
	 * Calculates the error in voltage calculation in the subset due to unaccounted
	 * for external changes.
	 * @param subset The subset whose pilot bus will be used to calculate the error.
	 * @return The value of epsilon.
	 */
	public double epsilon(int subset)
	{
		if(decomposition(subset) == null)
			return  -1;
		
		int unitCount = gridlessADP().getGridlessData().dgCount+gridlessADP().getGridlessData().storageCount+gridlessADP().getGridlessData().loadCount;
		
		// Get the true voltage change:
		List<Integer> B = decomposition[subset].getBusIndeces();
//		List<Integer> B = Bs.get(subset);
		int unitCount_B = B.size();
		RealVector DeltaS = maxForecastPowerChange();
		RealVector trueVoltageChange = gridlessADP().getGridlessData().sensitivities.operate(DeltaS);
		
		// Get the local changes in power (DeltaS):
		RealVector DeltaS_small = new ArrayRealVector(B.size()*2);
		int i = 0;
		for (Integer i_original : B)
		{
			double P = DeltaS.getEntry(i_original);
			double Q = DeltaS.getEntry(unitCount+i_original);
			DeltaS_small.setEntry(i, P);
			DeltaS_small.setEntry(i+unitCount_B, Q);
			++i;
		}
		
		// Get the local approximate change in voltage:
		RealVector approximateVoltageChange_small = decomposition(subset).getGridlessData().sensitivities.operate(DeltaS_small);
		
		// Expand the local approximate change in voltage to a
		// full sized vector padded with zeros:
		RealVector approximateVoltageChange = new ArrayRealVector(trueVoltageChange.getDimension());
		i = 0;
		for (Integer i_original : B)
		{
			double delta = approximateVoltageChange_small.getEntry(i);
			double abs = approximateVoltageChange_small.getEntry(unitCount_B+i);
			approximateVoltageChange.setEntry(i_original, delta);
			approximateVoltageChange.setEntry(unitCount+i_original, abs);
			++i;
		}
		
		// Compare the true voltage change with the approximate voltage change:
		RealVector error = trueVoltageChange.subtract(approximateVoltageChange);
		double absError = error.getEntry(pilots.get(subset));
		double argError = error.getEntry(unitCount+pilots.get(subset));
		return max(Math.abs(absError), Math.abs(argError));
	}

	public void setADPTemplate(GridlessADP adp)
	{
		if(adp instanceof CentralTemplateGridlessADP)
		{
//			final List<GridlessTransformerRegulator> transformers = ((CentralTemplateGridlessADP)adp).getTransformers();
			addLoopFinishedCallback(new Runnable()
			{
				int j = 0;
				@Override
				public void run()
				{
					++j; // increment here to avoid initial passing of time
//					for (GridlessTransformerRegulator transformer : transformers) // This was changing the wrong OLTCs
//					{
//						transformer.passTime(transformer.getTapChangeDelay()); // One tap change per central iteration.
//					}
					if(j%passTimeRate == 0)
					{
						int transformerCount = centralTemplate().getTransformerCount();
						GridlessData gridlessData = gridlessADP().getGridlessData();
						int transformerStartIndex = gridlessData.dgCount+gridlessData.storageCount;
						for(int x = 0; x < transformerCount; ++x)
						{
							int index = transformerStartIndex+x;
							GridlessTransformerRegulator transformer = (GridlessTransformerRegulator)decomposition[index];
							transformer.passTime(transformer.getTapChangeDelay());
						}
					}
				}
			});
		}
		
		this.localControllerTemplate = adp;
	}
	
	
	//// Optimization ////

	@Override
	protected double[][] makeBandwidthArray(ADPInterface adp, double bw)
	{
		GridlessADP gridlessADP = (GridlessADP)adp;
		return makeBandwidthArray(
				gridlessADP.getGridlessData().storageCount, 
				gridlessADP.getGridlessData().dgCount, 
				gridlessADP.getGridlessData().loadCount, 
				bw, bw, bw, bw);
	}

	/**
	 * 
	 * @return [time][state element]
	 */
	public double[][] makeBandwidthArray(int storageCount, int dgCount, int loadCount, double socBW, double dgBW, double storageBW, double loadBW)
	{
		double[][] bandwidth;
		int T = localControllerTemplate.getT();
		bandwidth = new double[T][storageCount];
		
		// SOC:
		int start = 0;
		int end = storageCount;
		for(int t = 0; t < T; ++t)
			for(int b = start; b < end; ++b)
				bandwidth[t][b] = socBW;

		return bandwidth;
	}

	/**
	 * Gathers power data from local controllers.
	 * Note: schedule lists schedules from the local controller's point of view,
	 *       previousControlSchedule stores the aggregated schedule.
	 */
	@Override
	public RealVector[] controlSchedule()
	{
//		int dgCount = gridlessADP().getGridlessData().dgCount;
		
		// Initialise:
		int T = localControllerTemplate.getT();
		RealVector[] u = new RealVector[T];
		int controlDimension = gridlessADP().controlDimension();
		int qOffset = controlDimension/2;
		if(previousUndampenedSchedule == null)
			previousUndampenedSchedule = new RealVector[T];
		if(previousDeltaS_hat == null)
			previousDeltaS_hat = new RealVector[T];
		
		// Fill in values:
		int qOffset_B = 0;
		int controlDimension_B = 0;
		RealVector[] undampenedSchedule = new RealVector[T];
		for(int t = 0; t < T; ++t)
		{
			// Initialise for time t:
			u[t] = new ArrayRealVector(controlDimension);
			undampenedSchedule[t] = localControllerTemplate.uZero();;
			if(previousUndampenedSchedule[t] == null)
				previousUndampenedSchedule[t] = localControllerTemplate.uZero();
			
			// Get relevant data from each local controller:
			for(int i = 0; i < decomposition.length; ++i)
			{
				// Check for available data:
				if(schedule[i] == null  && previousControlSchedule == null) // May be null if information has been delayed (see optimise()).
					continue; // No previous value, so assume 0 (default from initialisation of u above).
				
				if(decomposition(i) != null)
				{
					controlDimension_B = decomposition(i).controlDimension();
					qOffset_B = controlDimension_B/2;
				}
				else
				{
					qOffset_B = 1; // transformers only have P & Q
					controlDimension_B = 2;//
				}
				
				int pilotIndex = pilots.get(i);
				int pilotIndex_B = pilots_i.get(i);
			
				double pilotControl_P, pilotControl_Q = 0;
				
				// Get P from schedule or default to previous schedule:
				if(schedule[i] != null && schedule[i][t] != null)
					pilotControl_P = schedule[i][t].getEntry(pilotIndex_B);
				else if(previousControlSchedule != null) // No new information this iteration so use previous control value.
					pilotControl_P = previousControlSchedule[t].getEntry(pilotIndex);
				else
					continue;
				undampenedSchedule[t].setEntry(pilotIndex, pilotControl_P);
				
				// Get Q from schedule or default to previous schedule:
				if(controlDimension_B == 2*qOffset_B) // include Q control
				{
					if(schedule[i] != null && schedule[i][t] != null)
						pilotControl_Q = schedule[i][t].getEntry(pilotIndex_B+qOffset_B);
					else // No new information this iteration so use previous control value.
						pilotControl_Q = previousControlSchedule[t].getEntry(pilotIndex+qOffset);
					undampenedSchedule[t].setEntry(pilotIndex+qOffset, pilotControl_Q);
				}
			} // End for t
		}// End for decomposition
		
		// Combine with previous schedule:	
		for(int t = 0; t < T; ++t)
		{
			// Calculate step size:
			double alpha = alpha(previousDeltaS_hat[t], undampenedSchedule[t], previousUndampenedSchedule[t]);
			
			// Combine for each decomposition:
			for(int i = 0; i < decomposition.length; ++i)
			{
				int pilotIndex = pilots.get(i);
				double pilotControl_P = undampenedSchedule[t].getEntry(pilotIndex);
				double pilotControl_Q = undampenedSchedule[t].getEntry(pilotIndex+qOffset);
				
				if(previousControlSchedule != null)
				{
					// P:
					double u_old_P = previousControlSchedule[t].getEntry(pilotIndex);
					double u_new_P = (1-alpha)*u_old_P + alpha*pilotControl_P;
					u[t].setEntry(pilotIndex, u_new_P);

					// Q:
					if(controlDimension_B == 2*qOffset_B)
					{
						double u_old_Q = previousControlSchedule[t].getEntry(pilotIndex+qOffset);
						double u_new_Q = (1-alpha)*u_old_Q + alpha*pilotControl_Q;
						u[t].setEntry(pilotIndex+qOffset, u_new_Q);
					}
				}
				else
				{
					u[t].setEntry(pilotIndex, pilotControl_P);

					if(controlDimension_B == 2*qOffset_B)
					{
						u[t].setEntry(pilotIndex+qOffset, pilotControl_Q);
					}
				}
			}
			
			// Prepare previousDeltaS_hat for next iteration:
			if(previousControlSchedule != null)
				previousDeltaS_hat[t] = u[t].subtract(previousControlSchedule[t]);
		}
		
		// Remember undampened control schedule:
		previousUndampenedSchedule = undampenedSchedule;
		
		return u;
	}

	/**
	 * FIXME This calculation is incorrect: The dampening is too severe and forces convergence
	 *       without allowing for the controllers the required freedom; i.e. sum{alpha} < infty.
	 * @param s 
	 * @param s_hat
	 * @return
	 */
	private double alpha(
			RealVector previousDeltaS_hat, // \Delta \hat{s}^{(j-1)} = \hat{s}^{(j-1)} - \hat{s}^{(j-2)}
			RealVector s,
			RealVector previousS
			)
	{
		// If an alpha has been explicitly set then use it:
		if(alpha > 0)
			return alpha;

		RealMatrix sensitivities = this.gridlessADP().getGridlessData().sensitivities;
		
		if(previousDeltaS_hat == null)
			previousDeltaS_hat = new ArrayRealVector(sensitivities.getColumnDimension()); // default to zero vector

		// \Delta v^{(j)} = \Lambda\Delta \hat{s}^{(j-1)}
		// Change in voltage due to smoothed change in power at previous iteration.
		RealVector deltaV = sensitivities.operate(previousDeltaS_hat);

		// \Delta s^{(j)} = s^{(j)} - s^{(j-1)}
		RealVector deltaS = s.subtract(previousS);

		if(deltaS == null)
			return 1;
		
		RealMatrix N = MatrixHelper.invert(sensitivities);
		double gamma = 0.8;

		// \hat{\alpha} = 
		//                min{ 
		//                     ((\gamma\Lambda^{-1}|\Deltav^{(j)}|)_i - |\Delta\hat{s}^{(j-1)}_i)
		//                     /
		//                     (|\Delta s^{(j)}_i| - |\Delta\hat{s}^{(j-1_}_i|)
		//                }
		RealVector alpha_hat_num = N.operate(abs(deltaV)).mapMultiply(gamma).subtract(abs(previousDeltaS_hat));
		RealVector alpha_hat_den = abs(deltaS).subtract(abs(previousDeltaS_hat));

		int dimension = alpha_hat_num.getDimension();
		double alpha_min = Double.MAX_VALUE;
		for (int i = 0; i < dimension; i++)
		{
			double alpha_hat_i = alpha_hat_den.getEntry(i)/alpha_hat_num.getEntry(i);
			if(alpha_hat_i < alpha_min)
			{
				alpha_min = alpha_hat_i;
			}
		}

		return alpha_min < 1 ? alpha_min : 1.0;
	}

	private RealVector abs(RealVector v)
	{
		int dimension = v.getDimension();
		ArrayRealVector abs = new ArrayRealVector(dimension);
		for (int i = 0; i < dimension; i++)
		{
			abs.setEntry(i, Math.abs(v.getEntry(i)));
		}
		return abs;
	}

	/**
	 * Sets {@link #adpTemplate}.power_0 according to the given control at t=0, u_0. Updates {@link #adpTemplate}.v_0 and
	 * all {@link #localADP}s' v_0s.
	 * @param u_0
	 */
	@Override
	public void update(RealVector[] u)
	{
		GridlessADP gridlessADP = gridlessADP();
		GridlessData gridlessData = gridlessADP.getGridlessData();
		
		// Reset global external voltage changes (v_0 should not include this):
		for(int t = 0; t < localControllerTemplate.getT(); ++t)
		{
			int dimension = gridlessData.deltaV_external[t].getDimension();
			gridlessData.deltaV_external[t] = new ArrayRealVector(dimension);
		}
		
		// Update global values:
		RealVector w = gridlessADP.wZero(0);
		RealVector newV_0 = gridlessADP.voltagesFromControlAndNoise(0, u[0], w);
		gridlessData.v_0 = newV_0;
		gridlessData.power_0 = gridlessADP.powerFromControlAndNoise(u[0], w);
		
		// Update each subset:
		int QOffset_full = gridlessADP.unitCount();
		int absOffset_full = gridlessADP.unitCount();
		for(int i = 0; i < decomposition.length; ++i)
		{
//			GridlessData gridlessData_B = decomposition(i).getGridlessData();
			
			// Update decomposition[i].power_0 and .v_0:
			int QOffset_B = decomposition[i].unitCount();
			int absOffset_B = decomposition[i].unitCount();
			List<Integer> B = decomposition[i].getBusIndeces();
			int k = 0;
			for (Integer b : B)
			{
				// Ps & Qs (SOC doesn't change):
				double P = gridlessData.power_0.getEntry(b);
				double Q = gridlessData.power_0.getEntry(b+QOffset_full);
				decomposition[i].setPower0(k, P);
				decomposition[i].setPower0(k+QOffset_B, Q);
				
				// Voltages:
				double arg = gridlessData.v_0.getEntry(b);
				double abs = gridlessData.v_0.getEntry(b+absOffset_full);
				decomposition[i].setV0(k, arg);
				decomposition[i].setV0(k+absOffset_B, abs);
				
				++k;
			}
		}

		// Update external voltages:
		if(!ignoreExternal)
		{
			RealVector[][] deltaV_external = externalVoltageChanges(u);
			for(int i = 0; i < decomposition.length; ++i)
			{
				decomposition[i].setDeltaVExternal(deltaV_external[i]);
			}
		}
	}
	
	/**
	 * @return External voltage changes for each subset and each time respectively.
	 */
	public RealVector[][] externalVoltageChanges(RealVector[] u)
	{
		GridlessData gridlessData = gridlessADP().getGridlessData();
		
		// Initialise return array:
		RealVector[][] externalVoltageChanges = new RealVector[decomposition.length][];
		int T = localControllerTemplate.getT();
		for(int i = 0; i < decomposition.length; ++i)
		{
			externalVoltageChanges[i] = new RealVector[T+1];
		}

		// Get external changes for each time and each local controller:
		RealVector PQ_0 = gridlessData.power_0;
		int globalSize = PQ_0.getDimension()/2;
//		double[] xfmrSchedule = transformerRegulator.getRatioSchedule(); // ratios: 1 => no change
		for(int t = 0; t <= T; ++t)
		{
			// Get [\Delta P \Delta Q]:
			RealVector wZero = gridlessADP().wZero(t);
			RealVector u_t = t == T ? localControllerTemplate.uZero() : u[t];
			
			RealVector PQ_t = gridlessADP().powerFromControlAndNoise(u_t, wZero);
			RealVector DeltaPQ_t = PQ_t.subtract(PQ_0);
			
			// Calculate Delta v_t = \Lambda [\Delta P \Delta Q]:
			RealMatrix Lambda_t = gridlessData.sensitivities;
			RealVector DeltaV_t = Lambda_t.operate(DeltaPQ_t);
			
			// For each local controller get voltage changes:
			for (int i = 0; i < decomposition.length; ++i)
			{
				RealVector DeltaPQ_B_t;
				RealMatrix Lambda_BB_t;
				RealVector DeltaV_B_t;
				
				List<Integer> B = decomposition[i].getBusIndeces();
				int Bsize = B.size();
				
				// Get [\Delta P \Delta Q]_B, \Lambda_{B,B} and \Lambda_B [\Delta P \Delta Q]:
				DeltaPQ_B_t = new ArrayRealVector(2*Bsize);
				Lambda_BB_t = new Array2DRowRealMatrix(2*Bsize, 2*Bsize);
				DeltaV_B_t = new ArrayRealVector(2*Bsize);
				
				{
					int j = 0;
					for (Integer b : B)
					{
						// [\Delta P \Delta Q]_B:
						double DeltaP_b_t = DeltaPQ_t.getEntry(b);
						double DeltaQ_b_t = DeltaPQ_t.getEntry(b+globalSize);
						DeltaPQ_B_t.setEntry(j, DeltaP_b_t);
						DeltaPQ_B_t.setEntry(j+Bsize, DeltaQ_b_t);
						
						// \Lambda_{B,B}:
						int k = 0;
						for (Integer c : B)
						{
							double Lambda_bc_t_11 = Lambda_t.getEntry(b, c);
							double Lambda_bc_t_12 = Lambda_t.getEntry(b, c+globalSize);
							double Lambda_bc_t_21 = Lambda_t.getEntry(b+globalSize, c);
							double Lambda_bc_t_22 = Lambda_t.getEntry(b+globalSize, c+globalSize);
							Lambda_BB_t.setEntry(j, k, Lambda_bc_t_11);
							Lambda_BB_t.setEntry(j, k+Bsize, Lambda_bc_t_12);
							Lambda_BB_t.setEntry(j+Bsize, k, Lambda_bc_t_21);
							Lambda_BB_t.setEntry(j+Bsize, k+Bsize, Lambda_bc_t_22);
							
							++k;
						}
	
						// \Delta v_B = \Lambda_B [\Delta P \Delta Q]:
						double DeltaV_t_arg = DeltaV_t.getEntry(b);
						double DeltaV_t_abs = DeltaV_t.getEntry(b+globalSize);
						DeltaV_B_t.setEntry(j, DeltaV_t_arg);
						DeltaV_B_t.setEntry(j+Bsize, DeltaV_t_abs);
						
						++j;
					}
				}
				
				// Calculate \DeltaV_{B,B} = \Lambda_{B,B} [\Delta P \Delta Q]_B:
				RealVector DeltaV_BB_t = Lambda_BB_t.operate(DeltaPQ_B_t);
				
				// Calculate \Delta v_BnotB = \Lambda_B [\Delta P \Delta Q] - \Lambda_{B,B} [\Delta P \Delta Q]_B:
				RealVector DeltaV_BnotB_t = DeltaV_B_t.subtract(DeltaV_BB_t);
				
				externalVoltageChanges[i][t] = DeltaV_BnotB_t;
				
				// Update magnitude portion according to transformer tap position:
//				if(t < T)
//				{
//					int dimension = externalVoltageChanges[i][t].getDimension();
//					for(int j = dimension/2; j < dimension; ++j)
//					{
//						double extDeltaV_i_t_j = externalVoltageChanges[i][t].getEntry(j);
//						externalVoltageChanges[i][t].setEntry(j, extDeltaV_i_t_j + xfmrSchedule[t]-1);
//					}
//				}
			} // decompositions i
		} // time t

		return externalVoltageChanges;
	}
	
//	public RealVector[][] externalVoltageChanges_old(RealVector[] u)
//	{
//		// Initialise return array:
//		RealVector[][] externalVoltageChanges = new RealVector[Bs.size()][];
//		for(int i = 0; i < Bs.size(); ++i)
//		{
//			externalVoltageChanges[i] = new RealVector[localControllerTemplate.getT()+1];
//		}
//		
//		// Fill in for each time and each subset:
////		RealVector x_t = adp.x_0;
////		RealVector PQ_0 = adp.powerFromState(adp.x_0);
//		RealVector PQ_0 = gridlessADP().getGridlessData().power_0;
//		for(int t = 0; t <= localControllerTemplate.getT(); ++t)
//		{
//			RealVector wZero = gridlessADP().wZero(t);
//			RealVector u_t = t == localControllerTemplate.getT() ? localControllerTemplate.uZero() : u[t];
//			
//			for(int i = 0; i < Bs.size(); ++i)
//			{
//				// Lambda_~B (zero out B's rows and columns):
//				RealMatrix Lambda_notB_t = new Array2DRowRealMatrix(gridlessADP().getGridlessData().sensitivities.getData());
//				List<Integer> B = Bs.get(i);
//				int colDimension = Lambda_notB_t.getColumnDimension();
//				int rowDimension = Lambda_notB_t.getRowDimension();
//				for (Integer b : B)
//				{
//					// Clear column of B's elements (Delta Ps & Qs will be ignored except for external nodes):
//					for(int j = 0; j < rowDimension; ++j)
//					{
//						Lambda_notB_t.setEntry(j, b, 0); // P column
//						Lambda_notB_t.setEntry(j, b+colDimension/2, 0); // Q column	
//					}
//				}
//				
//				// DeltaPQ_~B_t:
//				RealVector PQ_t = gridlessADP().powerFromControlAndNoise(u_t, wZero);
//				RealVector DeltaPQ_t = PQ_t.subtract(PQ_0);
//				RealVector DeltaPQ_notB_t = new ArrayRealVector(DeltaPQ_t);
//				int QOffset = DeltaPQ_notB_t.getDimension()/2;
//				for (Integer b : B)
//				{
//					DeltaPQ_notB_t.setEntry(b, 0);
//					DeltaPQ_notB_t.setEntry(b+QOffset, 0);
//				}
//				
//				// Lambda_~B_t * DeltaPQ_~B_t (shrinking to local subset dimension):
//				RealVector DeltaV = Lambda_notB_t.operate(DeltaPQ_notB_t);
//				int voltageDimension_B = decomposition(i).voltageDimension();
//				externalVoltageChanges[i][t] = new ArrayRealVector(voltageDimension_B);
//				int j = 0;
//				int absOffset_B = voltageDimension_B/2;
//				int absOffset_full = DeltaV.getDimension()/2;
//				for (Integer b : B)
//				{
//					double arg = DeltaV.getEntry(b);
//					double abs = DeltaV.getEntry(b+absOffset_full);
//					externalVoltageChanges[i][t].setEntry(j, arg);
//					externalVoltageChanges[i][t].setEntry(j+absOffset_B, abs);
//					++j;
//				}
//			}
//			
////			x_t = adp.f(x_t, u_t, wZero);
//		}
//		
//		return externalVoltageChanges;
//	}
	
	
	//// Debug ////
	
	public Object showDecompositions = new Object()
	{
		public String toString() 
		{
			StringBuffer sb = new StringBuffer();
			for(int i = 0; i < decomposition.length; ++i)
			{
				Integer pilotIndex = pilots.get(i);
				String pilotName = gridlessADP().getGridlessData().controlNames.get(pilotIndex);
				
				sb.append("Subset for pilot ");
				sb.append(pilotName);
				sb.append("(");
				sb.append(pilotIndex);
				sb.append(")");
				sb.append(":\n");
				
				for (Integer b : decomposition[i].getBusIndeces())
				{
					String name = gridlessADP().getGridlessData().voltageNames.get(b);
					sb.append("\t");
					sb.append(name);
					sb.append("(");
					sb.append(b);
					sb.append(")");
					sb.append("\n");
				}
			}
			return sb.toString();
		};
	};
}
