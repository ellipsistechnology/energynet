package ellipsis.energy.util;

import static ellipsis.util.VectorHelper.vector;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Child;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Transformer;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.util.NonParametricApproximator;

/**
 * 
 * @author bmillar
 */
public class GridlessADPFactory
{
	public static GridlessADP fromGrid(Grid grid, double baseVoltage, double basePower, boolean failIfNotConvergent)
	{
		DiscreteGridlessADP adp = new DiscreteGridlessADP();
		return fromGridToGridless(grid, adp, baseVoltage, basePower, failIfNotConvergent);
	}

	public static GridlessADP fromGridToGridless(Grid grid, GridlessADP adp, double baseVoltage, double basePower, boolean failIfNotConvergent)
	{
		double baseImpedance = baseVoltage*baseVoltage/basePower;
		
		// Setup counts from grid:
		Collection<DistributedSource> dgs = grid.getDistributedSources(false);
		Collection<ControllableDemand> storages = grid.get(ControllableDemand.class);
		Collection<Load> loads = grid.getLoads(false);
		Collection<Transformer> transformers = grid.get(Transformer.class);
		
		adp.getGridlessData().loadCount = loads.size();
		adp.getGridlessData().dgCount = dgs.size();
		adp.getGridlessData().storageCount = storages.size();
		
		// Setup transformer:
		if(adp instanceof CentralTemplateGridlessADP)
		{
			((CentralTemplateGridlessADP)adp).setTransformerCount(transformers.size());
			for (Transformer transformer : transformers)
			{
				GridlessTransformerRegulator gtr = new GridlessTransformerRegulator();
				
				gtr.setMinStep((int)Math.round((transformer.getMinimumT()-1)/transformer.getTStepSize()));
				gtr.setMaxStep((int)Math.round((transformer.getMaximumT()-1)/transformer.getTStepSize()));
				gtr.setStepSize(transformer.getTStepSize());

//				gtr.setBusCount(1); // Only supporting a single bus at the moment.
				int T = adp.getT();
				gtr.setT(T);
				
				double tap[] = new double[T];
				for(int t = 0; t < T; ++t)
					tap[t] = 1.0;
				gtr.setRatioSchedule(tap);
				
				gtr.setTapChangeDelay(5); // seconds
				gtr.setVoltageReference(1.0);
				gtr.setAdmittance(transformer.admittance().multiply(baseImpedance));
				
				((CentralTemplateGridlessADP)adp).addTransformer(gtr);
			}
		}
		
		// Setup name mappings.
		// Control and noise are setup here, state and voltage are set up below:
		adp.getGridlessData().stateNames = new HashMap<>();
		adp.getGridlessData().controlNames = new HashMap<>();
		adp.getGridlessData().noiseNames = new HashMap<>();
		adp.getGridlessData().voltageNames = new HashMap<>();
		int i = 0;
		for (DistributedSource dg : dgs)
		{
			adp.getGridlessData().controlNames.put(i, dg.getName());
			adp.getGridlessData().noiseNames.put(i, dg.getName());
			++i;
		}
		for (ControllableDemand storage : storages)
		{
			adp.getGridlessData().controlNames.put(i, storage.getName());
			++i;
		}
		if(adp instanceof CentralTemplateGridlessADP)
		{
			for (Transformer transformer : transformers)
			{
				adp.getGridlessData().controlNames.put(i, transformer.getName());
				++i;
			}
		}
		i = adp.getGridlessData().noiseNames.size();
		for (Load load : loads)
		{
			adp.getGridlessData().noiseNames.put(i, load.getName());
			++i;
		}
		
		// Setup storage parameters:
		adp.getGridlessData().storageChargeRate = new double[adp.getGridlessData().storageCount];
		adp.getGridlessData().storageDischargeRate = new double[adp.getGridlessData().storageCount];
		adp.getGridlessData().storageMaxCapacity = new double[adp.getGridlessData().storageCount];
		i = 0;
		for (ControllableDemand storage : storages)
		{
			adp.getGridlessData().stateNames.put(i, storage.getName());
			
			double chargeRate = -storage.getMaxChargeRate();
			double dischargeRate = storage.getMaxDischargeRate();
			double maxCapacity = storage.getMaxCapacity();
			adp.getGridlessData().storageChargeRate[i] = chargeRate/basePower;
			adp.getGridlessData().storageDischargeRate[i] = dischargeRate/basePower;
			adp.getGridlessData().storageMaxCapacity[i] = maxCapacity/basePower;
			++i;
		}
		
		// Setup x_0 from grid:
		adp.setX0(new ArrayRealVector(adp.getGridlessData().storageCount));
		i = 0;
		for (ControllableDemand storage : storages)
		{
			adp.getX0().setEntry(i, storage.getCapacity()/basePower); // SOC
			++i;
		}

		// Setup power [P_DG P_S P_L Q_DG Q_S Q_L] @ t=0 :
		int unitCount = adp.unitCount();
		adp.getGridlessData().power_0 = new ArrayRealVector(2*unitCount); // P Q
		i = 0;
		for (DistributedSource dg : dgs)
		{
			Complex out = dg.getPowerOutput();
			adp.getGridlessData().power_0.setEntry(i, out.getReal()/basePower);
			adp.getGridlessData().power_0.setEntry(i+unitCount, out.getImaginary()/basePower);
			++i;
		}
		for (ControllableDemand storage : storages)
		{
			Complex in = storage.getChargeRate();
			adp.getGridlessData().power_0.setEntry(i, -in.getReal()/basePower);
			adp.getGridlessData().power_0.setEntry(i+unitCount, -in.getImaginary()/basePower);
			++i;
		}
		i += transformers.size(); // Skip transformers.
		for (Load load : loads)
		{
			Complex in = load.getLoad();
			adp.getGridlessData().power_0.setEntry(i, -in.getReal()/basePower);
			adp.getGridlessData().power_0.setEntry(i+unitCount, -in.getImaginary()/basePower);
			++i;
		}
//		for (Transformer transformer : transformers)
//		{
//			adp.getGridlessData().power_0.setEntry(i, 0);
//			adp.getGridlessData().power_0.setEntry(i+unitCount, 0);
//			++i;
//		}

		// Initialise deltaV_external:
		adp.getGridlessData().setDeltaVExternal(vector(unitCount*2, 0.0));

		// Analyse grid:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(basePower);
        lfa.setBaseVoltage(baseVoltage);
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
        {
        	if(failIfNotConvergent)
        		throw new RuntimeException("Alanysis did not converge!");
        	else
        		System.out.println("WARNING: Alanysis did not converge!\n");
        }
		
		// Setup v_0 from analysis results:
        adp.setV0(new ArrayRealVector(2*unitCount));
        i = 0;
		i = setV_0Voltages(adp, unitCount, i, results, dgs);
		i = setV_0Voltages(adp, unitCount, i, results, storages);
		i = setV_0VoltagesTransformer(adp, unitCount, i, results, transformers);
		i = setV_0Voltages(adp, unitCount, i, results, loads);

//        if(adp instanceof CentralTemplateGridlessADP)
//        {
//	        for (GridlessTransformerRegulator gtr : ((CentralTemplateGridlessADP)adp).getTransformers())
//			{
//				gtr.setV0(adp.getV0());
//	        }
//	    }
		
		// Admittances and sensitivities to slack bus:
		int slackIndex = results.getSlackIndex();
		Map<String, Integer> busNumbers = results.getBusNumbers();
		adp.getGridlessData().slackAdmittance = results.getAdmittanceMatrix().Y.getEntry(slackIndex, slackIndex);

		RealMatrix jacobian = results.getJacobian();
		DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
        RealMatrix sensitivities = solver.getInverse();
		adp.getGridlessData().slackSensitivities = new Array2DRowRealMatrix(2, unitCount*2); // 2 rows: arg and abs of voltages with respect to power
		{
			double[] argRow = new double[2*unitCount];
			double[] absRow = new double[2*unitCount];
			
			int i_new = 0;
			i_new = fillSlackSensitivitesRow(dgs, sensitivities, unitCount, argRow, absRow, slackIndex, busNumbers, i_new);
			i_new = fillSlackSensitivitesRow(storages, sensitivities, unitCount, argRow, absRow, slackIndex, busNumbers, i_new);
			i_new = fillSlackSensitivitesRow(loads, sensitivities, unitCount, argRow, absRow, slackIndex, busNumbers, i_new);

			adp.getGridlessData().slackSensitivities.setRow(0, argRow);
			adp.getGridlessData().slackSensitivities.setRow(1, absRow);
		}
        
        // Set sensitivities from analysis results:
        int oldOffset = sensitivities.getRowDimension()/2;
		adp.getGridlessData().sensitivities = new Array2DRowRealMatrix(unitCount*2, unitCount*2);
		{
			int i_new = 0;
			i_new = convertSensitivities_rows(adp, dgs, storages, loads, transformers, sensitivities, oldOffset, unitCount, busNumbers, slackIndex, i_new, dgs);
			i_new = convertSensitivities_rows(adp, dgs, storages, loads, transformers, sensitivities, oldOffset, unitCount, busNumbers, slackIndex, i_new, storages);
			i_new = convertSensitivitiesTransformer_rows(adp, dgs, storages, loads, transformers, sensitivities, oldOffset, unitCount, busNumbers, slackIndex, i_new, transformers);
			i_new = convertSensitivities_rows(adp, dgs, storages, loads, transformers, sensitivities, oldOffset, unitCount, busNumbers, slackIndex, i_new, loads);
		}
		
		return adp;
	}

	private static int fillSlackSensitivitesRow(Collection<? extends Child> children, RealMatrix sensitivities, int qOffset, double[] argRow, double[] absRow, int slackIndex, Map<String, Integer> busNumbers, int i_new)
	{
		int qOffset_old = sensitivities.getRowDimension()/2;
		
		for (Child child : children)
		{
			String busName_i = child.getParent().getName();
			int i_old = busNumbers.get(busName_i);
			if(i_old > slackIndex)
				--i_old;
			
			argRow[i_new] = sensitivities.getEntry(slackIndex, i_old); // d arg/d P
			argRow[i_new+qOffset] = sensitivities.getEntry(slackIndex, qOffset_old+i_old); // d arg/d Q
			absRow[i_new] = sensitivities.getEntry(qOffset_old+slackIndex, i_old); // d abs/d P
			absRow[i_new+qOffset] = sensitivities.getEntry(qOffset_old+slackIndex, qOffset_old+i_old); // d abs/d Q
			
			++i_new;
		}
		
		return i_new;
	}

//	public static int fillSlackAdmittances(GridlessADP adp,
//			AnalysisResults results, int slackIndex,
//			Map<String, Integer> busNumbers,
//			Collection<? extends Child> children,
//			int i_new)
//	{
//		for (Child child : children)
//		{
//			String busName_i = child.getParent().getName();
//			int i_old = busNumbers.get(busName_i);
//
//			Complex entry = results.getAdmittanceMatrix().Y.getEntry(slackIndex, i_old);
//			adp.slackAdmittances.setEntry(i_new, entry);
//			
//			++i_new;
//		}
//		return i_new;
//	}

	private static int convertSensitivities_rows(
			GridlessADP adp,
			Collection<DistributedSource> dgs, 
			Collection<ControllableDemand> storages,
			Collection<Load> loads, 
			Collection<Transformer> transformers,
			RealMatrix sensitivities, 
			int oldOffset,
			final int qOffset, 
			Map<String, Integer> busNumbers,
			int slackIndex,
			int i_new, Collection<? extends Child> children)
	{
		for (Child child : children)
		{
			String busName_i = child.getParent().getName();
			int i_old = busNumbers.get(busName_i);
			if(i_old > slackIndex)
				--i_old;
			int j_new = 0;
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, dgs);
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, storages);
			j_new = convertSensitivitiesTransformer_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, transformers);
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, loads);
			
			++i_new;
		}
		return i_new;
	}

	private static int convertSensitivitiesTransformer_rows(
			GridlessADP adp,
			Collection<DistributedSource> dgs, 
			Collection<ControllableDemand> storages,
			Collection<Load> loads, 
			Collection<Transformer> transformers,
			RealMatrix sensitivities, 
			int oldOffset,
			final int qOffset, 
			Map<String, Integer> busNumbers,
			int slackIndex,
			int i_new, Collection<? extends Transformer> children)
	{
		for (Transformer child : children)
		{
			String busName_i = child.getToBus().getName();
			int i_old = busNumbers.get(busName_i);
			if(i_old > slackIndex)
				--i_old;
			int j_new = 0;
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, dgs);
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, storages);
			j_new = convertSensitivitiesTransformer_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, transformers);
			j_new = convertSensitivities_columns(adp, qOffset, sensitivities, oldOffset, busNumbers, slackIndex, i_new, i_old, j_new, loads);
			
			++i_new;
		}
		return i_new;
	}

	private static int convertSensitivities_columns(
			GridlessADP adp, 
			final int qOffset,
			RealMatrix sensitivities, 
			int oldOffset,
			Map<String, Integer> busNumbers,
			int slackIndex,
			int i_new, 
			int i_old, 
			int j_new,
			Collection<? extends Child> children)
	{
		for (Child child : children)
		{
			String busName_j = child.getParent().getName();
			int j_old = busNumbers.get(busName_j);
			if(j_old > slackIndex)
				--j_old;
			
			double dDel_dP = sensitivities.getEntry(i_old, j_old);
			double dDel_dQ = sensitivities.getEntry(i_old, j_old+oldOffset);
			double dV_dP = sensitivities.getEntry(i_old+oldOffset, j_old);
			double dV_dQ = sensitivities.getEntry(i_old+oldOffset, j_old+oldOffset);
			
			adp.getGridlessData().sensitivities.setEntry(i_new, j_new, dDel_dP);
			adp.getGridlessData().sensitivities.setEntry(i_new, j_new+qOffset, dDel_dQ);
			adp.getGridlessData().sensitivities.setEntry(i_new+qOffset, j_new, dV_dP);
			adp.getGridlessData().sensitivities.setEntry(i_new+qOffset, j_new+qOffset, dV_dQ);
			
			++j_new;
		}
		return j_new;
	}

	private static int convertSensitivitiesTransformer_columns(
			GridlessADP adp, 
			final int qOffset,
			RealMatrix sensitivities, 
			int oldOffset,
			Map<String, Integer> busNumbers,
			int slackIndex,
			int i_new, 
			int i_old, 
			int j_new,
			Collection<? extends Transformer> children)
	{
		for (Transformer child : children)
		{
			String busName_j = child.getToBus().getName();
			int j_old = busNumbers.get(busName_j);
			if(j_old > slackIndex)
				--j_old;
			
			double dDel_dP = sensitivities.getEntry(i_old, j_old);
			double dDel_dQ = sensitivities.getEntry(i_old, j_old+oldOffset);
			double dV_dP = sensitivities.getEntry(i_old+oldOffset, j_old);
			double dV_dQ = sensitivities.getEntry(i_old+oldOffset, j_old+oldOffset);
			
			adp.getGridlessData().sensitivities.setEntry(i_new, j_new, dDel_dP);
			adp.getGridlessData().sensitivities.setEntry(i_new, j_new+qOffset, dDel_dQ);
			adp.getGridlessData().sensitivities.setEntry(i_new+qOffset, j_new, dV_dP);
			adp.getGridlessData().sensitivities.setEntry(i_new+qOffset, j_new+qOffset, dV_dQ);
			
			++j_new;
		}
		return j_new;
	}

	private static int setV_0Voltages(GridlessADP adp, final int qOffset, int i, AnalysisResults results, Collection<? extends Child> children)
	{
		for (Child child : children)
		{
			adp.getGridlessData().voltageNames.put(i, ((Unit)child).getName());
			
			String busName = child.getParent().getName();
			Complex v = results.getBusVoltage(busName);
			adp.getGridlessData().v_0.setEntry(i, v.getArgument());
			adp.getGridlessData().v_0.setEntry(i+qOffset, v.abs());
			++i;
		}
		return i;
	}

	private static int setV_0VoltagesTransformer(GridlessADP adp, final int qOffset, int i, AnalysisResults results, Collection<? extends Transformer> xfmrs)
	{
		for (Transformer xfmr : xfmrs)
		{
			adp.getGridlessData().voltageNames.put(i, xfmr.getName());
			
			String busName = xfmr.getToBus().getName();
			Complex v = results.getBusVoltage(busName);
			adp.setV0(i, v.getArgument());
			adp.setV0(i+qOffset, v.abs());
			
			++i;
		}
		return i;
	}

	/**
	 * Warning: Does not clone {@link DiscreteGridlessADP#approxJ} or other parameters
	 * that are considered constants. However a call to {@link DiscreteGridlessADP#train(double[])} will overwrite this anyway.
	 * @param from
	 * @param to
	 * @return
	 */
	public static DiscreteGridlessADP copyInto(DiscreteGridlessADP from, DiscreteGridlessADP to)
	{
		to.setX0(new ArrayRealVector(from.getX0()));
		to.getGridlessData().power_0 = new ArrayRealVector(from.getGridlessData().power_0);
		to.getGridlessData().v_0 = new ArrayRealVector(from.getGridlessData().v_0);
		if(from.approxJ != null)
		{
			to.approxJ = new NonParametricApproximator[from.approxJ.length];
			for (int i = 0; i < to.approxJ.length; i++)
			{
				to.approxJ[i] = from.approxJ[i];
			}
		}
		to.costFunction = from.costFunction;
		to.getGridlessData().deltaV_external = new RealVector[from.getGridlessData().deltaV_external.length];
		for (int i = 0; i < to.getGridlessData().deltaV_external.length; i++)
		{
			to.getGridlessData().deltaV_external[i] = new ArrayRealVector(from.getGridlessData().deltaV_external[i]);
		}
		to.getGridlessData().dgCount = from.getGridlessData().dgCount;
		to.gamma = from.gamma;
		to.getGridlessData().loadCount = from.getGridlessData().loadCount;
		to.getGridlessData().meanP_DG = from.getGridlessData().meanP_DG;
		to.getGridlessData().meanP_L = from.getGridlessData().meanP_L;
		to.getGridlessData().meanQ_L = from.getGridlessData().meanQ_L;
		to.N_r = from.N_r;
		to.getGridlessData().sdP_DG = from.getGridlessData().sdP_DG;
		to.getGridlessData().sdP_L = from.getGridlessData().sdP_L;
		to.getGridlessData().sdQ_L = from.getGridlessData().sdQ_L;
		to.getGridlessData().sensitivities = new Array2DRowRealMatrix(from.getGridlessData().sensitivities.getData());
		to.getGridlessData().storageChargeRate = from.getGridlessData().storageChargeRate;
		to.getGridlessData().storageCount = from.getGridlessData().storageCount;
		to.getGridlessData().storageDischargeRate = from.getGridlessData().storageDischargeRate;
		to.getGridlessData().storageMaxCapacity = from.getGridlessData().storageMaxCapacity;
		to.T = from.T;
		
		return to;
	}
}
