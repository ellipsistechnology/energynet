package ellipsis.energy.test;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;

public class TransformerTest
{
	public static void main(String[] args)
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		
		double oldM = 0.9;
		for(double m = 0.9; m <= 1.1; m += 0.01)
		{
			System.out.println("m "+oldM+" => "+m);
			AnalysisResults results = analyse(grid);
			testWithRatio(m, grid, results);
			oldM = m;
		}
		
		for(double m = 0.9; m <= 1.1; m += 0.01)
		{
			System.out.println("m 1.0 => "+m);
			grid.xfm650.setT(1.0);
			AnalysisResults results = analyse(grid);
			testWithRatio(m, grid, results);
		}
	}

	private static AnalysisResults analyse(IEEE13BusGrid grid)
	{
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
		lfa.setBasePower(IEEE13BusGrid.BASE_POWER);
		lfa.setBaseVoltage(IEEE13BusGrid.BASE_VOLTAGE);
		AnalysisResults results = lfa.analyse();
		return results;
	}

	@SuppressWarnings("unused")
	private static void testWithRatio(double tapRatio, IEEE13BusGrid grid, AnalysisResults results)
	{
		// Init:
		RealMatrix sensitivities = results.sensitivities();
		RealVector v0 = results.voltages();
		int dimension = v0.getDimension()/2;
		String bus650xfmName = grid.bus650xfm.getName();
		Complex xfmrShunt0 = grid.xfm650.getToShunt().getAdmittance().multiply(IEEE13BusGrid.BASE_IMPEDANCE);
		double tapRatio0 = grid.xfm650.voltageRatio();
		grid.xfm650.setT(tapRatio);
		
		// Baseline with new results:
		RealVector v_true;
		Complex newV;
		Complex newI;
		Complex newS;
		Complex trueS0 = results.getBusPower(bus650xfmName);
		Complex truev0 = results.getBusVoltage(bus650xfmName);
		Complex trueI0 = results.busCurrent(grid.bus650xfm);
		{
			AnalysisResults newResults = analyse(grid);
			v_true = newResults.voltages();
			
			// Results for comparison with pi method:
			newV = newResults.getBusVoltage(bus650xfmName);
			newI = newResults.busCurrent(grid.bus650xfm);
			newS = newResults.getBusPower(bus650xfmName);
		}
		
		// Transformer as PI equivalent with fixed sensitivities:
		RealVector v_pi;
		{
			Complex xfmrShunt = grid.xfm650.getToShunt().getAdmittance().multiply(IEEE13BusGrid.BASE_IMPEDANCE);
			int index650xfm = results.getBusNumbers().get(bus650xfmName);
			if(index650xfm > results.getSlackIndex())
				--index650xfm;
			int deltaOffset = dimension;
			
			// Get power before tap change:
			double v0_abs = v0.getEntry(deltaOffset+index650xfm);
			double v0_arg = v0.getEntry(index650xfm);
			Complex v0_650xfm = ComplexUtils.polar2Complex(v0_abs, v0_arg);
			Complex s0 = xfmrShunt0.conjugate().multiply(v0_650xfm.abs()*v0_650xfm.abs()).negate();
			
			// Power after tap change:
			Complex v_650xfm = new Complex(tapRatio); // arg is irrelevant since only abs is used below
			Complex s = xfmrShunt.conjugate().multiply(v_650xfm.abs()*v_650xfm.abs()).negate(); // s = |v|^2y*
			Complex i = v_650xfm.multiply(xfmrShunt).negate(); // For testing/reference only.
			
			// Get change in power:
			RealVector deltaPQ = new ArrayRealVector(v0.getDimension());
			deltaPQ.setEntry(index650xfm, s.getReal()-s0.getReal());
			deltaPQ.setEntry(index650xfm+deltaOffset, s.getImaginary()-s0.getImaginary());
			
			// Network change in voltages:
			RealVector deltaV = sensitivities.operate(deltaPQ);
			v_pi = v0.add(deltaV);
		}
		
		// Voltages all shifted equally according to ratio:
		RealVector v_m;
		{
			// Convert to complex:
			Complex[] v = new Complex[dimension];
			for (int i = 0; i < v.length; i++)
			{
				v[i] = ComplexUtils.polar2Complex(v0.getEntry(i+dimension), v0.getEntry(i));
				v[i] = v[i].add(tapRatio-tapRatio0);
			}
			
			// Convert back:
			v_m = new ArrayRealVector(dimension*2);
			for (int i = 0; i < v.length; i++)
			{
				v_m.setEntry(i+dimension, v[i].abs());
				v_m.setEntry(i, v[i].getArgument());
			}
		}
		
		// Display results:
		System.out.println(",abs,,,arg");
		System.out.println("bus,v_true,v_pi,v_m,v_true,v_pi,v_m");
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(i);
			System.out.print(",");
			System.out.print(v_true.getEntry(i+dimension));
			System.out.print(",");
			System.out.print(v_pi.getEntry(i+dimension));
			System.out.print(",");
			System.out.print(v_m.getEntry(i+dimension));
			System.out.print(",");
			System.out.print(v_true.getEntry(i));
			System.out.print(",");
			System.out.print(v_pi.getEntry(i));
			System.out.print(",");
			System.out.print(v_m.getEntry(i));
			System.out.print(",");
			System.out.println();
		}
	}
}