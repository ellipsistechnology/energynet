package ellipsis.energy.casestudy;

import java.io.PrintStream;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessCentralCoordinator;
import ellipsis.util.EigenValueHelper;

/**
 * Logs various information useful for analysing and reporting 
 * case study results.
 * @author bmillar
 *
 */
public class TestCaseLogger
{
	public static PrintStream out = System.out;

	public static void logDataWithNamesAndTimes(int T, RealVector[] timeData, Map<Integer, String> names)
	{
		// Header (times):
		out.print("-,");
		for(int t = 0; t < T; ++t)
		{
			out.print(t);
			out.print(',');
		}
		out.println();
		
		// Rows:
		for (Integer i : names.keySet())
		{
			// First column (unit name):
			out.print(names.get(i));
			out.print(',');
			for(int t = 0; t < T; ++t)
			{
				if(timeData != null && timeData.length > t && timeData[t] != null && timeData[t].getDimension() > i)
					out.print(timeData[t].getEntry(i));
				out.print(',');
			}
			out.println();
		}
		
		int dimension = names.size();
		if(timeData[0].getDimension() == dimension*2)
		{
			for (Integer i : names.keySet())
			{
				// First column (unit name):
				out.print(names.get(i));
				out.print(',');
				for(int t = 0; t < T; ++t)
				{
					if(timeData != null && timeData.length > t && timeData[t] != null)
						out.print(timeData[t].getEntry(i+dimension));
					out.print(',');
				}
				out.println();
			}
		}
	}
	
	public static void logControlStateNoiseVoltage(RealVector[] x, RealVector[] u,
			RealVector[] w, Map<Integer, String> controlNames, Map<Integer, String> stateNames, Map<Integer, String> noiseNames, Map<Integer, String> voltageNames,
			GridlessADP adp)
	{
		RealVector[] voltages = new RealVector[adp.getT()+1];
		if(u != null)
		{
			for(int t = 0; t < adp.getT(); ++t)
			{
				voltages[t] = adp.voltagesFromControlAndNoise(t, u[t], adp.wZero(t));
				voltages[t] = voltages[t].getSubVector(voltages[t].getDimension()/2, voltages[t].getDimension()/2); // Get magnitude only.
			}
		}
		logControlStateNoiseVoltage(x, u, w, voltages, controlNames, stateNames, noiseNames, voltageNames, adp);
	}

	public static void logControlStateNoiseVoltage(RealVector[] x, RealVector[] u, RealVector[] w, RealVector[] v, 
			Map<Integer, String> controlNames, Map<Integer, String> stateNames, Map<Integer, String> noiseNames, Map<Integer, String> voltageNames,
			GridlessADP adp)
	{
		// Log controls:
		out.println("Controls:");
		logDataWithNamesAndTimes(adp.getT(), u, controlNames);
		
		// Log states:
		out.println("States:");
		logDataWithNamesAndTimes(adp.getT(), x, stateNames);
		
		// Log noise:
		out.println("Noise:");
		logDataWithNamesAndTimes(adp.getT(), w, noiseNames);
		
		// Log voltages:
		out.println("Voltages:");
		logDataWithNamesAndTimes(adp.getT(), v, voltageNames);
	}

	public static void logCCConfig(GridlessCentralCoordinator cc, int iterations, double probabilityOfUpdate, double bw)
	{
		RealMatrix sensitivities = ((GridlessADP)cc.localControllerTemplate).getGridlessData().sensitivities;
		Complex[] eigens = EigenValueHelper.eigenValues(sensitivities);
		out.println(
				"N_r="+((GridlessADP)cc.localControllerTemplate).getNr()+
				",C_max="+cc.C_max+
				",U_max="+cc.U_max+
				",alpha="+cc.alpha+
				",dgCount="+cc.gridlessADP().getGridlessData().dgCount+
				",storageCount="+cc.gridlessADP().getGridlessData().storageCount+
				",loadCount="+cc.gridlessADP().getGridlessData().loadCount+
				",iterations="+iterations+
				",probabilityOfUpdate="+probabilityOfUpdate+
				",bw="+bw+
				",passTimeRate="+cc.passTimeRate+
				",ignoreExternal="+cc.ignoreExternal+
				",largesEigen="+EigenValueHelper.maxValue(eigens));
	}
}