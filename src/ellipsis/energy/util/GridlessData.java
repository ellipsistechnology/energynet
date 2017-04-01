package ellipsis.energy.util;

import static ellipsis.util.VectorHelper.vector;

import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class GridlessData
{
	// Parameters:
	public RealVector power_0; // [P_DG P_S P_L Q_DG Q_S Q_L] @ t=0
	public Complex slackAdmittance;
	public RealMatrix slackSensitivities;
	public RealMatrix sensitivities;
	public RealVector v_0;
	public RealVector[] deltaV_external; // [time]
	
	public int storageCount = 1;
	public int dgCount = 1;
	public int loadCount = 1;
	
	public double[][] meanP_L; // [time][unit]
	public double[][] meanQ_L; // [time][unit]
	public double[][] sdP_L; // [time][unit]
	public double[][] sdQ_L; // [time][unit]
	
	public double[][] meanP_DG; // [time][unit]
	public double[][] sdP_DG; // [time][unit]
	public double rating_DG[]; // [item]

	public double[] storageDischargeRate = new double[]{1}; // [unit]
	public double[] storageChargeRate = new double[]{-1}; // [unit]
	public double[] storageMaxCapacity = new double[]{2};
	
	// Map of original names to indices
	// (only relevant if constructed from a Grid):
	public Map<Integer, String> stateNames;
	public Map<Integer, String> controlNames;
	public Map<Integer, String> noiseNames;
	public Map<Integer, String> voltageNames;
	
	public List<Integer> B;
	
	private int T;
	
	public GridlessData(ADP adp)
	{
		meanP_L = new double[adp.T][];
		meanQ_L = new double[adp.T][];
		sdP_L = new double[adp.T][];
		sdQ_L = new double[adp.T][];
		meanP_DG = new double[adp.T][];
		sdP_DG = new double[adp.T][];

		for (int t = 0; t < adp.T; t++)
		{
			meanP_L[t] = new double[]{-1};
			meanQ_L[t] = new double[]{0};
			sdP_L[t] = new double[]{1};
			sdQ_L[t] = new double[]{0};
			meanP_DG[t] = new double[]{1};
			sdP_DG[t] = new double[]{0.1};
		}
		
		setDeltaVExternal(vector(6, 0.0));
		
		this.T = adp.T;
	}

	/**
	 * Sets the scheduled standard deviation DG real power output for all DG assuming that all DG have the same schedule.
	 * @param sdP_DG The schedule to apply to each unit.
	 */
	public void setsdP_DG(double[] sdP_DG)
	{
		this.sdP_DG = new double[T][];
		fillSchedules(sdP_DG, this.sdP_DG, dgCount);
	}

	/**
	 * Sets the scheduled mean DG real power output for all DG assuming that all DG have the same schedule.
	 * @param meanP_DG The schedule to apply to each unit.
	 */
	public void setMeanP_DG(double[] meanP_DG)
	{
		this.meanP_DG = new double[T][];
		fillSchedules(meanP_DG, this.meanP_DG, dgCount);
	}

	/**
	 * Sets the scheduled standard deviation load real power output for all DG assuming that all load have the same schedule.
	 * @param sdP_DG The schedule to apply to each unit.
	 */
	public void setsdP_L(double[] sdP_L)
	{
		this.sdP_L = new double[T][];
		fillSchedules(sdP_L, this.sdP_L, loadCount);
	}

	/**
	 * Sets the scheduled mean load real power output for all DG assuming that all load have the same schedule.
	 * @param meanP_DG The schedule to apply to each unit.
	 */
	public void setMeanP_L(double[] meanP_L)
	{
		this.meanP_L = new double[T][];
		fillSchedules(meanP_L, this.meanP_L, loadCount);
	}

	/**
	 * Sets the scheduled standard deviation load reactive power output for all DG assuming that all load have the same schedule.
	 * @param sdP_DG The schedule to apply to each unit.
	 */
	public void setsdQ_L(double[] sdQ_L)
	{
		this.sdQ_L = new double[T][];
		fillSchedules(sdQ_L, this.sdQ_L, loadCount);
	}

	/**
	 * Sets the scheduled mean load reactive power output for all DG assuming that all load have the same schedule.
	 * @param meanP_DG The schedule to apply to each unit.
	 */
	public void setMeanQ_L(double[] meanQ_L)
	{
		this.meanQ_L = new double[T][];
		fillSchedules(meanQ_L, this.meanQ_L, loadCount);
	}

	private void fillSchedules(double[] schedule, double[][] target, int unitCount)
	{
		for(int t = 0; t < T; ++t)
		{
			target[t] = new double[unitCount];
			for(int i = 0; i < unitCount; ++i)
				target[t][i] = schedule[t];
		}
	}
	
	public void setDeltaVExternal(RealVector template)
	{
		deltaV_external = new RealVector[T+1];
		for (int i = 0; i <= T; i++)
		{
			deltaV_external[i] = new ArrayRealVector(template);
		}
	}

	private double[][] sensitivitiesArray;
	public double[][] sensitivitiesArray()
	{
		if(sensitivitiesArray == null)
		{
			sensitivitiesArray = sensitivities.getData();
		}
		return sensitivitiesArray;
	}
}