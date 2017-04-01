package ellipsis.energy.util;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.QuadraticEstimator;
import ellipsis.util.SampleBasedEstimator;

public class LCGridlessADP extends LinearConvexADP implements GridlessADP
{
	protected GridlessData gridlessData;
	
	public LCGridlessADP()
	{
		gridlessData = new GridlessData(this);
	}
	
	private Random rand = new Random(0);

	@Override
	protected RealVector randomControl(int t, RealVector x_t, RealVector w_t)
	{
		RealVector u = new ArrayRealVector(controlDimension); // TODO this is only for active power control of DG and storage.
		
		int j = 0;
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			// DG \in [0, max]:
			u.setEntry(j, rand.nextDouble()*gridlessData.meanP_DG[t][i]); // FIXME change this back - it works better with it, so not sure...
																	      //       After all, this is for the purpose of learning the cost
			                                                              //       curve, which itself is not limited by constraints.
//			u.setEntry(j, rand.nextDouble()*gridlessData.rating_DG[i]);
			++j;
		}
		
		// UNTESTED!
		for(int i = 0; i < gridlessData.storageCount; ++i)
		{
			double soc = x_t.getEntry(i);
			
			// Storage \in [-max, 0, max] subject to capacity constraints:
			double maxDischarge = soc; // always +ve
			double maxCharge = soc - gridlessData.storageMaxCapacity[i]; // always -ve
			maxDischarge = Math.min(maxDischarge, gridlessData.storageDischargeRate[i]);
			maxCharge = Math.max(maxCharge, gridlessData.storageChargeRate[i]);
			
			u.setEntry(j, rand.nextDouble()*(maxDischarge - maxCharge) + maxCharge);
			
			++j;
		}
		
		return u;
	}

	@Override
	public ADP makeNew()
	{
		return new LCGridlessADP();
	}

	@Override
	public SampleBasedEstimator newEstimator(double[] defaultBandwidth)
	{
		QuadraticEstimator qe = new QuadraticEstimator();
		return qe;
	}

	protected int powerDimension;
	protected int unitCount;
	protected int controlDimension;
	protected int stateDimension;
	protected int voltageDimension;
	protected int dgPowerOffset;
	protected int loadPowerOffset;
	/**
	 * Sets some variables that can be used instead of method calls to speed up processing.
	 * Must be called before training or any other operations.
	 */
	public void finaliseData()
	{
		unitCount = unitCount();
		powerDimension = powerDimension();
		controlDimension = controlDimension();
		stateDimension = stateDimension();
		voltageDimension = voltageDimension();
		dgPowerOffset = dgPowerOffset();
		loadPowerOffset = loadPowerOffset();
	}
	
	
	//// GridlessADP Implemention ////

	@Override
	public GridlessData getGridlessData()
	{
		return gridlessData;
	}

	@Override
	public int stateDimension()
	{
		return gridlessData.dgCount;
	}

	@Override
	public int voltageDimension()
	{
		return 2*unitCount; // arg, abs
	}

	@Override
	public int powerDimension()
	{
		return 2*unitCount; // P, Q
	}

	@Override
	public int controlDimension()
	{
		return 2*(gridlessData.dgCount+gridlessData.storageCount); // P & Q
	}

	@Override
	public int unitCount()
	{
		return gridlessData.dgCount+gridlessData.storageCount+gridlessData.loadCount;
	}

	@Override
	public void setPower0(int i, double p)
	{
		gridlessData.power_0.setEntry(i, p);
	}

	@Override
	public void setV0(int i, double v)
	{
		gridlessData.v_0.setEntry(i, v);
	}
	
	@Override
	public RealVector getV0()
	{
		return gridlessData.v_0;
	}

	@Override
	public void setPower0(RealVector p0)
	{
		gridlessData.power_0 = p0;
	}

	@Override
	public void setV0(RealVector v0)
	{
		gridlessData.v_0 = v0;
	}

	@Override
	public void setDeltaVExternal(int t, int i, double deltaV)
	{
		gridlessData.deltaV_external[t].setEntry(i, deltaV);
	}
	
	@Override
	public void setDeltaVExternal(RealVector[] deltaVExt)
	{
		gridlessData.deltaV_external = deltaVExt;
	}
	
	@Override
	public RealVector[] getDeltaVExternal()
	{
		return gridlessData.deltaV_external;
	}
	
	/**
	 * [P Q] = [P_DG P_L Q_DG Q_L] = [u(0-dgCount) w(0-loadCount) u(dgCount-2dgCount) w(loadCount-2loadCount)]
	 * @param u
	 * @return
	 */
	public RealVector powerFromControlAndNoise(RealVector u, RealVector w)
	{
		double[] pqArray = new double[powerDimension];
		double[] uArray = u.toArray();
		double[] wArray = w.toArray();
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			pqArray[i+dgPowerOffset] = uArray[i]; // P_DG
			pqArray[i+dgPowerOffset+unitCount] = uArray[i+controlDimension/2]; // Q_DG
		}
		// TODO include storage
		for(int i = 0; i < gridlessData.loadCount; ++i)
		{
			pqArray[i+loadPowerOffset] = wArray[i+gridlessData.dgCount]; // P_L
			pqArray[i+loadPowerOffset+unitCount] = wArray[i+gridlessData.dgCount+gridlessData.loadCount]; // Q_L
		}
		RealVector power = new ArrayRealVector(pqArray);
		
		return power;
	}
	
	public int dgPowerOffset()
	{
		return 0;
	}
	
	public int storagePowerOffset()
	{
		return gridlessData.dgCount;
	}
	
	public int loadPowerOffset()
	{
		return gridlessData.dgCount+gridlessData.storageCount;
	}
	
	public RealVector voltagesFromControlAndNoise(int t, RealVector u, RealVector w)
	{
		RealVector power = powerFromControlAndNoise(u, w);
		return voltageFromPower(t, power);
	}

	public RealVector voltageFromPower(int t, RealVector power)
	{
		RealVector deltaS = power.subtract(gridlessData.power_0);

		RealVector deltaV = gridlessData.sensitivities.operate(deltaS);
		RealVector v = gridlessData.v_0.add(deltaV).add(gridlessData.deltaV_external[t]);

		return v;
	}
	
	
	//// MultiBus impl ////

	@Override
	public List<Integer> getBusIndeces()
	{
		return gridlessData.B;
	}

	@Override
	public void setBusIndeces(List<Integer> B)
	{
		gridlessData.B = B;
	}
}