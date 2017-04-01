package ellipsis.energy.util;

import static ellipsis.util.VectorHelper.vector;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/**
 * @author bmillar
 *
 */
public class DiscreteGridlessADP extends DiscreteADP implements GridlessADP
{
	private GridlessData gridlessData;
	
	public DiscreteGridlessADP()
	{
		gridlessData = new GridlessData(this);
	}

	@Override
	public RealVector randomVariations(int t)
	{
		RealVector w_t = new ArrayRealVector(gridlessData.dgCount+2*gridlessData.loadCount);
		
		// Random shift in DG output:
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			double P = rand.nextGaussian()*gridlessData.sdP_DG[t][i]; // zero mean
			w_t.setEntry(i, P);
		}
		
		// Loads:
		for(int i = 0; i < gridlessData.loadCount; ++i)
		{
			int j = gridlessData.dgCount + i;
			double sdP = rand.nextGaussian()*gridlessData.sdP_L[t][i];
			double sdQ = rand.nextGaussian()*gridlessData.sdQ_L[t][i];
			double P = gridlessData.meanP_L[t][i] + sdP;
			double Q = gridlessData.meanQ_L[t][i] + sdQ;
			w_t.setEntry(j, P);
			w_t.setEntry(j+gridlessData.loadCount, Q);
		}
		
		return w_t;
	}

	@Override
	public RealVector wZero(int t)
	{
		RealVector wZero = vector(gridlessData.dgCount+2*gridlessData.loadCount, 0.0);
		if(t == T)
			return wZero;
		int wDimension = wZero.getDimension();
		for(int i = 0; i < gridlessData.loadCount; ++i)
		{
			wZero.setEntry(wDimension-2*gridlessData.loadCount+i, this.gridlessData.meanP_L[t][i]);
		}
		for(int i = 0; i < gridlessData.loadCount; ++i)
		{
			wZero.setEntry(wDimension-gridlessData.loadCount+i, this.gridlessData.meanQ_L[t][i]);
		}
		return wZero;
	}

	@Override
	public RealVector uZero()
	{
		return vector(gridlessData.dgCount+gridlessData.storageCount, 0);
	}

	@Override
	public Set<RealVector> enumerateControls(int t, RealVector x_t)
	{
		Set<RealVector> U_t = new HashSet<>();
			
		double[][] options = new double[gridlessData.dgCount + gridlessData.storageCount][];
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			// DG \in [0, max]:
			options[i] = new double[]{0.0, gridlessData.meanP_DG[t][i]};
		}
		
		for(int i = 0; i < gridlessData.storageCount; ++i)
		{
			double soc = x_t.getEntry(i);
			
			// Storage \in [-max, 0, max] subject to capacity constraints:
			if(soc == 0) // Empty.
				options[gridlessData.dgCount+i] = new double[]{0.0, gridlessData.storageChargeRate[i]};
			else if(soc - gridlessData.storageDischargeRate[i] < 0) // Nearly empty.
				options[gridlessData.dgCount+i] = new double[]{soc, 0.0, gridlessData.storageChargeRate[i]};
			else if(soc == gridlessData.storageMaxCapacity[i]) // Full.
				options[gridlessData.dgCount+i] = new double[]{gridlessData.storageDischargeRate[i], 0.0};
			else if(soc - gridlessData.storageChargeRate[i] > gridlessData.storageMaxCapacity[i]) // Almost full.
				options[gridlessData.dgCount+i] = new double[]{gridlessData.storageDischargeRate[i], 0.0, soc-gridlessData.storageMaxCapacity[i]};
			else
				options[gridlessData.dgCount+i] = new double[]{gridlessData.storageDischargeRate[i], 0.0, gridlessData.storageChargeRate[i]};
		}
		
		RealVector u = uZero();
		int i = gridlessData.dgCount+gridlessData.storageCount-1;
		enumerateControls(U_t, options, u, i);
		return U_t;
	}

	private void enumerateControls(Set<RealVector> U, double[][] options, RealVector u, int i)
	{
		for(int j = 0; j < options[i].length; ++j)
		{
			u.setEntry(i, options[i][j]);
			
			// If iterating over the last unit's options then add this control combination:
			if(i == 0)
			{
				U.add(new ArrayRealVector(u));
			}
			else
			{
				enumerateControls(U, options, u, i-1);
			}
		}
	}

	@Override
	public void prepareNoise(RealVector u, RealVector w)
	{
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			if(u.getEntry(i) == 0)
				w.setEntry(i, 0);
		}
	}

	@Override
	public boolean breached(int t, RealVector u, RealVector w)
	{
		RealVector v = voltagesFromControlAndNoise(t, u, w);
		RealVector v_abs = v.getSubVector(v.getDimension()/2, v.getDimension()/2);
		boolean breached = v_abs.getMaxValue() > 1.05 || v_abs.getMinValue() < 0.95;
		return breached;
	}

	/**
	 *   x_t = [ SOC_t ]    u_t = [ P_{DG,t} ]    w_t = [ \Delta P_{DG,t} ]
		                          [ P_{S,t}  ]          [ P_{L,t}         ]
		                                                [ Q_{L,t}         ]

		 x_{t+1} = f(x_t, u_t, w_t) = Ax_t + Bu_t + Cw_t
		                            = x_t + P_{S,t}
		 
		 i.e. Noise only affects the constraints and not the next state.

	 * @param x
	 * @param u
	 * @param w
	 * @return
	 */
	@Override
	public RealVector f(int t, RealVector x_t, RealVector u_t, RealVector w_t)
	{
		// Split x_t into sub-vectors:
		int i = 0;
		RealVector SOC_t = x_t.getSubVector(i, gridlessData.storageCount);
		
		// Split u_t into sub-vectors:
		i = 0;
		i += gridlessData.dgCount;
		RealVector P_S_t1 = u_t.getSubVector(i, gridlessData.storageCount);
		
		// Next SOC:
		RealVector SOC_next = SOC_t.subtract(P_S_t1);
		
		// Next x_t:
		RealVector x_next = SOC_next;
		
		return x_next;
	}

	@Override
	public ADP makeNew()
	{
		return new DiscreteGridlessADP();
	}

	@Override
	public GridlessData getGridlessData()
	{
		return gridlessData;
	}

	@Override
	public int stateDimension()
	{
		return gridlessData.storageCount;
	}

	@Override
	public int voltageDimension()
	{
		return 2*unitCount(); // arg, abs
	}

	@Override
	public int powerDimension()
	{
		return 2*unitCount(); // P, Q
	}

	@Override
	public int controlDimension()
	{
		return gridlessData.dgCount+gridlessData.storageCount; // P
	}

	@Override
	public int unitCount()
	{
		return gridlessData.dgCount+gridlessData.storageCount+gridlessData.loadCount;
	}
	
	@Override
	public int dgPowerOffset()
	{
		return 0;
	}

	@Override
	public int storagePowerOffset()
	{
		return gridlessData.dgCount;
	}

	@Override
	public int loadPowerOffset()
	{
		return gridlessData.dgCount+gridlessData.storageCount;
	}
	
	/**
	 * Gives a vector of the form [P_DG P_S P_L Q_DG Q_S Q_L].
	 * @param u
	 * @return
	 */
	public RealVector powerFromControlAndNoise(RealVector u, RealVector w)
	{
		RealVector pq = 
				u.append( // P_DG P_S
				w.getSubVector(gridlessData.dgCount, gridlessData.loadCount)).append( // P_L
				vector(u.getDimension(), 0)).append( // Q_DG Q_S
				w.getSubVector(gridlessData.dgCount+gridlessData.loadCount, gridlessData.loadCount)); // Q_L
		
		// Add DG power output noise terms:
		for(int i = 0; i < gridlessData.dgCount; ++i)
		{
			double P_DG = u.getEntry(i)+w.getEntry(i);
			pq.setEntry(i, P_DG);
		}
		
		return pq;
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
	
	
	//// LocalController implementation ////

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
	public void setDeltaVExternal(RealVector[] deltaVExt)
	{
		gridlessData.deltaV_external = deltaVExt;
	}

	@Override
	public void setDeltaVExternal(int t, int i, double deltaV)
	{
		gridlessData.deltaV_external[t].setEntry(i, deltaV);
	}
	
	@Override
	public RealVector[] getDeltaVExternal()
	{
		return gridlessData.deltaV_external;
	}

	@Override
	public void setV0(RealVector v_0)
	{
		gridlessData.v_0 = v_0;
	}

	@Override
	public void setPower0(RealVector power_0)
	{
		gridlessData.power_0 = power_0;
	}

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
