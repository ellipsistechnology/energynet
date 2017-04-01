package ellipsis.energy.util;

import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class GridlessCapacitorRegulator implements LocalController
{
	private int T;
	private double voltageReference;
	private RealVector sensitivities; // d|v|/dQ
	private RealVector power0;
	private double Qmax;
	private int unitCount;
	private RealVector v0;
	private RealVector[] deltaVExternal;
	private int busIndex;
	private List<Integer> B;

	@Override
	public RealVector[] schedule(double[] ctgs)
	{
		RealVector[] schedule = new RealVector[T];
		for(int t = 0; t < T; ++t)
		{
			ArrayRealVector u_on = new ArrayRealVector(new double[]{Qmax});
			double costOn = cost(u_on);
			ArrayRealVector u_off = new ArrayRealVector(new double[]{0.0});
			double costOff = cost(u_off);
			if(costOn < costOff)
				schedule[t] = u_on;
			else
				schedule[t] = u_off;
		}
		return schedule;
	}

	@Override
	public void train(double[][] defaultBandwidth)
	{
		// Nothing to do.
	}

	@Override
	public RealVector uZero()
	{
		return new ArrayRealVector(new double[]{0.0});
	}
	
	protected double cost(RealVector u)
	{
		double sum = 0;
		for(int t = 0; t < T; ++t)
		{
			RealVector v = voltageMagnitudeFromControl(t, u);
			RealVector vDiff = v.mapSubtract(voltageReference);
			sum += vDiff.dotProduct(vDiff);
		}
		return sum;
	}
	
	private RealVector voltageMagnitudeFromControl(int t, RealVector u)
	{
		double q = u.getEntry(0);
		double deltaq = q - power0.getEntry(busIndex+unitCount);
		RealVector deltaV = sensitivities.mapMultiply(deltaq);
		return v0.add(deltaV).add(deltaVExternal[t]);
	}
	
	
	//// Local Controller Implementation ////

	@Override
	public int unitCount()
	{
		return unitCount;
	}

	@Override
	public void setPower0(RealVector pq0)
	{
		power0 = pq0;
	}

	@Override
	public void setPower0(int i, double p)
	{
		power0.setEntry(i, p);
	}

	@Override
	public void setV0(RealVector v0)
	{
		this.v0 = v0;
	}

	@Override
	public void setV0(int i, double v)
	{
		v0.setEntry(i, v);
	}
	
	@Override
	public RealVector getV0()
	{
		return v0;
	}

	@Override
	public void setDeltaVExternal(RealVector[] deltaVExt)
	{
		this.deltaVExternal = deltaVExt;
	}

	@Override
	public void setDeltaVExternal(int t, int i, double deltaV)
	{
		deltaVExternal[t].setEntry(i, deltaV);
	}
	
	@Override
	public RealVector[] getDeltaVExternal()
	{
		return deltaVExternal;
	}

	@Override
	public int getT()
	{
		return T;
	}

	
	//// Accessors ////
	
	public void setT(int t)
	{
		T = t;
	}

	public double getVoltageReference()
	{
		return voltageReference;
	}

	public void setVoltageReference(double voltageReference)
	{
		this.voltageReference = voltageReference;
	}

	public RealVector getSensitivities()
	{
		return sensitivities;
	}

	public void setSensitivities(RealVector sensitivities)
	{
		this.sensitivities = sensitivities;
	}

	public double getQmax()
	{
		return Qmax;
	}

	public void setQmax(double qmax)
	{
		Qmax = qmax;
	}

	public int getBusIndex()
	{
		return busIndex;
	}

	public void setBusIndex(int busIndex)
	{
		this.busIndex = busIndex;
	}
	
	public void setUnitCount(int unitCount)
	{
		this.unitCount = unitCount;
	}

	@Override
	public List<Integer> getBusIndeces()
	{
		return B;
	}

	@Override
	public void setBusIndeces(List<Integer> B)
	{
		this.B = B;
	}
}
