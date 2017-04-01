package ellipsis.energy.util;

import static ellipsis.util.VectorHelper.sum;
import static ellipsis.util.VectorHelper.vector;

import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * Transformer with subset including all loads, and externals including others such as DG.
 * @author bmillar
 *
 */
public class GridlessTransformerRegulator implements LocalController
{
	// Parameters:
//	private int busCount;
	private int T;
	private double voltageReference;
	private int minStep, maxStep;
	private double stepSize;
    private int tapChangeDelay;
    private Complex admittance;
	private List<Integer> B;
	private RealMatrix sensistivities;
	private RealVector[] loadPower;
	
	// Variables:
    private RealVector v0;
	private RealVector power_0;
	private double ratioSchedule[];
	private RealVector[] deltaVExternal;

	private double ratio(int step)
	{
		return 1+step*stepSize;
	}

	public double cost(int t, double m)
	{
		RealVector v[] = calculateVoltages();
		RealVector vDiff = v[t].mapAdd(m-ratioSchedule[t]-voltageReference);
		return vDiff.dotProduct(vDiff);
	}
	
	private RealVector[] calculateVoltages()
	{
		RealVector v[] = new RealVector[T];
		for(int t = 0; t < T; ++t)
		{
			// Initial voltages:
//			RealVector v_abs = v0.getSubVector(v0.getDimension()/2, v0.getDimension()/2);
//			v[t] = v_abs;
			
			int absOffset = v0.getDimension()/2;
			if(deltaVExternal == null || deltaVExternal[t] == null)
				v[t] = localBaseVoltage(t).getSubVector(absOffset, absOffset);
			else
				v[t] = localBaseVoltage(t).add(deltaVExternal[t]).getSubVector(absOffset, absOffset); // abs only
		}
		
		return v;
	}

	private RealVector localBaseVoltage(int t)
	{
		// Make power vector from OLTC grid-side power and load powers:
		Complex s = power(t);
		loadPower[t].setEntry(0, s.getReal());
		loadPower[t].setEntry(loadPower[t].getDimension()/2, s.getImaginary());
		
		// Get change in voltage from change in power:
		RealVector deltaPQ = loadPower[t].subtract(power_0);
		RealVector deltaV = sensistivities.operate(deltaPQ);
		
		return v0.add(deltaV);
	}

	public LocalController copy()
	{
		GridlessTransformerRegulator gtr = new GridlessTransformerRegulator();
//		gtr.setBusCount(busCount);
		gtr.setMaxStep(maxStep);
		gtr.setMinStep(minStep);
		double[] ratioScheduleCopy = copyRatioSchedule();
		gtr.setRatioSchedule(ratioScheduleCopy);
		gtr.setStepSize(stepSize);
		gtr.setT(T);
		gtr.setTapChangeDelay(tapChangeDelay);
		gtr.setVoltageReference(voltageReference);
		gtr.setAdmittance(admittance);
		if(v0 != null)
			gtr.setV0(v0.copy());
		
		return gtr;
	}

	public double[] copyRatioSchedule()
	{
		double[] ratioScheduleCopy = new double[ratioSchedule.length];
		System.arraycopy(ratioSchedule, 0, ratioScheduleCopy, 0, ratioSchedule.length);
		return ratioScheduleCopy;
	}
	

	//// Tap changing ////
    
    // Variables:
    private int currentTime;
    private int nextTapChangeTime;
    private double targetTapPosition;

	public void applySchedule(double[] schedule)
    {
    	// Update the schedule for time 0 with appropriate delay(s):
    	if(schedule[0] != targetTapPosition)
    	{
    		currentTime = 0;
	    	targetTapPosition = schedule[0];
	    	if(targetTapPosition != ratioSchedule[0])
	    	{
	    		nextTapChangeTime = tapChangeDelay;
	    	}
    	}
    	
    	// Copy the rest of the schedule:
    	for(int t = 1; t < schedule.length; ++t)
    		ratioSchedule[t] = schedule[t];
    }
	
    public void passTime(int timeDelta)
    {
        if(currentTime+timeDelta >= nextTapChangeTime)
        {
	        if(targetTapPosition > ratioSchedule[0])
	            ratioSchedule[0] += stepSize;
	        else
	            ratioSchedule[0] -= stepSize;
	        
	        if(ratioSchedule[0] == targetTapPosition)
	        {
	        	nextTapChangeTime = Integer.MAX_VALUE;
	        }
	        else
	        {
	        	nextTapChangeTime += tapChangeDelay;
	        }
        }
        
        currentTime += timeDelta;
    }
	
	
	//// LocalController implementation ////

    public double[] mSchedule; // For testing only.
	@Override
	public RealVector[] schedule(double[] ctgs)
	{
		RealVector v[] = calculateVoltages();
		mSchedule = new double[T];
		RealVector[] pqSchedule = new RealVector[T];
		for(int t = 0; t < T; ++t)
		{
			// Theoretically optimal tap ratio m:
			RealVector vDiff = v[t].mapAdd(-voltageReference);
			double m = -(1.0/v[t].getDimension())*sum(vDiff) + ratioSchedule[t];
			
			// Allowable tap positions around m:
			int tapUp = (int)((m-1)/stepSize+1);
			int tapDown = (int)((m-1)/stepSize);
			
			// Tap ratios around m:
			double mUp = ratio(tapUp);
			double mDown = ratio(tapDown);
			
			// Cost of each possible tap ratio:
			double costUp = cost(t, mUp);
			double costDown = cost(t, mDown);
			
			// Schedule according to the best option:
			if(costUp < costDown)
			{
				double maxM = ratio(maxStep);
				mSchedule[t] = mUp > maxM ? maxM : mUp;
			}
			else
			{
				double minM = ratio(minStep);
				mSchedule[t] = mDown < minM ? minM : mDown;
			}
		}
		
		// Apply the schedule (this will calculate ratioSchedule allowing for delayed operation):
		applySchedule(mSchedule);

		// Calculate power schedule:
		for(int t = 0; t < T; ++t)
		{
			Complex s = power(t);
			pqSchedule[t] = vector(s.getReal(), s.getImaginary()); // [P Q]
		}
		
		return pqSchedule;
	}

	private Complex power(int t)
	{
		// s = -|v|^2y*, |v| = m (approx.)
		// Note: this is negative because the power is defined as positive when flowing out of the bus.
		Complex s = shunt(ratioSchedule[t]).conjugate().multiply(-ratioSchedule[t]*ratioSchedule[t]);
		return s;
	}

	private Complex shunt(double m)
	{
        return admittance.multiply((1-m)/m/m); // Yj = Yt(1-m)/m^2
	}

	@Override
	public void train(double[][] defaultBandwidth)
	{
		// Nothing to do.
	}

	@Override
	public RealVector uZero()
	{
		return vector(1.0, 0.0); // [P Q]
	}

	@Override
	public int unitCount()
	{
		return B.size();//busCount;
	}

	@Override
	public void setPower0(RealVector power_0)
	{
		this.power_0 = power_0;
	}

	@Override
	public void setPower0(int k, double p)
	{
		this.power_0.setEntry(k, p);
	}

	@Override
	public void setV0(int i, double v)
	{
		this.v0.setEntry(i, v);
	}

	@Override
	public void setDeltaVExternal(RealVector[] deltaVExt)
	{
		this.deltaVExternal = deltaVExt;
	}

	@Override
	public void setDeltaVExternal(int t, int i, double deltaV)
	{
		this.deltaVExternal[t].setEntry(i, deltaV);
	}
	
	@Override
	public RealVector[] getDeltaVExternal()
	{
		return deltaVExternal;
	}

	@Override
	public void setV0(RealVector v0)
	{
		this.v0 = v0;
	}
	
	@Override
	public RealVector getV0()
	{
		return v0;
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
	
	
	//// Accessors ////

//	public RealVector[] getVoltages()
//	{
//		return v;
//	}
//
//	public void setVoltages(RealVector[] v)
//	{
//		this.v = v;
//	}

//	public int getBusCount()
//	{
//		return busCount;
//	}
//
//	public void setBusCount(int busCount)
//	{
//		this.busCount = busCount;
//	}

	public int getT()
	{
		return T;
	}

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

	public int getMinStep()
	{
		return minStep;
	}

	public void setMinStep(int minStep)
	{
		this.minStep = minStep;
	}

	public int getMaxStep()
	{
		return maxStep;
	}

	public void setMaxStep(int maxStep)
	{
		this.maxStep = maxStep;
	}

	public double getStepSize()
	{
		return stepSize;
	}

	public void setStepSize(double stepSize)
	{
		this.stepSize = stepSize;
	}

	public double[] getRatioSchedule()
	{
		return ratioSchedule;
	}

	public void setRatioSchedule(double[] ratio)
	{
		this.ratioSchedule = ratio;
	}
    
    public int getTapChangeDelay()
	{
		return tapChangeDelay;
	}

	public void setTapChangeDelay(int tapChangeDelay)
	{
		this.tapChangeDelay = tapChangeDelay;
	}

	public Complex getAdmittance()
	{
		return admittance;
	}

	public void setAdmittance(Complex admittance)
	{
		this.admittance = admittance;
	}

	public void setSensitivities(RealMatrix sensitivities)
	{
		this.sensistivities = sensitivities;
	}

	public RealVector[] getLoadPower()
	{
		return loadPower;
	}

	public void setLoadPower(RealVector[] loadPower)
	{
		this.loadPower = loadPower;
	}
}
