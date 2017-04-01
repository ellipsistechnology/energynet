package ellipsis.energy.util.constraint;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;

import ellipsis.energy.util.GridlessADP;
import ellipsis.energy.util.GridlessData;
import ellipsis.energy.util.TimeDependent;
import ellipsis.util.TwiceDifferentiableFunction;

public class VoltageConstraint implements TwiceDifferentiableFunction, TimeDependent, ConvexMultivariateRealFunction
{
	// Parameters (read-only):
	private GridlessADP adp;
	private double a; // a|v_i| + b <= 0
	private double b;
	private GridlessData gridlessData;
	private int dim;

	int i; // index of bus
	
	// Caches:
	private int unitCount;
	private RealVector[] wZeros;
	
	// Caches:
	private Map<RealVector, Double> valueCache;
	
	// Variables:
	private int t;

	public VoltageConstraint(GridlessADP adp, int dimension)
	{
		this.adp = adp;
		this.dim = dimension;
		gridlessData = adp.getGridlessData();
		
		// Fill in caches:
		unitCount = adp.unitCount();
		
		wZeros = new RealVector[adp.getT()+1];
		for(int t = 0; t < wZeros.length; ++t)
			wZeros[t] = adp.wZero(t);
	}
	
//public static int count_g = 0;
	/**
	 * gradient = a*d|v|/du.
	 */
	@Override
	public RealVector gradient(RealVector u)
	{
//++count_g;
//Timer.getGlobalTimer("VC.g").start();
		int loadCount = gridlessData.loadCount;
		int dgCount = gridlessData.dgCount;
		int unitCount = dgCount+loadCount;

		int dimension = 2*dgCount;
		double[] dvdu = new double[dimension];
		double[][] sensitivitiesArray = gridlessData.sensitivitiesArray();
		for(int j = 0; j < dgCount; ++j)
		{
			dvdu[j] = sensitivitiesArray[i+unitCount][j]*a; // I think this is transposed (thus i & j swapped).
			dvdu[j+dgCount] = sensitivitiesArray[i+unitCount][j+dgCount]*a;
		}

//Timer.getGlobalTimer("VC.g").stop();
		return new ArrayRealVector(dvdu, false);
	}
	
	@Override
	public RealMatrix hessian(RealVector u)
	{
		return null; // null means [0]
	}

//public static int count_v = 0;
	@Override
	public double value(RealVector u)
	{
//++count_v;
//Timer.getGlobalTimer("VC.v").start();
		if(valueCache.containsKey(u))
		{
//Timer.getGlobalTimer("VC.v").stop();
			return valueCache.get(u);
		}
		
		int unitCount = gridlessData.dgCount+gridlessData.loadCount;

		RealVector w = adp.wZero(t);
		RealVector power = adp.powerFromControlAndNoise(u, w);
		RealVector deltaS = power.subtract(gridlessData.power_0);

		// Calculate |v_i|:
		// FIXME make this more efficient:
		RealVector sensitivities_i = new ArrayRealVector(gridlessData.sensitivities.getRow(i + unitCount));
		double v_abs_i = gridlessData.v_0.getEntry(unitCount+i) + sensitivities_i.dotProduct(deltaS);
		
		double v = a*v_abs_i + b;
		valueCache.put(u, v);
//Timer.getGlobalTimer("VC.v").stop();
		return v;
	}

	@Override
	public void setTime(int t)
	{
		valueCache = new HashMap<>();
		this.t = t;
	}
	
	@Override
	public String toString()
	{
		return a+"|v["+i+"]| + "+b+" <= 0";
	}
	
	
	//// JOptimizer ////

	@Override
	public double value(double[] u)
	{
		RealVector power = adp.powerFromControlAndNoise(new ArrayRealVector(u), wZeros[t]);
		RealVector deltaS = power.subtract(gridlessData.power_0);

		// Calculate |v_i|:
		RealVector sensitivities_i = new ArrayRealVector(gridlessData.sensitivities.getRow(unitCount+i));
		double v_0_i = gridlessData.v_0.getEntry(unitCount+i);
		double deltaV_i = sensitivities_i.dotProduct(deltaS);
		double deltaVExt_i = gridlessData.deltaV_external[t].getEntry(unitCount+i);
		double v_abs_i = v_0_i + deltaV_i + deltaVExt_i;
		
		double v = a*v_abs_i + b;

		return v;
	}

	@Override
	public double[] gradient(double[] X)
	{
		int dgCount = gridlessData.dgCount;

		int dimension = 2*dgCount;
		double[] dvdu = new double[dimension];
		double[][] sensitivitiesArray = gridlessData.sensitivitiesArray();
		for(int j = 0; j < dgCount; ++j)
		{
			dvdu[j] = sensitivitiesArray[i+unitCount][j]*a; // d|v|/dP
			dvdu[j+dgCount] = sensitivitiesArray[i+unitCount][j+unitCount]*a; // d|v|/dQ
		}

		return dvdu;
	}

	@Override
	public double[][] hessian(double[] X)
	{
		return new double[dim][dim];
	}

	@Override
	public int getDim()
	{
		return dim;
	}
	
	
	//// Accessors ////

	public double getA()
	{
		return a;
	}

	public void setA(double a)
	{
		this.a = a;
	}

	public double getB()
	{
		return b;
	}

	public void setB(double b)
	{
		this.b = b;
	}

	public int getI()
	{
		return i;
	}

	public void setI(int i)
	{
		this.i = i;
	}
}