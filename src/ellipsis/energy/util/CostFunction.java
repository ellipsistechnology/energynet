package ellipsis.energy.util;

import org.apache.commons.math3.linear.RealVector;

public interface CostFunction
{
	public double g(ADP adp, int t, RealVector x, RealVector u);
	/*{ NOTE I think test will need this - should implement test specific cost function though.
		if(u == null)
			return 0;
		
		// Calculate voltages for state x:
		RealVector w = adp.wZero(t);
		RealVector v = adp.voltagesFromControlAndNoise(t, u, w);
		
		// Fuzzy cost of voltage abs deviations from 1p.u.:
		double cost = 0;
		int size = v.getDimension()/2;
		for(int i = size; i < v.getDimension(); ++i)
		{
			cost += GridlessADP.voltageCost(v.getEntry(i));
		}
		
		return cost/size;
	}*/
}
