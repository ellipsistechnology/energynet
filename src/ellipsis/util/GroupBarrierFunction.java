package ellipsis.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class GroupBarrierFunction implements BarrierFunction
{
	private List<BarrierFunction> barriers = new ArrayList<>();
	
	@Override
	public RealMatrix hessian(RealVector x)
	{
		int dimension = x.getDimension();
		RealMatrix m = new Array2DRowRealMatrix(dimension, dimension);
		for (TwiceDifferentiableFunction b : barriers)
		{
			m = m.add(b.hessian(x));
		}
		return m;
	}

	@Override
	public RealVector gradient(RealVector x)
	{
		RealVector g = new ArrayRealVector(x.getDimension());
		for (TwiceDifferentiableFunction b : barriers)
		{
			g = g.add(b.gradient(x));
		}
		return g;
	}

	@Override
	public double value(RealVector x)
	{
		double v = 0;
		for (TwiceDifferentiableFunction b : barriers)
		{
			v += b.value(x);
		}
		return v;
	}

	public boolean addBarrier(BarrierFunction e)
	{
		return barriers.add(e);
	}
	
	public List<BarrierFunction> getBarriers()
	{
		return barriers;
	}

	@Override
	public List<TwiceDifferentiableFunction> getConstraints()
	{
		List<TwiceDifferentiableFunction> constraints = new ArrayList<>();
		for (BarrierFunction barrier : barriers)
		{
			constraints.addAll(barrier.getConstraints());
		}
		return constraints;
	}
}