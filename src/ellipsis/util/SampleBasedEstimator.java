package ellipsis.util;

import java.util.List;

import org.apache.commons.math3.linear.RealVector;

public interface SampleBasedEstimator
{
	void addSample(RealVector x, double y);
	double value(RealVector x);
	void setDefaultValue(double def);
	List<Pair<RealVector, Double>> getSamplePoints();
}