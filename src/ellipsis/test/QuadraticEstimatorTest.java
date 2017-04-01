package ellipsis.test;

import java.util.Random;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.LinearLeastSquaresEstimator;
import ellipsis.util.QuadraticEstimator;
import ellipsis.util.VectorHelper;

public class QuadraticEstimatorTest extends LinearLeastSquaresTest
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		QuadraticEstimatorTest test = new QuadraticEstimatorTest();
		test.parseArgs(args);

		test.test_univariateGeneralQuadraticFunction();
		test.test_univariateQuadraticFunction();
		test.test_univariateQuadraticFunctionWithNoisySamples();
		test.test_univariateGeneralQuadraticFunctionWithNoisySamples();
		test.test_multivariateQuadraticFunctionWithNoisySamples();
		test.test_multivariateGeneralQuadraticFunctionWithNoisySamples();
	}
	
	@Override
	protected void test_univariateQuadraticFunction(double a, double b, double xShift, double noise, double tolerance, Random rand)
	{
		// Generate samples:
		LinearLeastSquaresEstimator estimator = new LinearLeastSquaresEstimator();
		QuadraticEstimator quad = new QuadraticEstimator();
		for(int i = 0; i < 30; ++i)
		{
			double x = 1-2*rand.nextDouble();
			double _x = x + xShift;
			double y = a*_x*_x + b + noise*(1-2*rand.nextGaussian());
			RealVector x2 = VectorHelper.vector(_x*_x);
			estimator.addSample(x2, y);
			quad.addSample(VectorHelper.vector(x), y);
		}

		// Test:
		for(double x = -1; x < 1; x += 0.1)
		{
			double _x = x+xShift;
			RealVector x2 = VectorHelper.vector(_x*_x);
			double estimate = estimator.value(x2);
			double quadEstimate = quad.value(VectorHelper.vector(x));
			assertTrue(equals(estimate, quadEstimate, tolerance), "LLS estimate did not match QuadraticEstimator estimate;  with a="+a+", b="+b+", at x="+x+"; LLS estimate ="+estimate+", quad estimate="+quadEstimate);
		}
	}

	@Override
	protected void test_multivariateQuadraticFunction(RealVector a, double b, double noise, double tolerance, Random rand)
	{
		// Generate samples:
		LinearLeastSquaresEstimator estimator = new LinearLeastSquaresEstimator();
		QuadraticEstimator quad = new QuadraticEstimator();
		for(int i = 0; i < 60; ++i)
		{
			RealVector x = VectorHelper.vector(1-2*rand.nextDouble(), 1-2*rand.nextDouble());
			double y = a.ebeMultiply(x).dotProduct(x) + b + noise*(1-2*rand.nextGaussian());
			RealVector x2 = x.ebeMultiply(x);
			estimator.addSample(x2, y);
			quad.addSample(x, y);
		}

		// Test:
		for(double x1 = -1; x1 < 1; x1 += 0.1)
		{
			for(double x2 = -1; x2 < 1; x2 += 0.1)
			{
				RealVector x = VectorHelper.vector(x1, x2);
				RealVector xsquared = x.ebeMultiply(x);
				double estimate = estimator.value(xsquared);
				double quadEstimate = quad.value(x);
				assertTrue(equals(estimate, quadEstimate, tolerance), "LLS estimate did not match QuadraticEstimator estimate;  with a="+a+", b="+b+", at x="+x+"; LLS estimate ="+estimate+", quad estimate="+quadEstimate);
			}
		}
	}
}