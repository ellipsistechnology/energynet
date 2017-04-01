package ellipsis.test;

import java.util.Random;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.LinearLeastSquaresEstimator;
import ellipsis.util.VectorHelper;

public class LinearLeastSquaresTest extends TestHelper
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		LinearLeastSquaresTest test = new LinearLeastSquaresTest();
		test.parseArgs(args);
		
		// Linear:
		test.test_univariateLinearFunction();
		test.test_univariateLinearFunctionWithNoisySamples();
		test.test_univariateGeneralLinearFunction();
		
		// Quadratic:
		test.test_univariateGeneralQuadraticFunction();
		test.test_univariateQuadraticFunction();
		test.test_univariateQuadraticFunctionWithNoisySamples();
		test.test_univariateGeneralQuadraticFunctionWithNoisySamples();
		test.test_multivariateQuadraticFunctionWithNoisySamples();
		test.test_multivariateGeneralQuadraticFunctionWithNoisySamples();
	}
	
	
	//// Linear ////
	
	private void test_univariateGeneralLinearFunction()
	{
		Random rand = new Random(0);
		startTest();
		
		// y = ax
		test_univariateLinearFunction(1, 1, 0, 1e-12, rand);
		test_univariateLinearFunction(-1, -1, 0, 1e-12, rand);
		test_univariateLinearFunction(0, 1, 0, 1e-12, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateLinearFunction(1-2*rand.nextDouble(), 1-2*rand.nextDouble(), 0, 1e-12, rand); // [-1, +1]
		
		endTest("test_univariateGeneralLinearFunction");
	}

	protected void test_univariateLinearFunction()
	{
		Random rand = new Random(0);
		startTest();
		
		// y = ax
		test_univariateLinearFunction(1, 0, 0, 1e-12, rand);
		test_univariateLinearFunction(-1, 0, 0, 1e-12, rand);
		test_univariateLinearFunction(0, 0, 0, 1e-12, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateLinearFunction(1-2*rand.nextDouble(), 0, 0, 1e-12, rand); // [-1, +1]
		
		endTest("test_univariateLinearFunction");
	}
	
	protected void test_univariateLinearFunctionWithNoisySamples()
	{
		Random rand = new Random(0);
		startTest();
		
		// y = ax
		test_univariateLinearFunction(1, 0, 0.025, 0.1, rand);
		test_univariateLinearFunction(-1, 0, 0.025, 0.1, rand);
		test_univariateLinearFunction(0, 0, 0.025, 0.1, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateLinearFunction(1-2*rand.nextDouble(), 0, 0.025, 0.12, rand); // [-1, +1]
		
		endTest("test_univariateLinearFunctionWithNoisySamples");
	}

	/**
	 * Test y = ax + b.
	 * @param a
	 * @param b
	 * @param noise
	 * @param tolerance
	 * @param rand
	 */
	protected void test_univariateLinearFunction(double a, double b, double noise, double tolerance, Random rand)
	{
		// Generate samples:
		LinearLeastSquaresEstimator estimator = new LinearLeastSquaresEstimator();
		for(int i = 0; i < 10; ++i)
		{
			double x = 1-2*rand.nextDouble();
			double y = a*x + b + noise*(1-2*rand.nextGaussian());
			estimator.addSample(VectorHelper.vector(x), y);
		}

		// Test:
		for(double x = -1; x < 1; x += 0.1)
		{
			double estimate = estimator.value(VectorHelper.vector(x));
			double y = a*x + b;
			assertTrue(equals(estimate, y, tolerance), "Estimate did not match with a="+a+", b="+b+", at x="+x+"; was "+estimate+", should be y="+y);
		}
	}
	
	
	//// Quadratic ////
	
	protected void test_univariateQuadraticFunction()
	{
		Random rand = new Random(0);
		startTest();

		// y = ax^2
		test_univariateQuadraticFunction(1, 0, 0, 0, 1e-12, rand);
		test_univariateQuadraticFunction(-1, 0, 0, 0, 1e-12, rand);
		test_univariateQuadraticFunction(0, 0, 0, 0, 1e-12, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateQuadraticFunction(1-2*rand.nextDouble(), 0, 0, 0, 1e-12, rand); // [-1, +1]

		endTest("test_univariateQuadraticFunction");
	}
	
	protected void test_univariateGeneralQuadraticFunction()
	{
		Random rand = new Random(0);
		startTest();

		// y = ax^2
		test_univariateQuadraticFunction(1, 1, 0.5, 0, 1e-12, rand);
		test_univariateQuadraticFunction(-1, -1, 0.25, 0, 1e-12, rand);
		test_univariateQuadraticFunction(0, 1, 0, 0, 1e-12, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateQuadraticFunction(1-2*rand.nextDouble(), 1-2*rand.nextDouble(), 1-2*rand.nextDouble(), 0, 1e-12, rand); // [-1, +1]

		endTest("test_univariateQuadraticFunction");
	}
	
	protected void test_univariateQuadraticFunctionWithNoisySamples()
	{
		Random rand = new Random(0);
		startTest();

		// y = ax^2
		test_univariateQuadraticFunction(1, 0, 0, 0.025, 0.1, rand);
		test_univariateQuadraticFunction(-1, 0, 0, 0.025, 0.1, rand);
		test_univariateQuadraticFunction(0, 0, 0, 0.025, 0.1, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateQuadraticFunction(1-2*rand.nextDouble(), 0, 0, 0.025, 0.15, rand); // [-1, +1]

		endTest("test_univariateQuadraticFunctionWithNoisySamples");
	}
	
	protected void test_univariateGeneralQuadraticFunctionWithNoisySamples()
	{
		Random rand = new Random(0);
		startTest();

		// y = ax^2
		test_univariateQuadraticFunction(1, 1, 0.5, 0.025, 0.1, rand);
		test_univariateQuadraticFunction(-1, -1, 0.25, 0.025, 0.1, rand);
		test_univariateQuadraticFunction(0, 1, 0, 0.025, 0.1, rand);
		for(int i = 0; i < 10; ++i)
			test_univariateQuadraticFunction(1-2*rand.nextDouble(), 1-2*rand.nextDouble(), 1-2*rand.nextDouble(), 0.025, 0.15, rand); // [-1, +1]

		endTest("test_univariateQuadraticFunctionWithNoisySamples");
	}

	// Test a*(x+xShift)^2 + b
	protected void test_univariateQuadraticFunction(double a, double b, double xShift, double noise, double tolerance, Random rand)
	{
		// Generate samples:
		LinearLeastSquaresEstimator estimator = new LinearLeastSquaresEstimator();
		for(int i = 0; i < 10; ++i)
		{
			double x = 1-2*rand.nextDouble();
			double _x = (x+xShift);
			double x2 = _x*_x;
			double y = a * x2 + b + noise*(1-2*rand.nextGaussian());
			estimator.addSample(VectorHelper.vector(x2), y);
		}

		// Test:
		for(double x = -1; x < 1; x += 0.1)
		{
			double _x = (x+xShift);
			double x2 = _x*_x;
			double estimate = estimator.value(VectorHelper.vector(x2));
			double y = a*x2 + b;
			assertTrue(equals(estimate, y, tolerance), "Estimate did not match with a="+a+", b="+b+", xShift="+xShift+", at x="+x+"; was "+estimate+", should be y="+y);
		}
	}
	
	
	//// Multivariate Quadratic ////

	protected void test_multivariateQuadraticFunctionWithNoisySamples()
	{
		Random rand = new Random(0);
		startTest();

		// y = x^TAx, for diagonal matrix A
		test_multivariateQuadraticFunction(VectorHelper.vector(1.0, 1.0), 0, 0.025, 0.15, rand);
		test_multivariateQuadraticFunction(VectorHelper.vector(-1.0, -1.0), 0, 0.025, 0.15, rand);
		test_multivariateQuadraticFunction(VectorHelper.vector(0.0, 0.0), 0, 0.025, 0.15, rand);
		for(int i = 0; i < 10; ++i)
		{
			RealVector a = VectorHelper.vector(1-2*rand.nextDouble(), 1-2*rand.nextDouble()); // 2D
			test_multivariateQuadraticFunction(a, 0, 0.025, 0.2, rand); // [-1, +1]
		}

		endTest("test_multivariateQuadraticFunctionWithNoisySamples");
	}

	protected void test_multivariateGeneralQuadraticFunctionWithNoisySamples()
	{
		Random rand = new Random(0);
		startTest();

		// y = x^TAx, for diagonal matrix A
		test_multivariateQuadraticFunction(VectorHelper.vector(1.0, 1.0), 1, 0.025, 0.15, rand);
		test_multivariateQuadraticFunction(VectorHelper.vector(-1.0, -1.0), -1, 0.025, 0.15, rand);
		test_multivariateQuadraticFunction(VectorHelper.vector(0.0, 0.0), 1, 0.025, 0.15, rand);
		for(int i = 0; i < 10; ++i)
		{
			RealVector a = VectorHelper.vector(1-2*rand.nextDouble(), 1-2*rand.nextDouble()); // 2D
			test_multivariateQuadraticFunction(a, 1-2*rand.nextDouble(), 0.025, 0.2, rand); // [-1, +1]
		}

		endTest("test_multivariateQuadraticFunctionWithNoisySamples");
	}

	// Test y = x^Tax+b
	protected void test_multivariateQuadraticFunction(RealVector a, double b, double noise, double tolerance, Random rand)
	{
		// Generate samples:
		LinearLeastSquaresEstimator estimator = new LinearLeastSquaresEstimator();
		for(int i = 0; i < 10; ++i)
		{
			RealVector x = VectorHelper.vector(1-2*rand.nextDouble(), 1-2*rand.nextDouble());
			double y = a.ebeMultiply(x).dotProduct(x) + b + noise*(1-2*rand.nextGaussian());
			estimator.addSample(x.ebeMultiply(x), y);
		}

		// Test:
		for(double x1 = -1; x1 < 1; x1 += 0.1)
		{
			for(double x2 = -1; x2 < 1; x2 += 0.1)
			{
				RealVector x = VectorHelper.vector(x1, x2);
				double estimate = estimator.value(x.ebeMultiply(x));
				double y = a.ebeMultiply(x).dotProduct(x) + b;
				assertTrue(equals(estimate, y, tolerance), "Estimate did not match with a="+a+", b="+b+", at x=["+x1+", "+x2+"]; was "+estimate+", should be y="+y);
			}
		}
	}
}