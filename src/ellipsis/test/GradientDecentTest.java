package ellipsis.test;

import static ellipsis.util.ArrayHelper.maxMagnitude;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.GradientDecent;
import ellipsis.util.GradientFunction;

public class GradientDecentTest extends TestHelper
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		GradientDecentTest test = new GradientDecentTest();
		test.test_univariateQuadratic();
		test.test_multivariateQuadratic();
	}
	
	// Basic quadratic function:
	static GradientFunction f= new GradientFunction()
	{
		@Override
		public RealVector gradient(RealVector x)
		{
			return x.mapMultiply(2);
		}
	};

	private void test_univariateQuadratic()
	{
		startTest();
		
		// Find minimum:
		GradientDecent gd = new GradientDecent();
		gd.setFunction(f);
		gd.setStep(0.2);
		gd.setMaxIterations(50);
		gd.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{1.0});
		RealVector x = gd.solve(x0);
		int k1 = gd.k;
		assertTrue(equals(x.getEntry(0), 0.0, 1e-3), "Minimum was not at zero, was at "+x.getEntry(0));
		assertTrue(maxMagnitude(gd.grad.toArray())<gd.getEpsilon(), "Didn't reduce to gradient less than "+gd.getEpsilon()+", last grad was "+maxMagnitude(gd.grad.toArray()));
		
		// Faster with larger epsilon?
		gd.setEpsilon(1e-3);
		x = gd.solve(x0);
		int k2 = gd.k;
		assertTrue(k2<k1, "Larger epsilon wasn't quicker: k1 = "+k1+", k2="+k2);
		assertTrue(maxMagnitude(gd.grad.toArray())<gd.getEpsilon(), "Didn't reduce to gradient less than "+gd.getEpsilon()+", last grad was "+maxMagnitude(gd.grad.toArray()));
		
		endTest("univariateQuadratic");
	}
	
	private void test_multivariateQuadratic()
	{
		startTest();
		
		// Find minimum:
		GradientDecent gd = new GradientDecent();
		gd.setFunction(f);
		gd.setStep(0.2);
		gd.setMaxIterations(50);
		gd.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{1.0, 1.0});
		RealVector x = gd.solve(x0);
		int k1 = gd.k;
		assertTrue(equals(x.getEntry(0), 0.0, 1e-3), "Minimum was not at zero for element 0, was at "+x.getEntry(0));
		assertTrue(equals(x.getEntry(1), 0.0, 1e-3), "Minimum was not at zero for element 1, was at "+x.getEntry(1));
		assertTrue(maxMagnitude(gd.grad.toArray())<gd.getEpsilon(), "Didn't reduce to gradient less than "+gd.getEpsilon()+", last grad was "+maxMagnitude(gd.grad.toArray()));
		
		// Faster with larger epsilon?
		gd.setEpsilon(1e-3);
		x = gd.solve(x0);
		int k2 = gd.k;
		assertTrue(k2<k1, "Larger epsilon wasn't quicker: k1 = "+k1+", k2="+k2);
		assertTrue(maxMagnitude(gd.grad.toArray())<gd.getEpsilon(), "Didn't reduce to gradient less than "+gd.getEpsilon()+", last grad was "+maxMagnitude(gd.grad.toArray()));
		
		endTest("multivariateQuadratic");
	}
}