package ellipsis.test;

import static ellipsis.util.ArrayHelper.maxMagnitude;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.HessianFunction;
import ellipsis.util.NewtonSolver;
import ellipsis.util.VectorHelper;

public class NewtonSolverTest extends TestHelper
{
	public static void main(String[] args)
	{
		assertAssertsOn();
		NewtonSolverTest test = new NewtonSolverTest();
		test.test_univariateQuadratic();
		test.test_multivariateQuadratic();
	}
	
	// Basic quadratic function y=(x1^-1)^2 + (x2-1)^2:
	static HessianFunction f= new HessianFunction()
	{
		@Override
		public RealVector gradient(RealVector x)
		{
			RealVector grad = x.subtract(VectorHelper.vector(x.getDimension(), 1.0)).mapMultiply(2);
			return grad;
		}

		@Override
		public RealMatrix hessian(RealVector x)
		{
			int dimension = x.getDimension();
			RealMatrix I = MatrixUtils.createRealIdentityMatrix(dimension);
			RealMatrix H = I.scalarMultiply(2);
			return H;
		}
	};

	private void test_univariateQuadratic()
	{
		startTest();
		
		// Find minimum:
		NewtonSolver ns = new NewtonSolver();
		ns.setFunction(f);
		ns.setMaxIterations(50);
		ns.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{0.0});
		RealVector x = ns.solve(x0);
//		int k1 = ns.k;
		assertTrue(equals(x.getEntry(0), 1.0, 1e-3), "Minimum was not at 1.0, was at "+x.getEntry(0));
		assertTrue(maxMagnitude(ns.grad.toArray())<ns.getEpsilon(), "Didn't reduce to gradient less than "+ns.getEpsilon()+", last grad was "+maxMagnitude(ns.grad.toArray()));
		
		// Faster with larger epsilon?
//		ns.setEpsilon(1e-3);
//		x = ns.solve(x0);
//		int k2 = ns.k;
//		assertTrue(k2<k1, "Larger epsilon wasn't quicker: k1 = "+k1+", k2="+k2);
//		assertTrue(maxMagnitude(ns.grad.toArray())<ns.getEpsilon(), "Didn't reduce to gradient less than "+ns.getEpsilon()+", last grad was "+maxMagnitude(ns.grad.toArray()));
		
		endTest("univariateQuadratic");
	}
	
	private void test_multivariateQuadratic()
	{
		startTest();
		
		// Find minimum:
		NewtonSolver ns = new NewtonSolver();
		ns.setFunction(f);
		ns.setMaxIterations(50);
		ns.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{0.0, 0.0});
		RealVector x = ns.solve(x0);
//		int k1 = ns.k;
		assertTrue(equals(x.getEntry(0), 1.0, 1e-3), "Minimum was not at 1.0 for element 0, was at "+x.getEntry(0));
		assertTrue(equals(x.getEntry(1), 1.0, 1e-3), "Minimum was not at 1.0 for element 1, was at "+x.getEntry(1));
		assertTrue(maxMagnitude(ns.grad.toArray())<ns.getEpsilon(), "Didn't reduce to gradient less than "+ns.getEpsilon()+", last grad was "+maxMagnitude(ns.grad.toArray()));
		
		// Faster with larger epsilon?
//		ns.setEpsilon(1e-3);
//		x = ns.solve(x0);
//		int k2 = ns.k;
//		assertTrue(k2<k1, "Larger epsilon wasn't quicker: k1 = "+k1+", k2="+k2);
//		assertTrue(maxMagnitude(ns.grad.toArray())<ns.getEpsilon(), "Didn't reduce to gradient less than "+ns.getEpsilon()+", last grad was "+maxMagnitude(ns.grad.toArray()));
		
		endTest("multivariateQuadratic");
	}
}
