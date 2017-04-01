package ellipsis.test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.StrictlyConvexMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import com.joptimizer.optimizers.OptimizationResponse;

import ellipsis.util.BarrierSolver;
import ellipsis.util.HessianFunction;
import ellipsis.util.LogBarrierFunction;
import ellipsis.util.TwiceDifferentiableFunction;
import ellipsis.util.VectorHelper;

public class BarrierSolverTest extends TestHelper
{
	private static final double CONSTRAINT = 1.2;

	public static void main(String[] args)
	{
		assertAssertsOn();
		BarrierSolverTest test = new BarrierSolverTest();
		test.test_quadratic();
		test.test_logBarrierFunction();
		test.test_JOptimizer();
//		test.test_LAML();
	}
	
//	private void test_LAML()
//	{
//		startTest();
//		
//		// Objective (f(x) = 0.5xQx + cx):
//		Matrix Q = new DenseMatrix(new double[][]{
//				{1,0},	
//				{0,1}
//		});
//		Matrix c = new DenseMatrix(new double[][]{
//				{0},	
//				{0}
//		});
//		
//		// Equality constraints (Ax = b):
//		Matrix A = new DenseMatrix(new double[][]{
//				{0,0},	
//				{0,0}
//		});
//		Matrix b = new DenseMatrix(new double[][]{
//				{0},	
//				{0}
//		});
//		
//		// Inequality constraints (Bx < d: 0.2 < x => -x < -0.2):
//		Matrix B = new DenseMatrix(new double[][]{
//				{-1,  0},	
//				{ 0, -1}
//		});
//		Matrix d = new DenseMatrix(new double[][]{
//				{-0.2},	
//				{-0.2}
//		});
//		
//		// Initial x, lambda, nu:
//		Matrix x = new DenseMatrix(new double[][]{
//				{0.5},	
//				{0.5}
//		});
//		
//		test_LAML_basicQuadratic(Q, c, A, b, B, d, x);
//		test_LAML_complexQuadratic(Q, c, A, b, B, d, x);
//		test_LAML_quadraticConstraints();
//		
//		endTest("LAML - PrimalDualInteriorPoint");
//	}
//
//	private void test_LAML_quadraticConstraints()
//	{
//		Matrix Q = new DenseMatrix(new double[][]{
//				{0.00687, 0},
//				{0, 0.00687}
//		});
//		Matrix c = new DenseMatrix(new double[][]{
//				{-0.00687*17.6},
//				{0}
//		});
//		double l = 0.00687*78.9;
//		Objective f = new QuadraticObjective(Q, c, l);
//		
//		
//		for(double x1 = 0; x1 < 10; x1 += 0.1)
//			System.out.println(x1+","+f.value(new DenseMatrix(new double[][]{
//					{x1},
//					{-1}
//			})));
//
//
////		final DeviceConstraint dc = new DeviceConstraint(0, new double[]{0}, null);
////		final ApparentPowerConstraint apc = new ApparentPowerConstraint(0, new double[]{0.01});
////		dc.setTime(0);
////		apc.setTime(0);
//		Constraints F = new Constraints()
//		{
//			@Override
//			public Matrix value(Matrix x)
//			{
////				RealVector v = toVector(x);
//				double x1 = x.getEntry(0, 0);
//				double x2 = x.getEntry(1, 0);
//				return new DenseMatrix(new double[][]{
//						{-x1}, // x1 > 0
//						{-x2}//{x1+x2-1000}//{x1*x1+x2*x2-1000.01}
//				});
//			}
//			
//			@Override
//			public Matrix derivative(Matrix x)
//			{
////				RealVector v = toVector(x);
//				double x1 = x.getEntry(0, 0);
//				double x2 = x.getEntry(1, 0);
//				return new DenseMatrix(new double[][]{
//						{-1,0},
//						{0,-1}
//				});
//			}
//			
//			@Override
//			public Matrix hessian(int i, Matrix x)
//			{
////				RealVector v = toVector(x);
//				if(i == 0)
//				{
////					RealMatrix hessian = dc.hessian(v);
////					if(hessian == null)
////						return new DenseMatrix(new double[][]{
////								{0, 0},
////								{0, 0}
////						});
////					return toMatrix(hessian);
//					return new DenseMatrix(new double[][]{
//							{0,0},
//							{0,0}
//					});
//				}
//				else
//				{
////					RealMatrix hessian = apc.hessian(v);
////					if(hessian == null)
////						return new DenseMatrix(new double[][]{
////								{0, 0},
////								{0, 0}
////						});
////					return toMatrix(hessian);
//					return new DenseMatrix(new double[][]{
//							{0,0},//{2,0},
//							{0,0}//{0,2}
//					});
//				}
//			}
//			
//			@Override
//			public int constrantCount()
//			{
//				return 2;
//			}
//		};
//		Matrix A = new DenseMatrix(new double[][]{
//				{0,0}
//		});
//		Matrix b = new DenseMatrix(new double[][]{
//				{0}
//		});
//		PrimalDualInteriorPointSolver solver = new PrimalDualInteriorPointSolver(f, F , A, b);
//		Matrix x = new DenseMatrix(new double[][]{
//				{1e-12},
//				{1e-12}
//		});
//		
//		
////		double fval = Matlab.innerProduct(x, Q.mtimes(x))/2 + Matlab.innerProduct(c, x); 
////		Matrix H_x = Q; 
////		DenseMatrix B = new DenseMatrix(new double[][]{
////				{-1,0},
////				{0,-1}
////		});
////		DenseMatrix d = new DenseMatrix(new double[][]{
////				{0},
////				{0}
////		});
////		Matrix DF_x = B;
////		Matrix F_x = B.mtimes(x).minus(d); 
////		Matrix G_f_x = Q.mtimes(x).plus(c); 
////		boolean flags[] = null; 
////		int k = 0; 
////		while (true) { 
////		  flags = PrimalDualInteriorPoint.run(A, b, H_x, F_x, DF_x, G_f_x, fval, x); 
////		  if (flags[0]) 
////		    break; 
////		  // Compute the objective function value, if flags[1] is true 
////		  // gradient will also be computed. 
////		  fval = Matlab.innerProduct(x, Q.mtimes(x)) + Matlab.innerProduct(c, x); 
////		  F_x = B.mtimes(x).minus(d); 
////		  if (flags[1]) { 
////		    k = k + 1; 
////		    // Compute the gradient 
////		    G_f_x = Q.mtimes(x).times(2).plus(c); 
////		  } 
////		}
////		System.out.println(PrimalDualInteriorPoint.x);
//		
//		solver.solve(x);
//		System.out.println(solver.getX());
//	}
//	
//	private RealVector toVector(Matrix x)
//	{
//		int dimension = x.getRowDimension();
//		double[][] xArray = x.getData();
//		double[] v = new double[dimension];
//		for(int i = 0; i < dimension; ++i)
//		{
//			v[i] = xArray[i][0];
//		}
//		return new ArrayRealVector(v);
//	}
//	
//	private Matrix toMatrix(RealMatrix rm)
//	{
//		double[][] m = new double[rm.getRowDimension()][rm.getColumnDimension()];
//		for(int i = 0; i < rm.getRowDimension(); ++i)
//		{
//			for(int j = 0; j < rm.getColumnDimension(); ++j)
//			{
//				m[i][j] = rm.getEntry(i, j);
//			}
//		}
//		return new DenseMatrix(m);
//	}
//
//	private void test_LAML_complexQuadratic(
//			Matrix Q, 
//			Matrix c, 
//			Matrix A,
//			Matrix b, 
//			Matrix B, 
//			Matrix d, 
//			Matrix x0)
//	{
//		// Objective is sum{(1-x)^2} = xIx -2x + 2
//		c = new DenseMatrix(new double[][]{
//				{-2},	
//				{-2}
//		});
//		
//		PrimalDualInteriorPointSolver.Objective f = new QuadraticObjective(Q, c, 2);
//		PrimalDualInteriorPointSolver.Constraints F = new LinearConstraints(B, d);
//		PrimalDualInteriorPointSolver solver = new PrimalDualInteriorPointSolver(f, F, A, b);
//		Matrix x = x0.copy();
//		boolean converged = solver.solve(x);
//		x = solver.getX();
//		
////System.out.println();
////for(double x1 = -1; x1 <= 1; x1 += 0.1)
////{
////	for(double x2 = -1; x2 <= 1; x2 += 0.1)
////	{
////		System.out.print(f.value(new DenseMatrix(new double[][]{{x1},{x2}})));
////		System.out.print(',');
////	}
////	System.out.println();
////}
////System.out.println();
//		
//		assertTrue(converged, "Did not converge");
//		assertTrue(equals(x.getEntry(0, 0), 1, 1e-4), "Found incorrect solution for f(x)=xIx-2x+2, 0.2-x<0, was "+x+", should be [1, 1]");
//		assertTrue(equals(x.getEntry(1, 0), 1, 1e-4), "Found incorrect solution for f(x)=xIx-2x+2, 0.2-x<0, was "+x+", should be [1, 1]");
//		
//
//
//		// Objective is 0.00708x^2 - 0.130x + 0.611:
//		Q = new DenseMatrix(new double[][]{
//				{0.00708}	
//		});
//		c = new DenseMatrix(new double[][]{
//				{-0.13},
//		});
//		
//		// Constraint 0 < x:
//		B = new DenseMatrix(new double[][]{
//				{-1}	
//		});
//		d = new DenseMatrix(new double[][]{
//				{0}	
//		});
//		
//		// Ax = b:
//		A = new DenseMatrix(new double[][]{
//				{0}	
//		});
//		b = new DenseMatrix(new double[][]{
//				{0}	
//		});
//		
//		f = new QuadraticObjective(Q, c, 0.611);
//		F = new LinearConstraints(B, d);
//		solver = new PrimalDualInteriorPointSolver(f, F, A, b);
//		x = new DenseMatrix(new double[][]{
//				{1e-9}	
//		});
//		converged = solver.solve(x);
//		x = solver.getX();
//		
//		assertTrue(converged, "Did not converge");
//		assertTrue(equals(x.getEntry(0, 0), 9.18079, 1e-5), "Found incorrect solution for f(x)=0.00708x^2 - 0.130x + 0.611, -x < 0, was "+x+", should be 9.18079");
//	}
//
//	private void test_LAML_basicQuadratic(
//			final Matrix Q, 
//			final Matrix c, 
//			final Matrix A, 
//			final Matrix b, 
//			final Matrix B, 
//			final Matrix d, 
//			final Matrix x0)
//	{
//		// The following approach is no longer supported:
//		// Initial state:
////		Matrix x = x0.copy();
////		Matrix Qx = Q.mtimes(x);
////		double cx = Matlab.innerProduct(c, x);
////		double fval = Matlab.innerProduct(x, Qx)/2 + cx; 
////		Matrix H_x = Q; 
////		Matrix DF_x = B; 
////		Matrix F_x = B.mtimes(x).minus(d);
////		Matrix G_f_x = Qx.plus(c);
//		
//		// Solve with manual loop:
////		boolean flags[] = null;
////		PrimalDualInteriorPointSolver solver = new PrimalDualInteriorPointSolver();
////		solver.init(x);
////		while (true)
////		{
//////System.out.println("1. x="+x);
////			flags = solver.step(A, b, H_x, F_x, DF_x, G_f_x, fval, x);
////			if (flags[0])
////				break;
////			// Compute the objective function value, if flags[1] is true
////			// gradient will also be computed.
////			fval = Matlab.innerProduct(x, Q.mtimes(x)) / 2 + cx;
////			F_x = B.mtimes(x).minus(d);
////			if (flags[1])
////			{
////				// Compute the gradient
////				G_f_x = Qx.plus(c);
////			}
////		}
////		x = solver.getX();
////		
////		assertTrue(equals(x.getEntry(0, 0), 0.2, 1e-6), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
////		assertTrue(equals(x.getEntry(1, 0), 0.2, 1e-6), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
//		
//		// Solve with solver's loop:
//		PrimalDualInteriorPointSolver.Constraints F1 = new PrimalDualInteriorPointSolver.Constraints()
//		{
//			@Override
//			public Matrix derivative(Matrix x)
//			{
//				return B;
//			}
//
//			@Override
//			public Matrix value(Matrix x)
//			{
//				return B.mtimes(x).minus(d);
//			}
//
//			@Override
//			public Matrix hessian(int i, Matrix x)
//			{
//				return new DenseMatrix(x.getRowDimension(), x.getRowDimension()); // zero
//			}
//
//			@Override
//			public int constrantCount()
//			{
//				return 1;
//			}
//		};
//		PrimalDualInteriorPointSolver.Objective f1 = new PrimalDualInteriorPointSolver.Objective()
//		{
//			@Override
//			public double value(Matrix x)
//			{
//				Matrix Qx = Q.mtimes(x);
//				double cx = Matlab.innerProduct(c, x);
//				double fval = Matlab.innerProduct(x, Qx)/2 + cx;
//				return fval;
//			}
//
//			@Override
//			public Matrix hessian(Matrix x)
//			{
//				return Q;
//			}
//
//			@Override
//			public Matrix gradient(Matrix x)
//			{
//				return Q.mtimes(x).plus(c);
//			}
//		};
//		PrimalDualInteriorPointSolver solver = new PrimalDualInteriorPointSolver(f1, F1, A, b);
//		Matrix x = x0.copy();
//		boolean converged = solver.solve(x);
//		x = solver.getX();
//		
//		assertTrue(converged, "Did not converge");
//		assertTrue(equals(x.getEntry(0, 0), 0.2, 1e-4), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
//		assertTrue(equals(x.getEntry(1, 0), 0.2, 1e-4), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
//		
//		// Use QuadraticObjective and LinearConstraints:
//		QuadraticObjective f2 = new QuadraticObjective(Q, c, 0);
//		LinearConstraints F2 = new LinearConstraints(B, d);
//		
////System.out.println();
////for(double x1 = -1; x1 <= 1; x1 += 0.1)
////{
////	for(double x2 = -1; x2 <= 1; x2 += 0.1)
////	{
////		System.out.print(f1.value(new DenseMatrix(new double[][]{{x1},{x2}})));
////		System.out.print(',');
////	}
////	System.out.println();
////}
////System.out.println();
////System.out.println();
////for(double x1 = -1; x1 <= 1; x1 += 0.1)
////{
////	for(double x2 = -1; x2 <= 1; x2 += 0.1)
////	{
////		System.out.print(f2.value(new DenseMatrix(new double[][]{{x1},{x2}})));
////		System.out.print(',');
////	}
////	System.out.println();
////}
////System.out.println();
//		solver = new PrimalDualInteriorPointSolver(f2, F2, A, b);
//		x = x0.copy();
//		converged = solver.solve(x);
//		x = solver.getX();
//		
//		assertTrue(converged, "Did not converge");
//		assertTrue(equals(x.getEntry(0, 0), 0.2, 1e-4), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
//		assertTrue(equals(x.getEntry(1, 0), 0.2, 1e-4), "Found incorrect solution for f(x)=xIx, 0.2 - x < 0");
//	}

	private void test_JOptimizer()
	{
		startTest();
		
		final RealMatrix Q = new Array2DRowRealMatrix(new double[][] { 
				{ 0.00687, 0 },
				{ 0, 0 } 
			});
		final RealVector c = new ArrayRealVector(new double[] { 
						-0.00687 * 17., 
						0
			});
		final double l = 0.00687 * 78.9;
		
		//you can implement the function definition using whatever linear algebra library you want, you are not tied to Colt
		StrictlyConvexMultivariateRealFunction objectiveFunction = new StrictlyConvexMultivariateRealFunction() {

			public double value(double[] X) 
			{
				ArrayRealVector x = new ArrayRealVector(X);
				return x.dotProduct(Q.operate(x)) + c.dotProduct(x) + l;
			}

			public double[] gradient(double[] X) 
			{
				ArrayRealVector x = new ArrayRealVector(X);
				return Q.operate(x).mapMultiply(2).add(c).toArray();
			}

			public double[][] hessian(double[] X) 
			{
				return Q.scalarMultiply(2).getData();
			}

			public int getDim() {
				return 2;
			}
		};

		// Inquality constraints
		ConvexMultivariateRealFunction[] inequalities = new ConvexMultivariateRealFunction[2];
		inequalities[0] = new ConvexMultivariateRealFunction()
		{
			@Override
			public double value(double[] X)
			{
				return X[0]*X[0]+X[1]*X[1]-0.01;
			}
			
			@Override
			public double[] gradient(double[] X)
			{
				return new double[] {2*X[0], 2*X[1]};
			}
			
			@Override
			public double[][] hessian(double[] X)
			{
				return new double[][]{
						{2,0},
						{0,2}
				};
			}
			
			@Override
			public int getDim()
			{
				return 2;
			}
		};
		inequalities[1] = new ConvexMultivariateRealFunction()
		{
			@Override
			public double value(double[] X)
			{
				return -X[0];
			}
			
			@Override
			public double[] gradient(double[] X)
			{
				return new double[] {-1, 0};
			}
			
			@Override
			public double[][] hessian(double[] X)
			{
				return new double[][]{
						{0,0},
						{0,0}
				};
			}
			
			@Override
			public int getDim()
			{
				return 2;
			}
		};

		OptimizationRequest or = new OptimizationRequest();
		or.setF0(objectiveFunction);
		or.setFi(inequalities);
		or.setInitialPoint(new double[]{1e-12, 1e-12});

		// optimization
		JOptimizer opt = new JOptimizer();
		opt.setOptimizationRequest(or);
		try
		{
			int returnCode = opt.optimize();
			assertTrue(returnCode==OptimizationResponse.SUCCESS, "Failed to optimise.");
			
			double[] sol = opt.getOptimizationResponse().getSolution();
			assertTrue(equals(sol[0], 0.1, 1e-6) && equals(sol[1], 0, 1e-6), "JOptimizer solution incorrect: "+sol[0]+", "+sol[1]);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			assertTrue(false, "Exception (see above).");
		}
		
		endTest("JOptimizer");
	}

	// Basic quadratic function y=(x1^-1)^2 + (x2-1)^2:
	static HessianFunction f = new HessianFunction()
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

	// Barrier function for 1.2 < x: 1.2 - x < 0
	TwiceDifferentiableFunction barrier = new TwiceDifferentiableFunction()
	{
		@Override
		public RealVector gradient(RealVector x)
		{
			double x0 = x.getEntry(0);
			return VectorHelper.vector(1/(CONSTRAINT - x0)); // grad -ln(-(1.2 - x)) = 1/(1.2-x)
		}

		@Override
		public RealMatrix hessian(RealVector x)
		{
			double x0 = x.getEntry(0);
			double d = CONSTRAINT - x0;
			return new Array2DRowRealMatrix(new double[][]{{ 1/d/d }}); // grad^2 -ln(-(1.2 - x)) = 1/(1.2-x)^2
		}

		@Override
		public double value(RealVector x)
		{
			return -Math.log(-(CONSTRAINT - x.getEntry(0)));
		}
	};

	private void test_quadratic()
	{
		startTest();
		
		BarrierSolver bs = new BarrierSolver();
		bs.setBarrier(barrier);
		bs.setFunction(f);
		bs.setConstraintDimension(1);
		bs.setStepSize(1.0);
		bs.setMultiplierZ(5);
		bs.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{2});
		RealVector x = bs.solve(x0);
		assertTrue(equals(x.getEntry(0), CONSTRAINT, 1e-3), "Minimum was not at 1.2, was at "+x.getEntry(0));
		
		endTest("quadratic");
	}
	
	private void test_logBarrierFunction()
	{
		startTest();
		
		BarrierSolver bs = new BarrierSolver();
		LogBarrierFunction barrier = new LogBarrierFunction();
		
		// Constraint 1.2 - x < 0:
		barrier.addConstraint(new TwiceDifferentiableFunction()
		{
			/**
			 * 1.2 - x
			 */
			@Override
			public double value(RealVector x)
			{
				return 1.2-x.getEntry(0);
			}
			
			/**
			 * -1
			 */
			@Override
			public RealVector gradient(RealVector x)
			{
				return VectorHelper.vector(-1.0);
			}
			
			/**
			 * 0
			 */
			@Override
			public RealMatrix hessian(RealVector x)
			{
				return new Array2DRowRealMatrix(new double[][]{{0.0}});
			}
		});

		bs.setBarrier(barrier);
		bs.setFunction(f);
		bs.setConstraintDimension(1);
		bs.setStepSize(1.0);
		bs.setMultiplierZ(5);
		bs.setEpsilon(1e-6);
		ArrayRealVector x0 = new ArrayRealVector(new double[]{2});
		RealVector x = bs.solve(x0);
		assertTrue(equals(x.getEntry(0), CONSTRAINT, 1e-3), "Minimum was not at 1.2, was at "+x.getEntry(0));
		
		endTest("quadratic");
	}
}