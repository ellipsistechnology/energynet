package ellipsis.util;

import la.matrix.Matrix;
import ml.utils.Matlab;
import ellipsis.util.PrimalDualInteriorPointSolver.Objective;

/**
 * For use with {@link PrimalDualInteriorPointSolver}.
 * f(x) = xQx + lx + s
 * @author bmillar
 *
 */
public class QuadraticObjective implements Objective
{
	private Matrix Q;
	private Matrix l;
	private double s;

	/**
	 * @param Q Quadratic matrix (xQx)
	 * @param l Linear vector (lx)
	 * @param s Scalar.
	 */
	public QuadraticObjective(Matrix Q, Matrix l, double s)
	{
		this.Q = Q;
		this.l = l;
		this.s = s;
	}

	@Override
	public double value(Matrix x)
	{
		Matrix Qx = Q.mtimes(x);
		double lx = Matlab.innerProduct(l, x);
		double fval = Matlab.innerProduct(x, Qx) + lx + s;
		return fval;
	}

	@Override
	public Matrix hessian(Matrix x)
	{
		return Q.times(2);
	}

	@Override
	public Matrix gradient(Matrix x)
	{
		return Q.mtimes(x).times(2).plus(l); // 2Qx - l
	}

}
