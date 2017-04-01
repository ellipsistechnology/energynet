package ellipsis.util;

import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import ellipsis.util.PrimalDualInteriorPointSolver.Constraints;

/**
 * Bx < d => Bx - d < 0
 * @author bmillar
 *
 */
public class LinearConstraints implements Constraints
{
	private Matrix B;
	private Matrix d;

	public LinearConstraints(Matrix B, Matrix d)
	{
		this.B = B;
		this.d = d;
	}

	@Override
	public Matrix derivative(Matrix x)
	{
		return B;
	}

	@Override
	public Matrix value(Matrix x)
	{
		return B.mtimes(x).minus(d);
	}

	@Override
	public Matrix hessian(int i, Matrix x)
	{
		int dim = x.getRowDimension();
		return new DenseMatrix(dim, dim); // zero
	}

	@Override
	public int constrantCount()
	{
		return 1;
	}
}