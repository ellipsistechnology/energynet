package ellipsis.util;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;

import com.mls.util.Util;

public class NewtonSolver extends GradientDecent
{
	public NewtonSolver()
	{
		setStep(1.0);
	}
	
	// Accessible for debug/test purposes only:
	public RealMatrix hessian;
	
	@Override
	public RealVector step(RealVector x)
	{
		grad = function.gradient(x);
		hessian = ((HessianFunction)function).hessian(x);
		RealVector step;
		try
		{
//			verbose(grad.getEntry(0)+","+hessian.getEntry(0, 0)+",");
			
			// Solve H(f(x))p=grad(f(x)); p = H^-1(f(x))*grad(f(x)) = step
			DecompositionSolver solver = new LUDecomposition(hessian).getSolver();
			RealVector fullStep = solver.solve(grad);
			step = fullStep.mapMultiply(stepSize);
		}
		catch(SingularMatrixException e)
		{
			if(verbose)
				System.err.println("Warning: Hessian was singular when calculating step size.");
			return null;
		}
		
		if(step.isNaN())
			Util.breakpoint();
		return step;
	}
	
	@Override
	public void setFunction(GradientFunction function)
	{
		if(!(function instanceof HessianFunction))
			throw new RuntimeException("Function must be a NewtonSolver.NSFunction.");
		
		super.setFunction(function);
	}
}
