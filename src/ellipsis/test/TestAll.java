package ellipsis.test;

public class TestAll
{
	public static void main(String[] args)
	{
		System.out.println("Gradient Decent:");
		GradientDecentTest.main(args);
		System.out.println("\nNewton's Method Test:");
		NewtonSolverTest.main(args);
		System.out.println("\nBarrier Method Test:");
		BarrierSolverTest.main(args);
		System.out.println("\nLinear Least Squares Test:");
		LinearLeastSquaresTest.main(args);
		System.out.println("\nQuadratic Estimator Test:");
		QuadraticEstimatorTest.main(args);
	}
}