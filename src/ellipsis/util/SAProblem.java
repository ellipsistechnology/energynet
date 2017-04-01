package ellipsis.util;

public interface SAProblem<X>
{
	/**
	 * Retrieves the initial solution - the starting point.
	 * @return x_0
	 */
	X getInitialSolution();
	
	/**
	 * Calculates the cost of solution x. Costs should be normalised:
	 * i.e. close to 1 for a high cost, clost to 0 for a low cost.
	 * @param x The solution.
	 * @return The cost of solution x.
	 */
	double cost(X x);
	
	/**
	 * Obtains a random solution in the neighbourhood of the current solution, x.
	 * This will be called repeatedly if given solutions are not being accepted.
	 * @param x The current solution.
	 * @param attempt The attempt number; will be 0 for the first attempt, > 0 for retries.
	 * @return The next solution.
	 */
	X nextSolution(X x, int attempt);

	/**
	 * Called when the 'temperature' T has been reduced and the SASolver
	 * is about to attempt another solution.
	 * This method can be used to set the standard deviation of any random
	 * elements used to calculate the next solution; i.e. reduce the neighbourhood
	 * before the next iteration.
	 * @param T
	 */
	void prepareNextSolution(double T);
	
	/**
	 * Performs any finalisation on the solution.
	 * @param x The final solution.
	 * @return The value to be returned by the SASolver.
	 */
	X finish(X x);
}
