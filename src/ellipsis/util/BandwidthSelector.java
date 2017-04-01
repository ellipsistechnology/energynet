package ellipsis.util;

import java.util.List;
import java.util.Set;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.direct.BaseAbstractMultivariateOptimizer;
import org.apache.commons.math3.optimization.direct.NelderMeadSimplex;
import org.apache.commons.math3.optimization.direct.SimplexOptimizer;

import com.mls.util.Util;

public class BandwidthSelector
{
	public static boolean DEBUG_COST = false;
	public static final double MIN_KERNEL = 1e-200;
	public static final int BOUNDARY = 2;
	
	private class CostFunction implements MultivariateFunction
	{
		private int nearestIndex;
		private int dimension;
		private RealVector minPs;
		public double[] minPoint;
		private RealVector maxPs;
		private Set<Double>[] ignore;
		public double minCost = Double.MAX_VALUE;

		private CostFunction(
				int nearestIndex, 
				int dimension, 
				RealVector minPs,
				RealVector maxPs, 
				Set<Double>[] ignore)
		{
			this.nearestIndex = nearestIndex;
			this.dimension = dimension;
			this.minPs = minPs;
			this.maxPs = maxPs;
			this.ignore = ignore;
		}

		@Override
		public double value(double[] params)
		{
			if(DEBUG_COST)
				for(int d = 0; d < params.length; ++d)
					System.out.print(params[d]+",");

			// Penalise parameters that are outside the min and max bounds:
			double penalty = 0;

			for (int d = 0; d < params.length; d++)
			{
				double min_j = minPs.getEntry(d);
				double max_j = maxPs.getEntry(d);
				double b_j = params[d];
				if(b_j < min_j)
				{
					penalty += Math.pow(b_j-min_j, 2);
					b_j = min_j;
				}
				else if(b_j > max_j)
				{
					penalty += Math.pow(b_j-max_j, 2);
					b_j = max_j;
				}
				
				if(Double.isInfinite(penalty))
					throw new RuntimeException("Infinite penalty!");
				
				params[d] = b_j;
			}
			
			// Set the bandwidth or bandwidth function parameters:
			setBandwidth(approximator, doubleSided, dimension, new ArrayRealVector(params));
			
			double ret;
			if(local)
			{
				if(queryPoint == null)
					throw new RuntimeException("Bandwidth selection in local mode must have query point set.");
				ret = singleCvCost(approximator, nearestIndex) + penalty;
			}
			else
			{
				ret = cvCost(approximator, ignore, queryPoint) + penalty;
			}

			if(ret < minCost)
			{
				minCost = ret;
				minPoint = params;
			}
			
			if(DEBUG_COST)
				System.out.println(","+ret);
			
			return ret;
		}
	}
	
	private NonParametricApproximator approximator;
	private RealVector queryPoint;
	private boolean local = false;
	private boolean doubleSided = false;
	public boolean tuneCorrectionRadius = false;

	public BandwidthSelector(NonParametricApproximator approximator)
	{
		this.approximator = approximator;
	}

	/**
	 * 
	 * @return true if the bandwidth was tuned, false otherwise.
	 */
	public boolean crossValidate()
	{
		if(approximator.sampleCount() <= 2*BOUNDARY)
			return false;
		
		// Check cache:
		final RealVector nearest;
		final int nearestIndex;
//		if(local && queryPoint != null)
//		{
//			nearestIndex = findNearestSample(approximator, queryPoint);
//			nearest = approximator.getSamplePoints().get(nearestIndex).getKey();
//			double[] cachedBandwidth = approximator.bandwidthCache.get(nearest);
//			if(cachedBandwidth != null)
//			{
//				approximator.setBandwidth(cachedBandwidth);
//				setBandwidth(approximator, doubleSided, nearest.getDimension(), new ArrayRealVector(cachedBandwidth));
//				return true;
//			}
//		}
//		else
		{
			nearest = null;
			nearestIndex = -1;
		}
		
		// Remember bandwidth in case we have to revert:
		double[] oldBW = approximator.getBandwidth();

		// Find min and max samples:
		final List<Pair<RealVector, Double>> samples = approximator.getSamplePoints();
		int size = samples.size();
		final RealVector max = new ArrayRealVector(samples.get(0).getKey());
		final RealVector min = new ArrayRealVector(samples.get(0).getKey());
		for (int i = 1; i < size; i++)
		{
			RealVector sample = samples.get(i).getKey();
			for(int d = 0; d < sample.getDimension(); ++d)
			{
				double entry = sample.getEntry(d);
				if(entry > max.getEntry(d))
					max.setEntry(d, entry);
				else if(entry < min.getEntry(d))
					min.setEntry(d, entry);
			}
		}
		
		// Find edge cases for each dimension:
//		@SuppressWarnings("unchecked")
		final Set<Double>[] edges = null;//new Set[min.getDimension()];
//		for (int d = 0; d < edges.length; d++)
//		{
//			// Get a sorted list of samples:
//			List<Double> list = new ArrayList<Double>();
//			for (Pair<RealVector,Double> pair : samples)
//			{
//				RealVector x_i = pair.getKey();
//				double x_i_d = x_i.getEntry(d);
//				list.add(x_i_d);
//			}
//			int termCount = list.size();
//			Double[] array = list.toArray(new Double[termCount]);
//			Arrays.sort(array);
//			
//			// Add boundary terms to the ignore set:
//			edges[d] = new HashSet<Double>();
//			if(approximator.sampleCount() > 2*BOUNDARY)
//			{
//				for(int i = 0; i < BOUNDARY; ++i)
//				{
//					edges[d].add(array[i]);
//					edges[d].add(array[termCount-i-1]);
//				}
//			}
//		}
		
		// Find min and max bandwidths:
		RealVector maxBandwidth = max.subtract(min); // Starting with all inclusive bandwidth.
		removeBadValues(oldBW, maxBandwidth);
		
		final int dimension = maxBandwidth.getDimension();
		RealVector minBandwidth;
//		if(local && queryPoint != null)
//		{
//			// Get min BW for nearest point to query point:
//			minBandwidth = nearest.subtract(queryPoint);//findMinDistance(samples, nearest, dimension, ignore);
//			minBandwidth = minBandwidth.mapMultiply(1/3.0);
//			VectorHelper.abs(minBandwidth);
//		}
//		else
//		{
			minBandwidth = findBandwidthLowerBound(samples, dimension, edges, min, max);
//		}
		removeBadValues(oldBW, minBandwidth);
		
		RealVector minParams;
		RealVector maxParams;
		if(doubleSided)
		{
			maxBandwidth = maxBandwidth.append(maxBandwidth); // left and right
			minBandwidth = minBandwidth.append(minBandwidth);
			maxParams = maxBandwidth;
			minParams = minBandwidth;
		}
		else
		{
			maxParams = maxBandwidth;
			minParams = minBandwidth;
		}
		
//		if(tuneCorrectionRadius())
//		{
//			maxParams = maxParams.append(Double.MAX_VALUE);
//			minParams = minParams.append(0);
//		}

		// Minimise bandwidth:
		CostFunction cost = new CostFunction(nearestIndex, dimension, minParams, maxParams, edges);

		/* Useful for plotting BW with Display.
		for(double b = -0.5; b < 2.5; b += 0.01)
			System.out.println(b+","+cost.value(new double[]{b}));
		 */
		
		RealVector optimalParameters = optimize(cost, minParams, maxParams, 10000);

		if(optimalParameters != null && paramsContainZero(approximator, optimalParameters)) // FIXME shouldn't this try cost.minPoint?
		{
			approximator.setBandwidth(oldBW);
			return false;
		}
		
		if(optimalParameters != null)
		{
			setBandwidth(approximator, doubleSided, dimension, optimalParameters);
			if(local && queryPoint != null)
				updateBandwidthCache(nearest, optimalParameters);
		}
		else if(cost.minPoint != null)
		{
//			System.err.println("WARNING: Optimum not found.");
			
			// Default to minimum found through iterations:
			setBandwidth(approximator, doubleSided, dimension, new ArrayRealVector(cost.minPoint));
		}
		else
		{
			approximator.setBandwidth(oldBW);
			return false;
		}
		
		return true;
	}

//	private boolean tuneCorrectionRadius()
//	{
//		return tuneCorrectionRadius && approximator.isCorrectInterpolation();
//	}

	private void updateBandwidthCache(RealVector nearest, RealVector optimalBandwidth)
	{
		double[] bwArray = optimalBandwidth.toArray();
		approximator.bandwidthCache.put(nearest, bwArray);
	}

	private void removeBadValues(double[] defaults, RealVector bandwidth)
	{
		for(int i = 0; i < bandwidth.getDimension(); ++i)
		{
			double b_i = bandwidth.getEntry(i);
			if(b_i == 0 || Double.isNaN(b_i) || b_i == Double.MAX_VALUE)
				bandwidth.setEntry(i, defaults[i]);
		}
	}

	private void setBandwidth(NonParametricApproximator approx, boolean doubleSided, int dimension, RealVector params)
	{
		RealVector bandwidth;
//		if(tuneCorrectionRadius())
//		{
//			bandwidth = params.getSubVector(0, params.getDimension()-1); // Remove tunable correction radius.
//			double radius = params.getEntry(params.getDimension()-1);
//			approx.setCorrectionRadius(radius);
//		}
//		else
//		{
			bandwidth = params;
//		}
		
//		if(doubleSided)
//		{
//			double[] left = bandwidth.getSubVector(0, dimension).toArray();
//			double[] right = bandwidth.getSubVector(dimension, dimension).toArray();
//			approx.setLeftBandwidth(left);
//			approx.setRightBandwidth(right);
//		}
//		else
//		{
			approx.setBandwidth(bandwidth.toArray());
//		}
	}

	private static RealVector optimize(MultivariateFunction costFunction, RealVector minBandwidth, RealVector maxBandwidth, int maxIterations)
	{
		BaseAbstractMultivariateOptimizer<MultivariateFunction> optimizer;
//		if(minBandwidth.getDimension() < 2)
//		{
			optimizer = new SimplexOptimizer();
			((SimplexOptimizer)optimizer).setSimplex(new NelderMeadSimplex(minBandwidth.getDimension()));
//		}
//		else
//		{
			//optimizer = new BOBYQAOptimizer(minBandwidth.getDimension()+2);
//			optimizer = new CMAESOptimizer();
//		}
			// Problems? Maybe try http://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/optimization/direct/CMAESOptimizer.html
		
		PointValuePair result;
		try
		{
			result = optimizer.optimize(maxIterations, costFunction, GoalType.MINIMIZE, minBandwidth.toArray());
		}
		catch(TooManyEvaluationsException e)
		{
//			e.printStackTrace();
			return null;
		}
		if(result.getValue().doubleValue() == 0 || result.getValue().doubleValue() == 1e6 || Double.isNaN(result.getValue().doubleValue()))
			return null;
		
		double[] b = result.getPoint();
		for (int j = 0; j < b.length; j++)
		{
			double min_j = minBandwidth.getEntry(j);
			double max_j = maxBandwidth.getEntry(j);
			if(b[j] < min_j)
				b[j] = min_j;
			else if(b[j] > max_j)
				b[j] = max_j;
		}
		return new ArrayRealVector(b);
	}

	/**
	 * Ref. "A Review and Comparison of Bandwidth Selection Methods for Kernel Regression"
	 * Max K\"{o}hler, Anja Schindler, Stefan Sperlich
	 * September 2011
	 * 3.1 The Corrected ASE, page 8
	 * Definition of h_{min,G}
	 * @param ignore 
	 * @param maxX 
	 * @param minX 
	 */
	protected static RealVector findBandwidthLowerBound(
			List<Pair<RealVector, Double>> samples, 
			int dimension, 
			Set<Double>[] ignore, 
			RealVector minX, 
			RealVector maxX)
	{
		// Find the max of the min distances:
//		RealVector max = new ArrayRealVector(dimension, 0);
		double[] maxDistance = new double[dimension];
		for(int i = 0; i < samples.size(); i++)
		{
			RealVector x_i = samples.get(i).getKey();

			// If this is an edge sample ignore it:
//			boolean skip = false;
//			if(ignore != null)
//				for(int d = 0; d < dimension; ++d)
//				{
//					double range = maxX.getEntry(d) - minX.getEntry(d); // In some cases the samples will have the same value for one dimension.
//					if(range > 0 && ignore[d].contains(x_i.getEntry(d)))
//					{
//						skip = true;
//						break;
//					}
//				}
//			if(skip)
//				continue;
			
			// Find min distance:
//			RealVector min = null;//new ArrayRealVector(dimension, Double.MAX_VALUE);
			double[] minDistance = Util.doubleArray(dimension, Double.MAX_VALUE);
			for(int j = 0; j < samples.size(); j++)
			{
				RealVector x_j = samples.get(j).getKey();
				for(int d = 0; d < dimension; ++d)
				{
					double x_i_d = x_i.getEntry(d);
					double x_j_d = x_j.getEntry(d);
					double distance = x_i_d - x_j_d;
					if(distance != 0 && distance < minDistance[d])
					{
						minDistance[d] = distance;
//						min = x_j.subtract(x_i);
//						VectorHelper.abs(min);
					}
				}
			}

			// Find the max of the mins:
//			if(min != null)
//			{
			for(int d = 0; d < dimension; ++d)
			{
				if(minDistance[d] != Double.MAX_VALUE && minDistance[d] > maxDistance[d])
				{
					maxDistance[d] = minDistance[d];
//					max = min;
				}
			}
//			}
		}

		RealVector max = new ArrayRealVector(maxDistance);
		return max.mapDivide(3.0); // Divide by 3 to get the Gaussian curve 99.9% point.
		
		/*
		 * Must have MIN_KERNEL <= (e^-((max/b)^2))^d,
		 * where max is the max of the min gaps as defined above,
		 * bw is the bandwidth, and d is the bandwidth dimension.
		 * This gives:
		 * b >= max*(ln(MIN_KERNEL^-(1/d)))^-(1/2)
		 */
//		double powerTerm = 1/Math.pow(MIN_KERNEL, 1.0/dimension);
//		double logTerm = Math.log(powerTerm);
//		double rootTerm = Math.sqrt(logTerm);
//		double inverseTerm = 1/rootTerm;
//		RealVector minBW = max.mapMultiply(inverseTerm);
		
		/*
		 * Take 2:
		 * Must have MIN_KERNEL <= e^-sum_{i=1}^d{(max_i/b_i)^2},
		 * where max is the max of the min gaps as defined above,
		 * bw is the bandwidth, and d is the bandwidth dimension.
		 * This gives:
		 * min{b_i} >= max{max_i}*(d/(ln(1/MIN_KERNEL)))^0.5
		 */
//		double max_max_i = 0;
//		for(int d = 0; d < dimension; ++d)
//			max_max_i = Math.max(max_max_i, max.getEntry(d));
//		double minBW = max_max_i*Math.sqrt(dimension/(Math.log(1/MIN_KERNEL)));
//		return VectorHelper.vector(Util.doubleArray(dimension, minBW));
	}

	private boolean paramsContainZero(NonParametricApproximator approximator, RealVector parameters)
	{
		int dimension = parameters.getDimension();
//		if(tuneCorrectionRadius())
//			--dimension;
		
		for (int i = 0; i < dimension; i++)
		{
			if(parameters.getEntry(i) == 0)
				return true;
		}
		return false;
	}

	public static double cvCost(NonParametricApproximator approx, Set<Double>[] ignore, RealVector queryPoint)
	{
		List<Pair<RealVector, Double>> samples = approx.getSamplePoints();
		if(samples.isEmpty())
			return 0;
		
		// Calculate cost:
		double cost = 0;
		int count = 0;
		for(int i = 0; i < samples.size(); ++i)
		{
			// NOTE: Given that I am using vectors there is no real
			//       way to order the sample points and so using a
			//       weighting function that trims the few samples
			//       near the edge is no longer practical.
			//
			// As it turns out, I really need to do this since larger
			// bandwidths result in a softening of the edge effects
			// and so larger bandwidths give a lower cost, reducing
			// accuracy.
			boolean onBoundary = false;
			Pair<RealVector, Double> sample2 = samples.get(i);
			RealVector x_i = sample2.getKey();
			if(ignore != null)
				for(int d = 0; d < x_i.getDimension(); ++d)
					if(ignore[d].contains(x_i.getEntry(d)))
					{
						onBoundary = true;
						break;
					}
			if(onBoundary)
				continue;

			double squaredError = singleCvCost(approx, i);
			cost += squaredError;//*weight;
			++count;
		}
		
		if(count == 0)
			return Double.NaN;

		return cost/count;
	}

	protected static double singleCvCost(NonParametricApproximator approx, int i)
	{
		List<Pair<RealVector, Double>> samples = approx.getSamplePoints();
		Pair<RealVector, Double> sample = samples.get(i);
		double y_i = sample.getValue();

		double approxY = approx.valueLeaveOneOut(i);

		double squaredError = Math.pow(approxY-y_i, 2);
		return squaredError;
	}
	
	
	//// LOCAL ////
	
	public static int findNearestSample(
			NonParametricApproximator approximator,
			RealVector queryPoint)
	{
		List<Pair<RealVector, Double>> samples = approximator.getSamplePoints();
		int nearest = -1;
		double minDistance = Double.MAX_VALUE;
		for(int i = 0; i < samples.size(); ++i)
		{
			RealVector x_i = samples.get(i).getKey();
			double distance = x_i.getDistance(queryPoint);
			if(distance < minDistance)
			{
				minDistance = distance;
				nearest = i;
			}
		}
		return nearest;
	}
	
	
	//// ACCESSORS ////
	
	public RealVector getQueryPoint()
	{
		return queryPoint;
	}

	public void setQueryPoint(RealVector queryPoint)
	{
		this.queryPoint = queryPoint;
	}

	public boolean isLocal()
	{
		return local;
	}

	public void setLocal(boolean local)
	{
		this.local = local;
	}

	public boolean isDoubleSided()
	{
		return doubleSided;
	}

	public void setDoubleSided(boolean doubleSided)
	{
		this.doubleSided = doubleSided;
	}

	public NonParametricApproximator getApproximator()
	{
		return approximator;
	}

	public void setApproximator(NonParametricApproximator approximator)
	{
		this.approximator = approximator;
	}
}
