package ellipsis.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

/**
 * Ref. Approximate Dynamic Programming: Solving the Curses of Dimensionality, 8.4 Nonparametric Models/8.4.2 Kernel Regression
 * @author bmillar
 */
public class NonParametricApproximator implements SampleBasedEstimator
{
	public static boolean displayWarnings = false;
	
	private LinkedList<Pair<RealVector, Double>> history = new LinkedList<Pair<RealVector, Double>>();
	private double defaultValue;
	private double[] bandwidth; // leftBandwidth is an alias for bandwidth
	private double[] rightBandwidth;
	private double forgettingFactor = 0;
	private int historySize = 30;
	private boolean correctInterpolation = false;
	private double correctionRadius;
	private boolean useTrueCorrection = false;
	
	Map<RealVector, double[]> bandwidthCache; // Bandwidths for each sample point for local BWs - used only by BandwidthSelector
	public Map<RealVector, Integer> neighbourCountCache; // Count of samples within correctionRadius for each sample.
	private Map<RealVector, Double> rangeCache;
	RealVector min;

	RealVector max;

	private double lastSumK;
	private double lastSumKy;

	public void addSample(RealVector x, double y)
	{
		// If we have too many samples already, remove the oldest:
		if(history.size() >= historySize)
			history.poll();
		
		// Add the new sample:
		history.add(new Pair<RealVector, Double>(x, y));

		// Reset and update info on the new sample set:
		int dimension = x.getDimension();
		resetForNewSample(dimension, history.size());
	}

	protected void resetForNewSample(int dimension, int sampleCount)
	{
		// Reset optimal bandwidth cache:
		bandwidthCache = new HashMap<>();
		
		// Reset neighbour count cache:
		neighbourCountCache = new HashMap<>();
		
		// Calculate min and max samples:
		min = new ArrayRealVector(dimension, Double.MAX_VALUE);
		max = new ArrayRealVector(dimension, -Double.MAX_VALUE);
		for (Pair<RealVector, Double> sample : history)
		{
			RealVector x_i = sample.getKey();
			for (int j = 0; j < dimension; j++)
			{
				double x_i_j = x_i.getEntry(j);
				double min_j = min.getEntry(j);
				double max_j = max.getEntry(j);
				
				min.setEntry(j, Math.min(min_j, x_i_j));
				max.setEntry(j, Math.max(max_j, x_i_j));
			}
		}
		
		// Update range cache:
		if(correctInterpolation && !useTrueCorrection && sampleCount > 1)
		{
			rangeCache = new HashMap<RealVector, Double>();
			RealVector cubeSize = max.subtract(min);
			RealVector x_i = new ArrayRealVector(dimension);
			int significantDimension = 0;
			for(int i = 0; i < dimension; ++i)
			{
				if(cubeSize.getEntry(i) > 0)
					++significantDimension;
			}
			int divisions = (int) Math.pow(sampleCount*1000, 1.0/significantDimension);
			RealVector divisionCube = cubeSize.mapDivide(divisions);
			double divisionVolume = 1.0;
			for(int i = 0; i < dimension; ++i)
				divisionVolume *= divisionCube.getEntry(i);
			
			updateRanges(x_i, 0, dimension, divisions, cubeSize, divisionVolume);
		}
	}
	
	private void updateRanges(RealVector x, int d, int dimension, int divisions, RealVector cubeSize, double divisionVolume)
	{
		if(d == dimension)
		{
			RealVector nearest = nearestSample(x).getKey();
			Double area = rangeCache.get(nearest);
			if(area == null)
				area = 0.0;
			area += divisionVolume;
			rangeCache.put(nearest, area);
		}
		else
		{
			double cubeSize_d = cubeSize.getEntry(d);
				
			for(int i = 0; i <= divisions; ++i)
			{
				// Set the value at this index:
				double x_d = min.getEntry(d) + cubeSize_d*i/(double)divisions;
				x.setEntry(d, x_d);
				
				// Recurse to next dimension:
				updateRanges(x, d+1, dimension, divisions, cubeSize, divisionVolume);
				
				// Only divide if the dimension is significant:
				if(cubeSize_d == 0)
					break;
			}
		}
	}

	public double value(RealVector x)
	{
		return value(x, -1);
	}

	public double valueLeaveOneOut(int i)
	{
		return value(history.get(i).getKey(), i);
	}

	/**
	 * Finds the approximate value using the given bandwidth for the kernel
	 * and ignoring the sample point at index leaveOut.
	 * @param x
	 * @param bw
	 * @param leaveOut The index of the sample to ignore. If invalid (e.g. -1) no samples will be ignored.
	 * @return
	 */
	public double value(RealVector x, int leaveOut)
	{
		if(history.isEmpty() || x.getDimension() == 0)
			return defaultValue;
		
		double weightedSum = 0;
		double sum = 0;
		
		int i = 0;
		double size = history.size();
		for (Pair<RealVector, Double> sample : history)
		{
			if(i == leaveOut)
			{
				++i;
				continue;
			}
			
			double ff = Math.pow((i+1)/size, forgettingFactor);//1 - Math.exp(-5*(i+1)/size)
//			double[] leftBW, rightBW;
//			leftBW = getLeftBandwidth();
//			rightBW = getRightBandwidth();

			RealVector x_i = sample.getKey();

			double K_h = kernel(x, x_i)*ff;
			double v_i = sample.getValue();
			double range = range(x_i);

			if(Double.isNaN(K_h))
				Util.breakpoint();
			if(Double.isNaN(v_i))
				Util.breakpoint();

			weightedSum += K_h*v_i*range;
			sum += K_h*range;
			
			++i;
		}
		
		lastSumK = sum;
		lastSumKy = weightedSum;
		
		if(sum == 0)
		{
			if(displayWarnings)
				System.err.println("WARNING: No samples in range for query point "+VectorHelper.printVector(x));

			Pair<RealVector, Double> nearest = nearestSample(x);

			double def = nearest.getValue();
			if(DEBUG)
			{
				boolean outside = false;
				int dimension = x.getDimension();
				for(int d = 0; d < dimension; ++d)
				{
					double x_d = x.getEntry(d);
					double min_d = min.getEntry(d);
					double max_d = max.getEntry(d);
					if(x_d > max_d || x_d < min_d)
					{
						outside = true;
						break;
					}
				}
				
				if(!outside) // Only care if defaulting while within bounds since this implies a bad bandwidth.
				{
					System.out.println("DEFAULTING TO "+def);
					++defaultCount;
				}
			}

			return def;
		}
		
		return weightedSum/sum;
	}

	public double range(RealVector x_i)
	{
		if(!correctInterpolation)
			return 1;
		
		double range;

		if(useTrueCorrection)
		{
			if(x_i.getDimension() != 1)
				throw new RuntimeException("Cannot use true correction for multi-dimensional states.");
			
			double sampleX_LnR[] = leftAndRightSampleX(x_i, 0);
			if(history.size() >= 2)
			{
				range = sampleX_LnR == null                ? 1 :
					sampleX_LnR[LEFT] == -Double.MAX_VALUE ? (sampleX_LnR[RIGHT] - x_i.getEntry(0))/2 :
				    sampleX_LnR[RIGHT] == Double.MAX_VALUE ? (x_i.getEntry(0) - sampleX_LnR[LEFT])/2 :
			    	   										 (sampleX_LnR[RIGHT] - sampleX_LnR[LEFT])/2;
			}
			else
			{
				range = 1;
			}
		}
		else
		{
			range = estimateRange(x_i);
		}
		
		return range;
	}
	
	private Pair<RealVector, Double> nearestSample(RealVector x)
	{
		Pair<RealVector, Double> closest = null;
		double closestDistance = Double.MAX_VALUE;
		for (Pair<RealVector, Double> sample : getSamplePoints())
		{
			double distance = sample.getKey().getDistance(x);
			if(distance < closestDistance)
			{
				closestDistance = distance;
				closest = sample;
			}
		}
		return closest;
	}
	
	private double estimateRange(RealVector x)
	{
		// Check that count is needed and possible:
		if(!correctInterpolation || useTrueCorrection || history.size() < 2)
			return 1;
		
		Double range = rangeCache.get(x);
		if(range == null)
			return 0;
		else
			return range;
	}

	private int countInRange(RealVector x)
	{
		// Check that count is needed and possible:
		if(!correctInterpolation || useTrueCorrection || correctionRadius == 0 || history.size() < 2)
			return 1;
		
		// Check cache:
		if(neighbourCountCache.get(x) != null)
			return neighbourCountCache.get(x);
		
		// Count while within radius:
		Set<RealVector> ignore = new HashSet<>();
		RealVector nearest_i = x;
		double previousNearest = 0;
		while(previousNearest < correctionRadius)
		{
			ignore.add(nearest_i);
			Pair<RealVector, Double> nearestSample = nearestSample(x, ignore);
			if(nearestSample == null) // all are within range
				break;
			nearest_i = nearestSample.getKey();
			previousNearest = nearest_i.getDistance(x);
		}
		int count = ignore.size();
		
		// Add to cache:
		neighbourCountCache.put(x, count);

		// Return count:
		return count;
	}

	public int getNeighbourCount(RealVector x_q)
	{
		Pair<RealVector, Double> nearestSample = nearestSample(x_q, null);
		RealVector x_i = nearestSample.getKey();
		return countInRange(x_i);
	}
	
	private Pair<RealVector, Double> nearestSample(RealVector x, Set<RealVector> x_ignore)
	{
		Pair<RealVector, Double> nearest = null;
		double nearestDistance = Double.MAX_VALUE;
		for (Pair<RealVector, Double> sample : history)
		{
			RealVector x_ = sample.getKey();
			boolean cont = false;
			if(x_ignore != null)
				for (RealVector ignore : x_ignore)
				{
					if(x_.equals(ignore))
						cont = true;
				}
			if(cont)
				continue;
			
			double distance = x.getDistance(x_);
			if(distance < nearestDistance)
			{
				nearestDistance = distance;
				nearest = sample;
			}
		}
		return nearest;
	}
	
public static int defaultCount = 0;
public static boolean DEBUG = false;

	/**
	 * 
	 * @param queryX Query point.
	 * @param sampleX Sample to compare with query point.
	 * @return
	 */
	public double kernel(RealVector queryX, RealVector sampleX)
	{
		double[] leftBandwidth = getLeftBandwidth(); 
		double[] rightBandwidth = getRightBandwidth();
		
		RealVector dx = sampleX.subtract(queryX);
		
		if(leftBandwidth.length != dx.getDimension() || leftBandwidth.length != dx.getDimension())
			throw new RuntimeException("Invalid bandwidth dimension: "+leftBandwidth.length+", expecting "+dx.getDimension());

		if(rightBandwidth == null)
			rightBandwidth = leftBandwidth;

		RealVector left = new ArrayRealVector(leftBandwidth);
		RealVector right = new ArrayRealVector(rightBandwidth);
		
		// Multiply kernels for each dimension:
		double product = 1;
		for (int d = 0; d < leftBandwidth.length; d++)
		{
			// Divide each element by the respective bandwidth:
			double x_d = dx.getEntry(d);
			
			double left_d = left.getEntry(d);
			double right_d = right.getEntry(d);
			
			double b_d = x_d < 0 ? left_d : right_d;
			double exp = x_d/b_d;

			// Gaussian kernel:
			double K_h = Math.exp(-exp*exp);

			product *= K_h;
			
			if(Double.isNaN(product))
				Util.breakpoint();
		}
		
		return product;//+KERNEL_BASE;
	}

	private double[] leftAndRightSampleX(RealVector sampleX, int d)
	{
		double left = -Double.MAX_VALUE;
		double right = Double.MAX_VALUE;
		double x = sampleX.getEntry(d);
		
		for (Pair<RealVector, Double> sample : history)
		{
			RealVector x_i = sample.getKey();
			double x_i_d = x_i.getEntry(d);
			
			if(x_i_d == x)
				continue;
			
			if(x_i_d < x) // Left
			{
				left = Math.max(left, x_i_d);
			}
			else // Right
			{
				right = Math.min(right, x_i_d);
			}
		}
		
		if(left == -Double.MAX_VALUE && right == Double.MAX_VALUE)
			return null;
		
		double leftAndRight[] = new double[2];
		leftAndRight[LEFT] = left;
		leftAndRight[RIGHT] = right;

//System.out.println("Left/Right @ "+sampleX.getEntry(d)+" = "+leftAndRight[LEFT]+"/"+leftAndRight[RIGHT]);
		
		return leftAndRight;
	}
	
	private static final int LEFT = 0;
	private static final int RIGHT = 1;
	
	
	//// Accessors ////

	public double getDefaultValue()
	{
		return defaultValue;
	}

	@Override
	public void setDefaultValue(double defaultValue)
	{
		this.defaultValue = defaultValue;
	}

	/**
	 * Alias for {@link #getLeftBandwidth()}.
	 * @return
	 */
	public double[] getBandwidth()
	{
		return bandwidth;
	}

	/**
	 * Alias for {@link #setLeftBandwidth(double[])}.
	 * @param bandwidth
	 */
	public void setBandwidth(double[] bandwidth)
	{
		if(bandwidth == null)
			throw new RuntimeException("Invalid banddwidth: null");

		this.bandwidth = bandwidth;
	}

	public double[] getLeftBandwidth()
	{
		return bandwidth;
	}

	public void setLeftBandwidth(double[] bandwidth)
	{
		if(bandwidth == null || bandwidth.length == 0 || bandwidth[0] <= 0.0)
			throw new RuntimeException("Invalid banddwidth: "+(bandwidth == null ? "null" : bandwidth.length == 0 ? "length=0" : bandwidth[0]));
		this.bandwidth = bandwidth;
	}

	public double[] getRightBandwidth()
	{
		return rightBandwidth;
	}

	public void setRightBandwidth(double[] bandwidth)
	{
		if(bandwidth == null)
		{
			this.rightBandwidth = null;
			return;
		}
		
		if(bandwidth.length == 0 || bandwidth[0] <= 0)
			throw new RuntimeException("Invalid banddwidth: "+(bandwidth == null ? "null" : bandwidth.length == 0 ? "length=0" : bandwidth[0]));
		this.rightBandwidth = bandwidth;
	}

	public double getForgettingFactor()
	{
		return forgettingFactor;
	}

	public void setForgettingFactor(double forgettingFactor)
	{
		this.forgettingFactor = forgettingFactor;
	}

	public int getHistorySize()
	{
		return historySize;
	}

	public void setHistorySize(int historySize)
	{
		this.historySize = historySize;
	}

	public List<Pair<RealVector, Double>> getSamplePoints()
	{
		return this.history;
	}

	public boolean isCorrectInterpolation()
	{
		return correctInterpolation;
	}

	public void setCorrectInterpolation(boolean correctInperpolation)
	{
		this.correctInterpolation = correctInperpolation;
	}

	public int sampleCount()
	{
		return history.size();
	}
	
	public double getCorrectionRadius()
	{
		return correctionRadius;
	}

	public void setCorrectionRadius(double correctionRadius)
	{
		this.correctionRadius = correctionRadius;
		this.neighbourCountCache = new HashMap<RealVector, Integer>();
	}

	public boolean isUseTrueCorrection()
	{
		return useTrueCorrection;
	}

	public void setUseTrueCorrection(boolean useTrueCorrection)
	{
		this.useTrueCorrection = useTrueCorrection;
	}
	
	
	//// DEBUG ////

	@SuppressWarnings("unused")
	private Object showSamples = new Object(){
		public String toString() 
		{
			StringBuffer sb = new StringBuffer();
			for (Pair<RealVector, Double> sample : history)
			{
				sb.append(VectorHelper.printVector(sample.getKey()));
				sb.append(" =, ");
				sb.append(sample.getValue());
				sb.append("\n");
			}
			return sb.toString();
		};
	};

	public double getLastSumK()
	{
		return lastSumK;
	}

	public void setLastSumK(double lastSumK)
	{
		this.lastSumK = lastSumK;
	}

	public double getLastSumKy()
	{
		return lastSumKy;
	}

	public void setLastSumKy(double lastSumKy)
	{
		this.lastSumKy = lastSumKy;
	}
	
	
	//// TEST ////
	
//	public static void main(String[] args)
//	{
//		int historySize = 600;
//		double inc = 0.1;
//		
//		NonParametricApproximator approximatorCI = approximator(historySize);
//		approximatorCI.setCorrectInterpolation(true);
//		approximatorCI.setCorrectionRadius(0.1);
//		approximatorCI.setUseTrueCorrection(true);
//		NonParametricApproximator approximatorNW = approximator(historySize);
//		
//		// Train:
//		Random r = new Random(0);
//		for (int i = 0; i < historySize; i++)
//		{
//			double x = gaussian(r)*10*inc;
//			double y = f1d(x);
//			approximatorCI.addSample(new ArrayRealVector(new double[]{x}), y);
//			approximatorNW.addSample(new ArrayRealVector(new double[]{x}), y);
//		}
//		
//		// Tune:
//		new BandwidthSelector(approximatorNW).crossValidate();
//		new BandwidthSelector(approximatorCI).crossValidate();
//		new CorrectionRadiusSelector(approximatorCI).tuneCorrectionRadius();
//		
//		// Log estimate:
//		System.out.println("x,y,NW(x),CI(x),Error(NW),Error(CI)");
//		for (Pair<RealVector, Double> sample : approximatorCI.getSamplePoints())
//		{
//			double x = sample.getKey().getEntry(0);
//			double y = sample.getValue();
//			double nwY = approximatorNW.value(sample.getKey());
//			double ciY = approximatorCI.value(sample.getKey());
//			double errorNW = Math.abs(nwY-y);
//			double errorCI = Math.abs(ciY-y);
//			System.out.println(x+","+y+","+nwY+","+ciY+","+errorNW+","+errorCI);
//		}
//	}
//
//	private static double f1d(double x)
//	{
//		return (Math.sin(5*x)+Math.sin(11*x)+2)/4;
//	}
//
//	private static NonParametricApproximator approximator(int historySize)
//	{
//		NonParametricApproximator approximator = new NonParametricApproximator();
//		approximator.setDefaultValue(1);
//		approximator.setForgettingFactor(0);
//		approximator.setBandwidth(new double[]{0.1});
//		approximator.setHistorySize(historySize);
//		return approximator;
//	}
//	
//	public static void main2d(String[] args)
//	{
//		NonParametricApproximator approximator = new NonParametricApproximator();
//		approximator.setDefaultValue(1);
//		approximator.setForgettingFactor(0);
//		approximator.setBandwidth(new double[]{0.1, 0.1});
//		approximator.setCorrectInterpolation(true);
//		approximator.setCorrectionRadius(0.1);
//		int historySize = 600;
//		approximator.setHistorySize(historySize);
//		double inc = 0.1;
//		
//		// Train:
//		Random r = new Random(0);
//		System.out.println("a,b,y");
//		for (int i = 0; i < historySize; i++)
//		{
//			double a = gaussian(r)*10*inc;
//			double b = gaussian(r)*10*inc;
//			double y = f2d(a,b);
//			approximator.addSample(new ArrayRealVector(new double[]{a,b}), y);
//			System.out.println(a+","+b+","+y);
//		}
//		
//		// Tune:
//		BandwidthSelector bandwidthSelector = new BandwidthSelector(approximator);
//		bandwidthSelector.crossValidate();
//		new CorrectionRadiusSelector(approximator).tuneCorrectionRadius();
//		
//		// Log estimate:
//		System.out.println();
//		System.out.println("a,b,y");
//		for (Pair<RealVector, Double> sample : approximator.getSamplePoints())
//		{
//			double a = sample.getKey().getEntry(0);
//			double b = sample.getKey().getEntry(1);
//			double y = approximator.value(sample.getKey());
//			System.out.println(a+","+b+","+y);
//		}
//		
//		// Log comparison:
//		System.out.println();
//		System.out.print("a\\b");
//		for (double b = 0; b < 10*inc; b += inc)
//			System.out.print(","+b+","+b);
//		System.out.println();
//		for (double a = 0; a < 10*inc; a += inc)
//		{
//			System.out.print(a);
//			for (double b = 0; b < 10*inc; b += inc)
//			{
//				double realValue = f2d(a,b);
//				double y = approximator.value(new ArrayRealVector(new double[]{a,b}));
//				System.out.print(","+y+","+realValue);
//			}
//			System.out.println();
//		}
//		
//		// Log errors:
//		System.out.println();
//		for (double b = 0; b < 10*inc; b += inc)
//			System.out.print(","+b);
//		System.out.println();
//		for (double a = 0; a < 10*inc; a += inc)
//		{
//			System.out.print(a+",");
//			for (double b = 0; b < 10*inc; b += inc)
//			{
//				double ap = approximator.value(new ArrayRealVector(new double[]{a,b}));
//				double realValue = f2d(a,b);
//				System.out.print(realValue-ap);
//				System.out.print(",");
//			}
//			System.out.println();
//		}
//	}
//
//	private static double gaussian(Random r)
//	{
//		double d;
//		do
//		{
//			d = r.nextGaussian()*0.1;
//		} while(d < 0 || d > 1);
//		return d;
//	}
//
////	private static final double A = 0.8;
////	private static final double B = 0.6;
//	private static double f2d(double a, double b)
//	{
//		//return A*a+B*b*b;
////		return 2*a*a*B*Math.sin(2*b);
//		return Math.sin(6*a)*Math.sin(11*b);
//	}
}
