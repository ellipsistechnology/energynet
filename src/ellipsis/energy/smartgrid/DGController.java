package ellipsis.energy.smartgrid;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayRealVector;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.util.KalmanFilter;

public class DGController
{
    private static final boolean TEST = false;
    
    private static final int DISCRETE_RESOLUTION = 5;
    private static final int HORIZON = 2;
    private static final int MONTE_CARLO_SAMPLE_SIZE = 100;
    private static final int GAMMA = 1;
    private static final double STANDARD_DEVIATION = 0.03;
    private static final double VOLTAGE_BREACH_PENALTY = 1e12;
    
    private static final int STATE_LENGTH = 6;
    private static final int S_re = 0; // power output
    private static final int S_im = 1; // power output
    private static final int D = 2; // demand
    private static final int T = 3; // temperature
    private static final int E = 4; // irradiance
    private static final int z = 5; // wind speed
    
    private DistributedSource generator;
    private ControllableDemand demand;
    
    // Forecasts:
    private double[] forecastDemand;
    private double[] forecastTemperature;
    private double[] forecastSolarIrradiance;
    private double[] forecastWindSpeed;
    private int timeIndex;
    
    // Current state:
    private Complex gridConnectionVoltage;
    
    // Constraints:
    private double minVoltage;
    private double maxVoltage;
    
    // Estimation:
    private KalmanFilter kf = new KalmanFilter();
    
    
    //// Control ////
    
    public void regulate()
    {
if(TEST) System.out.println("\n//// Regulate(t="+timeIndex+"): ////");
        double[] x = currentState();
        Complex minU = costToGo(x, timeIndex, HORIZON).minControl;
        
        if(minU != null)
            generator.setPowerOutput(generator.getPowerOutput().add(minU));

if(TEST) System.out.println("\n//// END Regulate ////");
    }

    private double[] currentState()
    {
        double[] x = new double[STATE_LENGTH];
        Complex powerOutput = generator.getPowerOutput();
        x[S_re] = powerOutput.getReal();
        x[S_im] = powerOutput.getImaginary();
        x[D] = forecastDemand[timeIndex];
        x[T] = forecastTemperature[timeIndex];
        x[E] = forecastSolarIrradiance[timeIndex];
        x[z] = forecastWindSpeed[timeIndex];
        return x;
    }
    
    private static class CostToGo
    {
        Complex minControl;
        double minCost;
        public CostToGo(Complex minControl, double minCost)
        {
            this.minControl = minControl;
            this.minCost = minCost;
        }
    }

    /**
     * 
     * @param x The current state.
     * @return The control providing the minimum cost to go.
     */
    private CostToGo costToGo(double[] x, int timeIndex, int horizon)
    {
if(TEST) System.out.println("\nCost-to-go(k="+(timeIndex-this.timeIndex)+"):");
        
        // Cost (g):
        double g = cost(x);

        // If we are at the horizon then terminate with this cost:
        if(horizon <= 0)
            return new CostToGo(null, g);

        // Control set (U):
        List<Complex> U = viablePowerOutputs(x); // The set of possible controls.

        // DP:
        double minCostToGo = Double.MAX_VALUE;
        Complex minU = null;
        double minNextCostToGo = Double.MAX_VALUE;
        for (Complex u : U)
        {
            double nextCostToGo = GAMMA*expectation(x, u, timeIndex, horizon);
            double costToGo = g + nextCostToGo;
            if(costToGo < minCostToGo)
            {
if(TEST) System.out.println("nextCostToGo("+(timeIndex-this.timeIndex)+")="+nextCostToGo);
                minCostToGo = costToGo;
                minNextCostToGo = nextCostToGo;
                minU = u;
            }
        }
if(TEST) System.out.println("END Cost-to-go: min{cost-to-go}="+g+" + "+minNextCostToGo+"="+minCostToGo);
        return new CostToGo(minU, minCostToGo);
    }

    private Random rand = new Random(4);
    
    private double expectation(double[] x, Complex u, int timeIndex, int horizon)
    {
        // Iterate over a sample of next states and calculate the cost to go:
        double mean = 0;
//        double variance = 0;
        for (int i = 0; i < MONTE_CARLO_SAMPLE_SIZE; i++)
        {
            // Generate random variance vector as percentage of current state:
            double w[] = new double[STATE_LENGTH];
            w[S_re] = 0;// NOTE gaussianVariance(x, S_re);
            w[S_im] = 0;// NOTE gaussianVariance(x, S_im);
            w[D] = gaussianVariance(x, D, forecastDemand[timeIndex%forecastDemand.length]);
            w[T] = gaussianVariance(x, T, forecastTemperature[timeIndex%forecastTemperature.length]);
            w[E] = gaussianVariance(x, E, forecastSolarIrradiance[timeIndex%forecastSolarIrradiance.length]);
            w[z] = gaussianVariance(x, z, forecastWindSpeed[timeIndex%forecastWindSpeed.length]);
            
            // Calculate next state:
            double[] x2 = nextState(x, u, w);
            
            // Add cost:
            double cost = costToGo(x2, timeIndex+1, horizon-1).minCost;
//System.out.println(cost);
            mean = (mean*i + cost)/(i+1);
//            variance = (variance*i + Math.pow(cost-mean, 2))/(i+1);
        }

//        double sd = Math.sqrt(variance);
///*if(TEST) */System.out.println(mean+",\t"+sd);
        return mean;
    }

    private double[] nextState(double[] x, Complex u, double[] w)
    {
        double[] x2 = new double[STATE_LENGTH];
        x2[S_re] = x[S_re] + u.getReal() + w[S_re];
        x2[S_im] = x[S_im] + u.getImaginary() + w[S_im];
        x2[D] = x[D] + w[D];
        x2[T] = x[T] + w[T];
        x2[E] = x[E] + w[E];
        x2[z] = x[z] + w[z];
        
        return x2;
    }

    private double gaussianVariance(double[] x, int index, double forecastValue)
    {
        double currentValue = x[index];
        double mean = forecastValue - currentValue; // expected change
        double standardDeviation = STANDARD_DEVIATION*(currentValue+forecastValue)/2; // Approximately +/-10% of average value.
        return mean + rand.nextGaussian()*standardDeviation;
    }

    private double cost(double[] x)
    {
        double v = estimateVoltage(x).abs();
        double vCost = 
                minVoltage < v && v < maxVoltage ? 0 :
                                   /* otherwise */ Math.pow(1-v, 2)*VOLTAGE_BREACH_PENALTY;
        return vCost - Complex.valueOf(x[S_re], x[S_im]).abs();
    }

    private List<Complex> viablePowerOutputs(double[] x)
    {
        // Possible controls (discretised):
        double minP = generator.getPmin();
        double maxP = generator.getAvailablePower();
        double magnitudes[] = new double[DISCRETE_RESOLUTION];
        
        for(int i = 0; i < DISCRETE_RESOLUTION; ++i)
        {
            magnitudes[i] = minP + i*((maxP-minP)/(DISCRETE_RESOLUTION-1));
        }
        
        double maxAngle = generator.getMaxAngle();
        double minAngle = generator.getMinAngle();
        double angles[] = new double[DISCRETE_RESOLUTION];
        
        for(int i = 0; i < DISCRETE_RESOLUTION; ++i)
        {
            angles[i] = minAngle + i*((maxAngle-minAngle)/(DISCRETE_RESOLUTION-1));
        }
        
        // Filter out viable controls based on constraints:
//        double[] _x = Arrays.copyOf(x, x.length);
        ArrayList<Complex> viableSs = new ArrayList<Complex>();
        for (int i = 0; i < magnitudes.length; i++)
        {
            for (int j = 0; j < angles.length; j++)
            {
                Complex s = ComplexUtils.polar2Complex(magnitudes[i], angles[j]);
//                _x[S_re] = s.getReal();
//                _x[S_im] = s.getImaginary();
//                double v = estimateVoltage(_x).abs();
//                
//System.out.println("v =\t"+v);
//                if(v < maxVoltage && v > minVoltage)
//                {
                    Complex deltaS = new Complex(s.getReal()-x[S_re], s.getImaginary()-x[S_im]);
                    viableSs.add(deltaS);
//                }
            }
        }
        
        return viableSs;
    }
    
    
    //// Voltage Approximation ////
    
    private Complex estimateVoltage(double[] x)
    {
        // Input vector:
        ArrayRealVector vX = trainingInput(x);
        
        // Kalman Filter:
        double[] approx = kf.evaluate(vX);
        
        // Reconstruct result from array of approximations:
        return new Complex(1-approx[0], approx[1]);
    }

    private ArrayRealVector trainingInput(double[] x)
    {
        ArrayRealVector vX = new ArrayRealVector(new double[]{
                x[S_re], 
                x[S_im], 
                x[D],
                x[T],
                x[E], 
                x[z],
                x[z]*x[z],
                x[z]*x[z]*x[z],  // wind power is proportional to the cube of the wind speed
                x[z]*x[z]*x[z]*x[z],
                });
        return vX;
    }
    
    public Complex train()
    {
        // Output array:
        double[] y = new double[]{
                1-gridConnectionVoltage.getReal(), // "1 - re{v2}" makes s2=0 (v2=1) => y_re=0
                gridConnectionVoltage.getImaginary()
        };
        
        // Input vector:
        ArrayRealVector vX = trainingInput(currentState());
        
        // Kalman Filter:
        double[] approx = kf.kalmanFilter(y, vX);
        
        // Reconstruct result from array of approximations:
        return new Complex(1-approx[0], approx[1]);
    }
    
    
    //// Accessors ////

    public DistributedSource getGenerator()
    {
        return generator;
    }
    public void setGenerator(DistributedSource generator)
    {
        this.generator = generator;
    }
    public double[] getForecastDemand()
    {
        return forecastDemand;
    }
    public void setForecastDemand(double[] forecastDemand)
    {
        this.forecastDemand = forecastDemand;
    }
    public double[] getForecastTemperature()
    {
        return forecastTemperature;
    }
    public void setForecastTemperature(double[] forecastTemperature)
    {
        this.forecastTemperature = forecastTemperature;
    }
    public double[] getForecastSolarIrradiance()
    {
        return forecastSolarIrradiance;
    }
    public void setForecastSolarIrradiance(double[] forecastSolarIrradiance)
    {
        this.forecastSolarIrradiance = forecastSolarIrradiance;
    }
    public double[] getForecastWindSpeed()
    {
        return forecastWindSpeed;
    }
    public void setForecastWindSpeed(double[] forecastWindSpeed)
    {
        this.forecastWindSpeed = forecastWindSpeed;
    }
    public Complex getGridConnectionVoltage()
    {
        return gridConnectionVoltage;
    }
    public void setGridConnectionVoltage(Complex gridConnectionVoltage)
    {
        this.gridConnectionVoltage = gridConnectionVoltage;
    }
    public double getMinVoltage()
    {
        return minVoltage;
    }
    public void setMinVoltage(double minVoltage)
    {
        this.minVoltage = minVoltage;
    }
    public double getMaxVoltage()
    {
        return maxVoltage;
    }
    public void setMaxVoltage(double maxVoltage)
    {
        this.maxVoltage = maxVoltage;
    }
    public ControllableDemand getDemand()
    {
        return demand;
    }
    public void setDemand(ControllableDemand demand)
    {
        this.demand = demand;
    }

    public int getTimeIndex()
    {
        return timeIndex;
    }

    public void setTimeIndex(int timeIndex)
    {
        this.timeIndex = timeIndex;
    }
    
    
    //// FAILED THREADING EXPERIMENT! ////
    
    /*private class ExpectationTask implements Runnable
    {
        private int horizon;
        private double[] x;
        private Complex u;
        private int timeIndex;
        private Sum sum;
        
        public ExpectationTask(double[] x, Complex u, int timeIndex, int horizon, Sum sum)
        {
            super();
            this.horizon = horizon;
            this.x = x;
            this.u = u;
            this.timeIndex = timeIndex;
            this.sum = sum;
        }
        
        @Override
        public void run()
        {
            // Generate random variance vector as percentage of current state:
            double w[] = new double[STATE_LENGTH];
            w[S_re] = 0;// TODO gaussianVariance(x, S_re);
            w[S_im] = 0;// TODO gaussianVariance(x, S_im);
            w[D] = gaussianVariance(x, D, forecastDemand[timeIndex%forecastDemand.length]);
            w[T] = gaussianVariance(x, T, forecastTemperature[timeIndex%forecastTemperature.length]);
            w[E] = gaussianVariance(x, E, forecastSolarIrradiance[timeIndex%forecastSolarIrradiance.length]);
            w[z] = gaussianVariance(x, z, forecastWindSpeed[timeIndex%forecastWindSpeed.length]);
            
            // Calculate next state:
            double[] x2 = nextState(x, u, w);
            
            // Add cost:
            sum.add(costToGo(x2, timeIndex+1, horizon-1).minCost);
        }
    }
    private static class Sum
    {
        double sum;
        synchronized void add(double d)
        {
            sum += d;
        }
    }

    private static Queue<Runnable> taskQueue = new LinkedList<Runnable>();
    private static int activeCount = 0;
    static
    {
        for(int i = 0; i < 1000; ++i)
        {
            new Thread("Thread "+i){
                public void run() 
                {
                    while(true)
                    {
                        Runnable task = null;
                        synchronized (taskQueue)
                        {
                            if(!taskQueue.isEmpty())
                                task = taskQueue.poll();
                        }
                        if(task != null)
                        {
                            ++activeCount;
System.out.println(getName()+":\tactiveCount="+activeCount);
                            task.run();
                            --activeCount;
                            System.out.println(getName()+":Complete\tactiveCount="+activeCount);
                        }
                        else
                            try
                            {
                                sleep(10);
                            } catch (InterruptedException e)
                            {
                                throw new RuntimeException(e);
                            }
                    }
                };
            }.start();
        }
    }

    public static void waitForQueue(Object o)
    {
        boolean complete = false;
        while(!complete)
        {
            synchronized(taskQueue)
            {
                complete = activeCount == 0 && taskQueue.isEmpty();
            }
            if(!complete)
                try
                {
                    synchronized (o)
                    {
                        o.wait(10);
                    }
                } catch (InterruptedException e)
                {
                    throw new RuntimeException(e);
                }
        }
    }*/
}
