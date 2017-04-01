package ellipsis.energy.calculation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Source;

public class BasicPowerGenerationOptimiser
{
    private Map<String, Source> sources = new HashMap<String, Source>();

    public void addSource(Source source)
    {
        sources.put(source.getName(), source);
    }
    
    public void addAllSources(Collection<Source> sources)
    {
        for (Source source : sources)
        {
            addSource(source);
        }
    }
    
    public Source getSource(String name)
    {
        return sources.get(name);
    }
    
    public Collection<Source> getSources()
    {
        return sources.values();
    }
    
    /**
     * Sets the power generation levels of the power sources.
     * An attempt is made to set all power levels such that the generation costs are the same.
     * However, source min and max values will be respected.
     * As such there is no guarantee that the target power level will be met.
     * @param targetTotalPowerOutput The desired total power output in Watts.
     */
    public void optimise(double targetTotalPowerOutput)
    {
        Collection<Source> sources = new ArrayList<Source>(getSources());
        double remainder = targetTotalPowerOutput;
        
        // Get target cost:
        double targetCost = targetCost(targetTotalPowerOutput, sources);

        // Set power outputs:
        for (Source source : sources)
        {
            // Get target power output:
            double targetPowerOutput = targetPowerOutput(targetCost, source);
            double pmin = source.getPmin();
            double pmax = source.getPmax();
            if(targetPowerOutput < pmin)
            {
                source.setPowerOutput(new Complex(pmin));
                remainder -= pmin;
            }
            else if(targetPowerOutput > pmax)
            {
                source.setPowerOutput(new Complex(pmax));
                remainder -= pmax;
            }
            else
            {
                source.setPowerOutput(new Complex(targetPowerOutput));
                remainder -= targetPowerOutput;
            }
        }
        
        // Assign remainder:
        final double PRECISION = 0.001;
        if(remainder <= PRECISION)
            return;
        
        final int MAX_ITERATIONS = 5;
        int iterationCount = 0;
        double previousRemainder;
        do
        {
            targetCost = targetCost(remainder, sources);

            previousRemainder = remainder;
            for (Iterator<Source> iterator = sources.iterator(); iterator.hasNext();)
            {
                Source source = (Source) iterator.next();

                double adjustment = targetPowerOutput(targetCost, source);
                double pmin = source.getPmin();
                double pmax = source.getPmax();
                if(source.getPowerOutput().abs() + adjustment < pmin)
                {
                    adjustment = pmin - source.getPowerOutput().abs();
                    source.setPowerOutput(new Complex(pmin));
                }
                else if(source.getPowerOutput().abs() + adjustment > pmax)
                {
                    adjustment = pmax - source.getPowerOutput().abs();
                    source.setPowerOutput(new Complex(pmax));
                }
                else
                {
                    source.setPowerOutput(new Complex(source.getPowerOutput().getReal() + adjustment));
                }
                
                remainder -= adjustment;
                
                // Remove any sources that can't help balance the remainder any more:
                if(remainder < 0 && pmin == source.getPowerOutput().abs())
                    iterator.remove();
                else if(remainder > 0 && pmax == source.getPowerOutput().abs())
                    iterator.remove();
            }

            ++iterationCount;
        } while(remainder > PRECISION && previousRemainder != remainder && iterationCount < MAX_ITERATIONS);
    }
    
    private double targetPowerOutput(double targetCost, Source source)
    {
        return (targetCost - source.getBeta())/(2*source.getAlpha());
    }

    private double targetCost(double targetTotalPowerOutput,
            Collection<Source> sources)
    {
        /* 
         * Ct = target cost of each source
         * PT = target total combined power output
         * alpha and beta are constants of the cost equation
         * 
         * Ct = (PT + sum[beta/2alpha])/sum[1/2alpha]
         */
        double sum1 = 0;
        double sum2 = 0;
        for (Source source : sources)
        {
            sum1 += source.getBeta()/(2*source.getAlpha());
            sum2 += 1/(2*source.getAlpha());
        }
        double targetCost = (targetTotalPowerOutput + sum1)/sum2;
        return targetCost;
    }
}
