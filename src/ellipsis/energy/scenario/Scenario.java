package ellipsis.energy.scenario;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Source;

public class Scenario
{
    /* Scenario parameters */
    
    private Grid grid;
    private int startTimeIndex;
    private int endTimeIndex;
    //private long period; // seconds between consecutive time indices
    private double basePower;
    private double baseVoltage;
    private int iterations;
    private double targetError;
    private boolean reanalyseAfterRegulation = true;
    private List<Regulator> regulators = new ArrayList<Regulator>();
    private List<PreAnalysisUpdater> preAnalysisUpdaters = new ArrayList<Scenario.PreAnalysisUpdater>();
    
    
    /* Scenario execution variables */
    
    private int currentTimeIndex;
    
    // NOTE replace these with a new class that will give the map for a given time rather than a given index
    //      by default this will simply be a wrapper for the list of maps returning the nearest map for the
    //      given time.
    private List<Map<Load, Complex>> loadProfile;
    private List<Map<Source, Complex>> sourceProfile;
    
    
    public static interface AnalysisCallback
    {
        void handleResults(int timeIndex, AnalysisResults results);
    }
    
    public static interface PreAnalysisUpdater
    {
        void update(int timeIndex);
    }
    
    
    /* Scenario execution */
    
    /**
     * Resets the scenario to the first time index.
     */
    public void init()
    {
        // Check all parameters are valid:
        if(startTimeIndex > endTimeIndex || startTimeIndex < 0)
            throw new RuntimeException("Start time index must be >= 0 and < end time index; start time index = "+
                                        startTimeIndex+", end time index = "+endTimeIndex);
        if(basePower <= 0)
            throw new RuntimeException("Base power must be > 0; base power = "+basePower);
        if(baseVoltage <= 0)
            throw new RuntimeException("Base voltage must be > 0; base voltage = "+baseVoltage);
        
        // Init:
        currentTimeIndex = startTimeIndex-1;
    }
    
    /**
     * Updates the loads and sources for the next time index.
     * @return true if the are more executions to follow, false otherwise
     */
    public boolean next()
    {
        if(isComplete())
            throw new RuntimeException("Call to next was made on completed scenario");
        ++currentTimeIndex;
        
        update();
        
        return !isComplete();
    }

    public boolean isComplete()
    {
        return currentTimeIndex >= endTimeIndex;
    }

    private void update()
    {
        int i = currentTimeIndex-startTimeIndex;
        
        // Update loads:
        if(loadProfile != null)
            for (Load load : loadProfile.get(i).keySet())
            {
                load.setLoad(loadProfile.get(i).get(load));
            }
        
        // Update sources:
        if(sourceProfile != null)
            for (Source source : sourceProfile.get(i).keySet())
            {
                source.setPowerOutput(sourceProfile.get(i).get(source));
            }
        
        // Execute pre-analysis updates:
        for (PreAnalysisUpdater updater : preAnalysisUpdaters)
        {
            updater.update(currentTimeIndex);
        }
    }
    
    /**
     * Executes the scenario from the start time index to the end time index.
     * Each time index will have a LoadFlowAnalyser.analyse() performed and
     * the results will be passed into the callback.
     * Callback execution is synchronous. 
     * @param callback The callback will have handleResults() called for each
     * time index.
     */
    public void execute(AnalysisCallback callback)
    {
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(basePower);
        lfa.setBaseVoltage(baseVoltage);
        if(iterations != 0)
            lfa.setIterations(iterations);
        if(targetError != 0)
            lfa.setTargetError(targetError);
        
        init();
        boolean more = true;
        while(more)
        {
            // Update and analyse:
            more = next();
            AnalysisResults results = lfa.analyse();

            // Regulate:
            if(!regulators.isEmpty())
            {
                // Apply regulators:
                for (Regulator regulator : regulators)
                {
                    regulator.regulate(results);
                }
                
                // Reanalyse:
                if(this.reanalyseAfterRegulation)
                    results = lfa.analyse();
            }
            
            // Handle results:
            callback.handleResults(currentTimeIndex, results);
        }
    }
    
    /**
     * Executes the scenario from the start time index to the end time index.
     * Each time index will have a LoadFlowAnalyser.analyse() performed and
     * the results will be added to the returned list.
     * @return A list of results, one form each time index.
     */
    public List<AnalysisResults> execute()
    {
        final List<AnalysisResults> resultsList = new ArrayList<AnalysisResults>();
        execute(new AnalysisCallback()
        {
            @Override
            public void handleResults(int timeIndex, AnalysisResults results)
            {
                resultsList.add(results);
            }
        });
        
        return resultsList;
    }
    
    
    /* Accessors */

    public Grid getGrid()
    {
        return grid;
    }

    public void setGrid(Grid grid)
    {
        this.grid = grid;
    }

    public int getStartTimeIndex()
    {
        return startTimeIndex;
    }

    public void setStartTimeIndex(int startTimeIndex)
    {
        this.startTimeIndex = startTimeIndex;
    }

    public int getEndTimeIndex()
    {
        return endTimeIndex;
    }

    public void setEndTimeIndex(int endTimeIndex)
    {
        this.endTimeIndex = endTimeIndex;
    }

    public double getBasePower()
    {
        return basePower;
    }

    public void setBasePower(double basePower)
    {
        this.basePower = basePower;
    }

    public double getBaseVoltage()
    {
        return baseVoltage;
    }

    public void setBaseVoltage(double baseVoltage)
    {
        this.baseVoltage = baseVoltage;
    }

    /**
     * @return The current time index as at the previous call to next.
     * After calling {@link #init()} current time index will be start index - 1.
     */
    public int getCurrentTimeIndex()
    {
        return currentTimeIndex;
    }

    public List<Map<Load, Complex>> getLoadProfile()
    {
        return loadProfile;
    }

    public void setLoadProfile(List<Map<Load, Complex>> loadProfile)
    {
        this.loadProfile = loadProfile;
    }

    public List<Map<Source, Complex>> getSourceProfile()
    {
        return sourceProfile;
    }

    public void setSourceProfile(List<Map<Source, Complex>> sourceProfile)
    {
        this.sourceProfile = sourceProfile;
    }

    public int getIterations()
    {
        return iterations;
    }

    public void setIterations(int iterations)
    {
        this.iterations = iterations;
    }

    public double getTargetError()
    {
        return targetError;
    }

    public void setTargetError(double targetError)
    {
        this.targetError = targetError;
    }

    public List<Regulator> getRegulators()
    {
        return regulators;
    }

    public void setRegulators(List<Regulator> regulators)
    {
        this.regulators = regulators;
    }

    public boolean isReanalyseAfterRegulation()
    {
        return reanalyseAfterRegulation;
    }

    public void setReanalyseAfterRegulation(boolean reanalyseAfterRegulation)
    {
        this.reanalyseAfterRegulation = reanalyseAfterRegulation;
    }

    public void addRegulator(Regulator dso)
    {
        regulators.add(dso);
    }

    public void removeRegulator(Regulator dso)
    {
        regulators.remove(dso);
    }

    public List<PreAnalysisUpdater> getPreAnalysisUpdaters()
    {
        return preAnalysisUpdaters;
    }

    public void setPreAnalysisUpdaters(List<PreAnalysisUpdater> preAnalysisUpdaters)
    {
        this.preAnalysisUpdaters = preAnalysisUpdaters;
    }

    public void addPreAnalysisUpdater(PreAnalysisUpdater updater)
    {
        preAnalysisUpdaters.add(updater);
    }

    public void removePreAnalysisUpdater(PreAnalysisUpdater updater)
    {
        preAnalysisUpdaters.remove(updater);
    }
}
