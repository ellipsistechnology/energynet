package ellipsis.energy.scenario;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Calendar;
import java.util.Date;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.scenario.Scenario.PreAnalysisUpdater;
import ellipsis.util.ContinuousTimeProperties;
import ellipsis.util.ConverterUtil;

public class ContinuousUpdater implements PreAnalysisUpdater
{
    private ContinuousTimeProperties loads;
    private Date startTime;
    private int secondsPerUnitTime;
    private Grid grid;
    
    public ContinuousUpdater(Grid grid)
    {
        this.grid = grid;
    }
    
    @Override
    public void update(int timeIndex)
    {
        Calendar cal = Calendar.getInstance();
        cal.setTime(startTime);
        cal.add(Calendar.SECOND, timeIndex*secondsPerUnitTime);
        Date time = cal.getTime();
        
        Complex loadPower = ConverterUtil.toComplex(loads.get(time));
        setLoads(loadPower);
    }

    public void setLoads(Complex loadPower)
    {
        for (Load load : grid.getLoads())
        {
            setLoad(load, loadPower);
        }
    }
    
    public void setLoad(Load load, Complex loadPower)
    {
        load.setLoad(loadPower);
    }

    
    //// Loads ////

    public void setLoads(String propertiesPath)
    {
        setLoads(new File(propertiesPath));
    }

    private void setLoads(File propertiesFile)
    {
        try
        {
            setLoads(new FileInputStream(propertiesFile));
        } 
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }

    private void setLoads(InputStream propertiesInputStream)
    {
        loads = new ContinuousTimeProperties();
        try
        {
            loads.load(propertiesInputStream);
        } 
        catch (IOException e)
        {
            throw new RuntimeException(e);
        }
        loads.sort();
    }

    public ContinuousTimeProperties getLoads()
    {
        return loads;
    }

    public void setLoads(ContinuousTimeProperties loads)
    {
        this.loads = loads;
    }
    
    
    //// Accessors ////

    /**
     * The time at timeIndex 1.
     * @return
     */
    public Date getStartTime()
    {
        return startTime;
    }

    public void setStartTime(Date startTime)
    {
        this.startTime = startTime;
    }

    /**
     * Multiplied by the current timeIndex in {@link #update(int)} and added to the
     * startTime to find the corresponding date.
     * @return
     */
    public int getSecondsPerUnitTime()
    {
        return secondsPerUnitTime;
    }

    public void setSecondsPerUnitTime(int secondsPerUnitTime)
    {
        this.secondsPerUnitTime = secondsPerUnitTime;
    }    
}
