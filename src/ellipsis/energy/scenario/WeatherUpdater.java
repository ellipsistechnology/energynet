package ellipsis.energy.scenario;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.LowVoltageGroup;
import ellipsis.energy.grid.LowVoltageGroup.VirtualSource;
import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.grid.res.WindTurbine;
import ellipsis.energy.scenario.Scenario.PreAnalysisUpdater;

/**
 * NOTE Allow data to be read straight from properties files.
 * @author bmillar
 *
 */
public class WeatherUpdater implements PreAnalysisUpdater
{
    private int startIndex = 0;
    private List<Double> windData;
    private List<Double> irradianceData;
    private Collection<SolarPanel> solarPanels = new HashSet<SolarPanel>();
    private Collection<WindTurbine> windTurbines = new HashSet<WindTurbine>();
    private Collection<VirtualSource> virtualSources = new HashSet<VirtualSource>();
    
    @Override
    public void update(int timeIndex)
    {
        Double irradiance = irradianceData.get(timeIndex-startIndex);
        if(irradiance == null && !solarPanels.isEmpty())
            throw new NullPointerException("Irradiance at time index "+timeIndex+" was null");
        
        Double windSpeed = windData.get(timeIndex-startIndex);
        if(windSpeed == null && !windTurbines.isEmpty())
            throw new NullPointerException("Wind speed at time index "+timeIndex+" was null");
        
        // Update Solar Panels:
        for (SolarPanel pv : solarPanels)
        {
            pv.setIrradiance(irradiance);
        }
        
        // Update Wind Turbines:
        for (WindTurbine wt : windTurbines)
        {
            wt.setWindSpeed(windSpeed);
        }
        
        // Update Low Voltage Groups:
        for (VirtualSource virtualSource : virtualSources)
        {
            LowVoltageGroup lvg = virtualSource.getLowVoltageGroup();
            if(irradiance != null)
                lvg.setIrradiance(irradiance);
            if(windSpeed != null)
                lvg.setWindSpeed(windSpeed);
        }
    }

    public List<Double> getWindData()
    {
        return windData;
    }

    public void setWindData(List<Double> windData)
    {
        this.windData = windData;
    }

    public List<Double> getIrradianceData()
    {
        return irradianceData;
    }

    public void setIrradianceData(List<Double> irradianceData)
    {
        this.irradianceData = irradianceData;
    }

    public boolean add(SolarPanel arg0)
    {
        return solarPanels.add(arg0);
    }

    public boolean add(WindTurbine arg0)
    {
        return windTurbines.add(arg0);
    }

    public boolean add(VirtualSource arg0)
    {
        return virtualSources.add(arg0);
    }

    /**
     * Adds all the {@link WindTurbine}s and {@link SolarPanel}s from the given collection.
     * @param sources
     */
    public void addAll(Collection<DistributedSource> sources)
    {
        for (DistributedSource source : sources)
        {
            if(source instanceof SolarPanel)
                add((SolarPanel)source);
            else if(source instanceof WindTurbine)
                add((WindTurbine)source);
            else if(source instanceof VirtualSource)
                add((VirtualSource)source);
        }
    }

    public int getStartIndex()
    {
        return startIndex;
    }

    public void setStartIndex(int startIndex)
    {
        this.startIndex = startIndex;
    }
}
