package ellipsis.energy.grid;

import java.util.Collection;
import java.util.HashSet;

import org.apache.commons.math3.complex.Complex;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.grid.res.WindTurbine;

public class LowVoltageGroup extends Unit implements Child
{
    /* Source and Load wrapper classes */
    
    public final class VirtualSource extends DistributedSource
    {
        private static final String EXPLICIT_ACCESS_EXCEPTION = "LowVoltageGroup source cannot be explicitly set; update individual sources via LowVoltageGroup.sources instead";
        
        private VirtualSource() {}
        
        public LowVoltageGroup getLowVoltageGroup()
        {
            return LowVoltageGroup.this;
        }

        void update()
        {
            // Update power output:
            Complex powerOutput = Complex.ZERO;
            for (Source source : sources)
            {
                powerOutput = powerOutput.add(source.getPowerOutput());
            }
            super.setPowerOutput(powerOutput);
            
            // Update maximum power:
            double pMax = 0;
            for (Source source : sources)
            {
                pMax += source.getPmax();
            }
            super.setPmax(pMax);
            
            // Update minimum power:
            double pMin = 0;
            for (Source source : sources)
            {
                pMin += source.getPmax();
            }
            super.setPmin(pMin);
            
            sourcesInvalid = false;
        }

        @Override
        public Complex getPowerOutput()
        {
            if(sourcesInvalid)
                update();
            return super.getPowerOutput();
        }

        @Override
        public double getPmax()
        {
            if(sourcesInvalid)
                update();
            return super.getPmax();
        }

        @Override
        public double getPmin()
        {
            if(sourcesInvalid)
                update();
            return super.getPmin();
        }

        @Override
        public void setPowerOutput(Complex powerOutput)
        {
            throw new UnsupportedOperationException(EXPLICIT_ACCESS_EXCEPTION);
        }
        
        @Override
        public void setPmax(double pmax)
        {
            throw new UnsupportedOperationException(EXPLICIT_ACCESS_EXCEPTION);
        }
        
        @Override
        public void setPmin(double pmin)
        {
            throw new UnsupportedOperationException(EXPLICIT_ACCESS_EXCEPTION);
        }
        
        @Override
        public void setReactivePowerOutput(double Q)
        {
            throw new UnsupportedOperationException("TODO: Not yet implemented."); // TODO Not yet implemented.
        }
    }
    
    public final class VirtualLoad extends Load
    {   
        private VirtualLoad(){}
        
        public LowVoltageGroup getLowVoltageGroup()
        {
            return LowVoltageGroup.this;
        }

        @Override
        public Complex getLoad()
        {
            if(loadsInvalid)
                update();
            return super.getLoad();
        }

        private void update()
        {
            Complex loadPower = Complex.ZERO;
            for (Load load : loads)
            {
                loadPower = loadPower.add(load.getLoad());
            }
            super.setLoad(loadPower);
            
            loadsInvalid = false;
        }
        
        @Override
        public void setLoad(Complex loadPower)
        {
            Complex powerPerLoad = loadPower.divide(loads.size());
            for (Load load : loads)
            {
                load.setLoad(powerPerLoad);
            }
            loadsInvalid = true;
        }
    }
    
    
    /* Member Variables */

    private Collection<Source> sources = new HashSet<Source>();
    private Collection<Load> loads = new HashSet<Load>();
    boolean sourcesInvalid = true;
    boolean loadsInvalid = true;
    
    private Bus parent;
    private Source source = new VirtualSource();
    private Load load = new VirtualLoad();
    
    
    /* Parent Bus Managers */
    
    @Override
    public void wasRemovedFromBus(Bus parent)
    {
        parent.removeChild(source);
        parent.removeChild(load);
    }
    
    @Override
    public void setParent(Bus parent)
    {
        this.parent = parent;
        parent.addChild(source);
        parent.addChild(load);
    }
    
    
    /* Source and Load Managers */

    public void setIrradiance(double irradiance)
    {
        for (Source source : sources)
        {
            if(source instanceof SolarPanel)
                ((SolarPanel)source).setIrradiance(irradiance);
        }
        
        sourcesInvalid = true;
    }

    public void setWindSpeed(double windSpeed)
    {
        for (Source source : sources)
        {
            if(source instanceof WindTurbine)
                ((WindTurbine)source).setWindSpeed(windSpeed);
        }
        
        sourcesInvalid = true;
    }

    public void setAverageLoad(Complex loadPower)
    {
        for (Load load : loads)
        {
            load.setLoad(loadPower);
        }
        
        loadsInvalid = true;
    }
    
    
    /* Accessors */
    
    public void add(Source source)
    {
        sources.add(source);
        sourcesInvalid = true;
    }
    
    public void add(Load load)
    {
        loads.add(load);
        loadsInvalid = true;
    }

    @Override
    public Bus getParent()
    {
        return parent;
    }
    
    public Source getSource()
    {
        return source;
    }
    
    public Load getLoad()
    {
        return load;
    }

    @Override
    public void setName(String name)
    {
        super.setName(name);
        source.setName(name+":source");
        load.setName(name+":load");
    }

    public Collection<Source> getSources()
    {
        return sources;
    }

    public Collection<Load> getLoads()
    {
        return loads;
    }
    
    public int solarPanelCount()
    {
        int count = 0;
        for (Source source : sources)
        {
            if(source instanceof SolarPanel)
                ++count;
        }
        return count;
    }
    
    public int windTurbineCount()
    {
        int count = 0;
        for (Source source : sources)
        {
            if(source instanceof WindTurbine)
                ++count;
        }
        return count;
    }
    
    
    /* Factory Methods */
    
    public static LowVoltageGroup makeGroup(String name, int solarPanelCount, int windTurbineCount, Complex loadPower, int houses)
    {
        LowVoltageGroup group = new LowVoltageGroup();
        group.setName(name);
        
        // PVs:
        for (int i = 0; i < solarPanelCount; i++)
        {
            SolarPanel pv = new SolarPanel();
            pv.setName(name+":PV:"+(i+1));
            pv.setCurrentToIrradianceRatio(0.0075);
            pv.setVoltage(24);
            pv.setPvUnitCount(8);
            pv.setIrradiance(500);
            group.add(pv);
        }
        
        // Add 3kW wind turbines:
        // Ref. http://www.energymatters.com.au/southwest-windpower-48volt-3000watt-wind-turbine-p-445.html?zenid=maugskr0tqv4p1uv7spnoecv65
        for(int i = 0; i < windTurbineCount; ++i)
        {
            WindTurbine turbine = new WindTurbine();
            turbine.setName(name+":Wind:"+(i+1));
            turbine.setMaximumWindSpeed(10.5);
            turbine.setMinimumWindSpeed(3.4);
            turbine.setPmax(3000);
            turbine.setPowerCoefficient(0.266);
            turbine.setRadius(4.5/2);
            turbine.calculateLamda(15, 10);
            group.add(turbine);
        }
        
        // Add a load per house:
        if(loadPower != null)
        {
            for (int i = 0; i < houses; i++)
            {
                Load load = new Load();
                load.setLoad(loadPower);
                group.add(load);
            }
        }
        
        return group;
    }

	@Override
	public Unit copy(Grid grid) 
	{
		throw new NotImplementedException();
	}

	@Override
	public Complex netPower() 
	{
		throw new NotImplementedException();
	}
    
    /* 2MW
     * Ref. http://swanenergy.com.au/1-65-2-5mw-commercial-wind-turbine/
    for(int i = 0; i < windTurbineCount; ++i)
    {
        WindTurbine turbine = new WindTurbine();
        turbine.setName("W675");
        turbine.setMaximumWindSpeed(20);
        turbine.setMinimumWindSpeed(3.0);
        turbine.setPmax(2000000);
        turbine.setPowerCoefficient(0.5);
        turbine.setRadius(45);
        turbine.calculateLamda(15, 10);
        group.add(turbine);
    }*/
}
