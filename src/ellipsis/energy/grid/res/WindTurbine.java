package ellipsis.energy.grid.res;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.Forecast;

public class WindTurbine extends RenewableGenerator
{
    // Parameters with defaults:
    private double airDensity = 1.225;
    private double powerCoefficient = 0.59;
    
    // Calculated parameters (read-only):
    private double area;
    
    // Configuration parameters:
    private double radius;
    private double lamda;
    private double minimumWindSpeed;
    private double maximumWindSpeed;
    
    // Variables:
    private double windSpeed;
    
    // Private variables:
    private boolean powerValid = false;
    
    @Override
    public Complex getPowerOutput()
    {
        if(!powerValid)
            updatePowerOutput();
        return super.getPowerOutput();
    }
    
    private void updatePowerOutput()
    {
        if(windSpeed < minimumWindSpeed)
        {
            setPowerOutput(new Complex(getPmin()));
        } 
        else
        {
            double pmax = getPmax();
            if(windSpeed > maximumWindSpeed)
            {
                setPowerOutput(new Complex(pmax));
            }
            else
            {
                double p = 0.5*airDensity*area*windSpeed*windSpeed*windSpeed*powerCoefficient;
                if(p > pmax)
                    setPowerOutput(new Complex(pmax));
                else
                    setPowerOutput(new Complex(p));
            }
        }
        powerValid = true;
    }

    public void calculateLamda(double ratedRPM, double ratedWindSpeed)
    {
        double ratedTipSpeed = ratedRPM/60*2*Math.PI*radius;
        setLamda(ratedTipSpeed/ratedWindSpeed);
    }

    public double getAirDensity()
    {
        return airDensity;
    }

    public void setAirDensity(double airDensity)
    {
        this.airDensity = airDensity;
        powerValid = false;
    }

    public double getPowerCoefficient()
    {
        return powerCoefficient;
    }

    public void setPowerCoefficient(double powerCoefficient)
    {
        this.powerCoefficient = powerCoefficient;
        powerValid = false;
    }

    public double getRadius()
    {
        return radius;
    }

    public void setRadius(double radius)
    {
        this.radius = radius;
        this.area = Math.PI*radius*radius;
        powerValid = false;
    }

    public double getLamda()
    {
        return lamda;
    }

    public void setLamda(double lamda)
    {
        this.lamda = lamda;
        powerValid = false;
    }
    
    @Override
    public void setPmax(double pmax)
    {
        super.setPmax(pmax);
        powerValid = false;
    }

    public double getMinimumWindSpeed()
    {
        return minimumWindSpeed;
    }

    public void setMinimumWindSpeed(double minimumWindSpeed)
    {
        this.minimumWindSpeed = minimumWindSpeed;
        powerValid = false;
    }

    public double getMaximumWindSpeed()
    {
        return maximumWindSpeed;
    }

    public void setMaximumWindSpeed(double maximumWindSpeed)
    {
        this.maximumWindSpeed = maximumWindSpeed;
        powerValid = false;
    }

    public double getArea()
    {
        return area;
    }

    public double getWindSpeed()
    {
        return windSpeed;
    }

    public void setWindSpeed(double windSpeed)
    {
        this.windSpeed = windSpeed;
        powerValid = false;
    }

    @Override
    public Complex getExpectedPower(Forecast forecast)
    {
        return null;
    }

    @Override
    public Complex getExpectedMinimumPower(Forecast forecast)
    {
        return null;
    }

    @Override
    public Complex getExpectedMaximumPower(Forecast forecast)
    {
        return null;
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	WindTurbine copy = new WindTurbine();
    	copyInto(copy);
    	copy.setAirDensity(getAirDensity());
    	copy.setLamda(getLamda());
    	copy.setMaximumWindSpeed(getMaximumWindSpeed());
    	copy.setMinimumWindSpeed(getMinimumWindSpeed());
    	copy.setPowerCoefficient(getPowerCoefficient());
    	copy.setRadius(getRadius());
    	copy.setWindSpeed(getWindSpeed());
    	return copy;
    }
}
