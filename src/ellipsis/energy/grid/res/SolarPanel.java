package ellipsis.energy.grid.res;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;

import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.Forecast;

/**
 * This simple solar panel model assumes a linear relationship between irradiance
 * and current.
 * i.e. I = aE
 * where I is current (A), a is a constant (the ratio of current to irradiance) and
 * E is the irradiance (W/m^2).
 * @author bmillar
 *
 */
public class SolarPanel extends RenewableGenerator
{   
    // Configuration parameters:
    private double voltage;
    private double currentToIrradianceRatio;
    private double pvUnitCount;
    
    // Variables:
    private double irradiance;

    // Private variables:
    private boolean powerValid = false;
    private Complex targetPowerOutput = null;
    
    @Override
    public Complex getPowerOutput()
    {
        if(!powerValid)
            updatePowerOutput();
        return super.getPowerOutput();
    }

    private void updatePowerOutput()
    {
        double availablePower = getAvailablePower();
        if(targetPowerOutput == null)
        {
            // No target set, just set as calculated:
            super.setPowerOutput(new Complex(availablePower));
        }
        else if(availablePower >= targetPowerOutput.abs())
        {
            // Got enough power to reach target:
            super.setPowerOutput(targetPowerOutput);
        }
        else
        {
            // Not enough power to reach target, so set power to available with target angle:
            double angle = targetPowerOutput.getArgument();
            super.setPowerOutput(new Complex(
                    availablePower*Math.cos(angle),
                    availablePower*Math.sin(angle)));
        }
        
        powerValid = true;
    }
    
    /**
     * Sets the target power output; actual power output will depend on
     * the available irradiance.
     * Note: Setting power output to null will set it to its
     * maximum value based on the current irradiance.
     */
    @Override
    public void setPowerOutput(Complex powerOutput)
    {
    	if(powerOutput != null)
    	{
	        double angle = powerOutput.getArgument();
	        if(angle > getMaxAngle())
	            throw new RuntimeException("angle of "+powerOutput+" exceeds maximum angle of "+getMaxAngle());
	        if(angle < getMinAngle())
	            throw new RuntimeException("angle of "+powerOutput+" exceeds minimum angle of "+getMinAngle());
    	}
        
        targetPowerOutput = powerOutput;
        powerValid = false;
    }
    
    public Complex getTargetPowerOutput()
    {
    	return targetPowerOutput;
    }
    
    @Override
    public double getAvailablePower()
    {
        return getAvailablePower(irradiance, 0);
    }

    public double getAvailablePower(double irradiance, double temperature)
    {
        return Math.min(pvUnitCount*voltage*currentToIrradianceRatio*irradiance, getPmax());
    }
    
    /**
     * Calculates the power output at the current target power's angle with the available power's magnitude
     * for the given irradiance and temperature.
     * @param irradiance
     * @param temperature
     * @return
     */
    @Override
    public Complex getExpectedPower(Forecast forecast)
    {
        double irradiance = forecast.getIrradianceMean();
        double temperature = forecast.getTemperatureMean();
        
        return getExpectedPower(irradiance, temperature);
    }

    private Complex getExpectedPower(double irradiance, double temperature)
    {
        double availablePower = getAvailablePower(irradiance, temperature);
        Complex powerOutput = targetPowerOutput;
        if(powerOutput == null)
            powerOutput = getPowerOutput();
        double angle = powerOutput.getArgument();
        return ComplexUtils.polar2Complex(availablePower,angle);
    }
    
    /**
     * Calculates the likely minimum power output based on the standard deviation of power
     * output for the given irradiance and temperature.
     * @param irradiance
     * @param temperature
     * @return The power at 3 standard deviations from the given irradiance and temperature mean values.
     */
    @Override
    public Complex getExpectedMinimumPower(Forecast forecast)
    {
        double irradianceMean = forecast.getIrradianceMean();
        double irradianceStandardDeviation = forecast.getIrradianceStandardDeviation();
//        double temperatureMean = forecast.getTemperatureMean();
//        double temperatureSD = forecast.getTemperatureStandardDeviation();
        
        double irradianceLowerBound = irradianceMean-3*irradianceStandardDeviation;
        if(irradianceLowerBound < 0)
            irradianceLowerBound = 0;
        return getExpectedPower(irradianceLowerBound, 0);
    }
    
    /**
     * Calculates the likely maximum power output based on the standard deviation of power
     * output for the given irradiance and temperature.
     * @param irradiance
     * @param temperature
     * @return The power at 3 standard deviations from the given irradiance and temperature mean values.
     */
    @Override
    public Complex getExpectedMaximumPower(Forecast forecast)
    {
        double irradianceMean = forecast.getIrradianceMean();
        double irradianceStandardDeviation = forecast.getIrradianceStandardDeviation();
//        double temperatureMean = forecast.getTemperatureMean();
//        double temperatureSD = forecast.getTemperatureStandardDeviation();
        
        double irradianceUpperBound = irradianceMean+3*irradianceStandardDeviation;
        return getExpectedPower(irradianceUpperBound, 0);
    }

    public double getVoltage()
    {
        return voltage;
    }

    public void setVoltage(double voltage)
    {
        this.voltage = voltage;
        powerValid = false;
    }

    public double getIrradiance()
    {
        return irradiance;
    }

    public void setIrradiance(double irradiance)
    {
        this.irradiance = irradiance;
        powerValid = false;
    }

    public double getCurrentToIrradianceRatio()
    {
        return currentToIrradianceRatio;
    }

    public void setCurrentToIrradianceRatio(double currentToIrradianceRatio)
    {
        this.currentToIrradianceRatio = currentToIrradianceRatio;
        powerValid = false;
    }

    public double getPvUnitCount()
    {
        return pvUnitCount;
    }

    public void setPvUnitCount(double pvUnitCount)
    {
        this.pvUnitCount = pvUnitCount;
        powerValid = false;
    }
    
    @Override
    public void setPmax(double pmax)
    {
        super.setPmax(pmax);
        powerValid = false;
    }
    
    @Override
    public void setPmin(double pmin)
    {
        super.setPmin(pmin);
        powerValid = false;
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	SolarPanel copy = new SolarPanel();
    	copyInto(copy);
    	copy.setVoltage(getVoltage());
    	copy.setCurrentToIrradianceRatio(getCurrentToIrradianceRatio());
    	copy.setPvUnitCount(getPvUnitCount());
    	copy.setIrradiance(getIrradiance());
    	copy.targetPowerOutput = this.targetPowerOutput;
    	return copy;
    }
}
