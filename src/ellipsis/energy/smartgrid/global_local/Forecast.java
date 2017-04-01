package ellipsis.energy.smartgrid.global_local;

public class Forecast
{
    private double irradianceMean;
    private double irradianceStandardDeviation;
    private double temperatureMean;
    private double temperatureStandardDeviation;
    private double windSpeedMean;
    private double windSpeedStandardDeviation;
    private double loadMean;
    private double loadStandardDeviation;
    private int time;
    
    public double getIrradianceMean()
    {
        return irradianceMean;
    }
    public void setIrradianceMean(double irradianceMean)
    {
        this.irradianceMean = irradianceMean;
    }
    public double getIrradianceStandardDeviation()
    {
        return irradianceStandardDeviation;
    }
    public void setIrradianceStandardDeviation(double irradianceStandardDeviation)
    {
        this.irradianceStandardDeviation = irradianceStandardDeviation;
    }
    public double getTemperatureMean()
    {
        return temperatureMean;
    }
    public void setTemperatureMean(double temperatureMean)
    {
        this.temperatureMean = temperatureMean;
    }
    public double getTemperatureStandardDeviation()
    {
        return temperatureStandardDeviation;
    }
    public void setTemperatureStandardDeviation(double temperatureStandardDeviation)
    {
        this.temperatureStandardDeviation = temperatureStandardDeviation;
    }
    public double getWindSpeedMean()
    {
        return windSpeedMean;
    }
    public void setWindSpeedMean(double windSpeedMean)
    {
        this.windSpeedMean = windSpeedMean;
    }
    public double getWindSpeedStandardDeviation()
    {
        return windSpeedStandardDeviation;
    }
    public void setWindSpeedStandardDeviation(double windSpeedStandardDeviation)
    {
        this.windSpeedStandardDeviation = windSpeedStandardDeviation;
    }
    public double getLoadMean()
    {
        return loadMean;
    }
    public void setLoadMean(double loadMean)
    {
        this.loadMean = loadMean;
    }
    public double getLoadStandardDeviation()
    {
        return loadStandardDeviation;
    }
    public void setLoadStandardDeviation(double loadStandardDeviation)
    {
        this.loadStandardDeviation = loadStandardDeviation;
    }
    
    
    //// Calculated fields ////
    
    public double getLoadLower()
    {
    	return getLoadMean() - getLoadStandardDeviation()*3;
    }

    public double getLoadUpper()
    {
    	return getLoadMean() + getLoadStandardDeviation()*3;
    }
    
    public double getIrradianceLower()
    {
    	return getIrradianceMean() - getIrradianceStandardDeviation()*3;
    }

    public double getIrradianceUpper()
    {
    	return getIrradianceMean() + getIrradianceStandardDeviation()*3;
    }
    
    
    //// PDFs ////
    
    /**
     * Assumes Gaussian distribution.
     */
	public double irradianceProbability(double irradiance) 
	{
		double d = 1/(irradianceStandardDeviation*Math.sqrt(2*Math.PI));
		double exponent = -0.5*Math.pow((irradiance-irradianceMean)/irradianceStandardDeviation, 2);
		return d*Math.exp(exponent);
	}
	public int getTime()
	{
		return time;
	}
	public void setTime(int time)
	{
		this.time = time;
	}
}
