package ellipsis.energy.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;

public class CentralTemplateGridlessADP extends LCGridlessADP
{
	private int transformerCount;
	private List<GridlessTransformerRegulator> transformers = new ArrayList<>();

	public List<GridlessTransformerRegulator> getTransformers()
	{
		return transformers;
	}

	public void setTransformers(List<GridlessTransformerRegulator> transformers)
	{
		this.transformers = transformers;
	}

	public int getTransformerCount()
	{
		return transformerCount;
	}

	public void setTransformerCount(int transformerCount)
	{
		this.transformerCount = transformerCount;
	}
	
	protected int transformerPowerOffset;
	@Override
	public void finaliseData()
	{
		super.finaliseData();
		transformerPowerOffset = transformerPowerOffset();
	}
	
	@Override
	public int unitCount()
	{
		return super.unitCount() + transformerCount;
	}

	@Override
	public int loadPowerOffset()
	{
		return super.loadPowerOffset()+transformerCount;
	}

	public int transformerPowerOffset()
	{
		return super.storagePowerOffset()+gridlessData.storageCount; // after storage
	}
	
	public void addTransformer(GridlessTransformerRegulator t)
	{
		transformers.add(t);
	}
	
	@Override
	public int controlDimension()
	{
		// Control for xfmr is actually the equivalent shunt power (i.e. P&Q).
		return super.controlDimension() + transformerCount*2;
	}
	
	@Override
	public RealVector powerFromControlAndNoise(RealVector u, RealVector w)
	{
		RealVector pq = super.powerFromControlAndNoise(u, w);
		int offset = transformerPowerOffset;
		for(int i = offset; i < offset+transformerCount; ++i)
		{
			pq.setEntry(i, u.getEntry(i)); // P
			pq.setEntry(i+powerDimension/2, u.getEntry(i+controlDimension/2)); // Q
		}
		return pq;
	}
}