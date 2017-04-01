package ellipsis.energy.grid;

/**
 * A switch behaves as a line when on, and as a line with zero inductance and
 * {@link Double#MAX_VALUE} resistance when off.
 * A switch will be on by default.
 * @author bmillar
 *
 */
public class Switch extends Line
{
    private boolean on = true;
    
    public Switch()
    {
        setLength(1); // Make per metre and absolute values the same
    }

    public boolean isOn()
    {
        return on;
    }

    public void setOn(boolean on)
    {
        this.on = on;
    }

    @Override
    public double getInductancePerMetre()
    {
        if(on)
            return super.getInductancePerMetre();
        else
            return 0;
    }
    
    @Override
    public double getResistancePerMetre()
    {
        if(on)
            return super.getResistancePerMetre();
        else
            return Double.POSITIVE_INFINITY;
    }
    
    public void setResistance(double resistance)
    {
        super.setResistancePerMetre(resistance);
    }
    
    public void setInductance(double inductance)
    {
        super.setInductancePerMetre(inductance);
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	Switch copy = new Switch();
    	copyInto(copy);
    	copy.setOn(isOn());
    	return copy;
    }
}
