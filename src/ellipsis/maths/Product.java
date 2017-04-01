package ellipsis.maths;

public class Product extends Variable
{
	private Variable v1, v2;
	
	public Product(Variable v1, Variable v2)
	{
		super(v1.getName()+"."+v2.getName());
		this.v1 = v1;
		this.v2 = v2;
	}
	
	@Override
	public boolean sameBase(Variable v)
	{
		if(!(v instanceof Product))
			return false;
		
		Product p = (Product)v;
		return 
				(v1.sameBase(p.v1) && v2.sameBase(p.v2)) ||
				(v1.sameBase(p.v2) && v2.sameBase(p.v1));
	}
	
	@Override
	public Variable copy()
	{
		return new Product(v1.copy(), v2.copy());
	}
	
	@Override
	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		
		if(getMultiple() != 1)
			sb.append(getMultiple());
		
		if(getExponent() != 1 || getMultiple() != 1)
			sb.append('(');
		
		if(v1 instanceof Sum)
		{
			sb.append("(");
			sb.append(v1.toString());
			sb.append(")");
		}
		else
		{
			sb.append(v1.toString());
		}
		sb.append('.');
		if(v2 instanceof Sum)
		{
			sb.append("(");
			sb.append(v2.toString());
			sb.append(")");
		}
		else
		{
			sb.append(v2.toString());
		}

		if(getExponent() != 1 || getMultiple() != 1)
			sb.append(')');
		
		if(getExponent() != 1)
		{
			sb.append('^');
			sb.append(getExponent());
		}
		
		return sb.toString();
	}

	public Variable getV1()
	{
		return v1;
	}

	public void setV1(Variable v1)
	{
		this.v1 = v1;
	}

	public Variable getV2()
	{
		return v2;
	}

	public void setV2(Variable v2)
	{
		this.v2 = v2;
	}
	
	@Override
	public Double getValue()
	{
		Double val1 = v1.getValue();
		Double val2 = v2.getValue();
		if(val1 != null && val2 != null)
			return val1*val2;
		else
			return Double.NaN;
	}
}