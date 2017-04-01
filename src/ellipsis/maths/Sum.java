package ellipsis.maths;

public class Sum extends Variable
{
	private Variable v1, v2;
	
	public Sum(Variable v1, Variable v2)
	{
		super(v1.getName()+'+'+v2.getName());
		this.v1 = v1;
		this.v2 = v2;
	}
	
	@Override
	public Variable copy()
	{
		return new Sum(v1.copy(), v2.copy());
	}
	
	@Override
	public boolean sameBase(Variable v)
	{
		if(!(v instanceof Sum))
			return false;
		
		Sum p = (Sum)v;
		return 
				(v1.sameBase(p.v1) && v2.sameBase(p.v2)) ||
				(v1.sameBase(p.v2) && v2.sameBase(p.v1));
	}
	
//	@Override
//	public Variable add(Variable v)
//	{
//		Variable vs[];
//		if(v instanceof Sum)
//			vs = new Variable[]{((Sum)v).v1, ((Sum)v).v2};
//		else
//			vs = new Variable[]{v};
//		
//		Variable sum = Variable.ZERO;
//		for (int i = 0; i < vs.length; i++)
//		{
//			if(vs[i].sameBase(v1) && vs[i].getExponent() == v1.getExponent())
//			{
//				sum = sum.add(vs[i].add(v1));
//			}
//			else if(vs[i].sameBase(v2) && vs[i].getExponent() == v2.getExponent())
//			{
//				sum = sum.add(vs[i].add(v2));
//			}
//			else
//			{
//				sum = sum.add(vs[i]);
//			}
//		}
//		
////		sum = consolodate(sum.v1, sum.v2);
//		
//		return sum;
//	}
	
	@Override
	public Variable multiply(Variable v)
	{
		Variable vs[];
		if(v instanceof Sum)
			vs = new Variable[]{((Sum)v).v1, ((Sum)v).v2};
		else
			vs = new Variable[]{v};
		
		Variable sum = Variable.ZERO;
		for (int i = 0; i < vs.length; i++)
		{
			sum = sum.add(vs[i].multiply(v1));
			sum = sum.add(vs[i].multiply(v2));
		}
		
//		sum = consolodate(sum.v1, sum.v2);
		
		return sum;
	}
	
//	private Variable consolodate(Variable v1, Variable v2)
//	{
//		if(v1 instanceof Sum && v2 instanceof Sum)
//			return consolodate((Sum)v1, (Sum)v2);
//		
//		if(v1 instanceof Sum)
//			return consolodate((Sum)v1, v2);
//		
//		if(v2 instanceof Sum)
//			return consolodate((Sum)v2, v1);
//		
//		if(v1.sameBase(v2) && v1.getExponent() == v2.getExponent())
//		{
//			Variable v = v1.copy();
//			v.setMultiple(v1.getMultiple()+v1.getMultiple());
//			return v;
//		}
//		else
//		{
//			return this;
//		}
//	}
//
//	private Variable consolodate(Sum v1, Variable v2)
//	{
//		return v2.multiply(v1.v1).add(v2.multiply(v1.v2));
//	}
//
//	private Variable consolodate(Sum v1, Sum v2)
//	{
//		return 
//	}

	@Override
	public String toString()
	{
		return v1.toString()+'+'+v2.toString();
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
			return val1+val2;
		else
			return Double.NaN;
	}
}