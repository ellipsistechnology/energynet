package ellipsis.maths;

import java.text.DecimalFormat;



public class Variable
{
	public static final Variable ZERO = new Variable(0);
	
	private String name;
	private Double value;
	private double exponent = 1;
	private double multiple = 1;
	
//	private static Hashtable<String, Variable> variables = new Hashtable<>();
	public static Variable variable(String name)
	{
		Variable v;
//		v = variables.get(name);
//		if(v != null)
//			return v;
		v = new Variable(name);
//		variables.put(name, v);
		return v;
	}
	
	public Variable(String name)
	{
		this.name = name;
	}
	
	public Variable(double d)
	{
		this(""+d);
		setValue(d);
	}
	
	public boolean sameBase(Variable v)
	{
		if(v instanceof Product || v instanceof Sum)
			return v.sameBase(this);
		
		return name.equals(v.name);
	}

	public boolean equals(Object obj) 
	{
		if(obj == this)
			return true;
		if(!(obj instanceof Variable))
			return false;
		
		Variable v = (Variable)obj;
		return 
				name.equals(v.name)
			 && exponent == v.exponent
			 && multiple == v.multiple;
	};
	
	public int hashCode() 
	{
		return name.hashCode();
	}
	
	public Variable copy()
	{
		Variable v = new Variable(name);
		v.value = value;
		v.exponent = exponent;
		return v;
	}
	
	@Override
	public String toString()
	{
		if(multiple == 0)
			return "0";
		
		if(exponent == 0)
			return "1";
		
		DecimalFormat format = new DecimalFormat("0.##");
		
		StringBuffer sb = new StringBuffer();
		if(multiple != 1)
			sb.append(format.format(multiple));
		sb.append(name);
		if(exponent != 1)
		{
			sb.append("^");
			sb.append(format.format(exponent));
		}
 
		return sb.toString();
	}
	
	
	//// Operations ////

	public Variable multiply(Variable v)
	{
		if(sameBase(v))
		{
			Variable product = copy();
			product.exponent += v.exponent;
			product.multiple *= v.multiple;
			return product;
		}
		else if(v instanceof Sum)
		{
			return ((Sum)v).multiply(this);
		}
		else
		{
			return new Product(this, v);
		}
	}

	public Variable add(Variable v)
	{
		if(equals(ZERO))
			return v.copy();
		else if(v.equals(ZERO))
			return copy();
		
//		if(v instanceof Sum)
//			return ((Sum)v).add(this);
		
		if(sameBase(v) && exponent == v.exponent)
		{
			Variable sum = copy();
			sum.multiple += v.multiple;
			return sum;
		}
		else
		{
			return new Sum(this, v);
		}
	}
	
	
	//// Accessors ////

	public Double getValue()
	{
		return value;
	}

	public void setValue(Double value)
	{
		this.value = value;
	}

	public String getName()
	{
		return name;
	}

	public void setName(String name)
	{
		this.name = name;
	}

	public double getExponent()
	{
		return exponent;
	}

	public void setExponent(double exponent)
	{
		this.exponent = exponent;
	}

	public double getMultiple()
	{
		return multiple;
	}

	public void setMultiple(double multiple)
	{
		this.multiple = multiple;
	}
}