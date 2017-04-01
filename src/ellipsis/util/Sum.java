package ellipsis.util;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;


public class Sum 
{
    //// Lambda expression interfaces ////
    
    public static interface IndexedFunction
    {
        double value(int i);
    }
    
    public static interface NumberFunction
    {
        double value(double d);
    }
    
    public static interface IndexedBooleanFunction
    {
        boolean value(int i);
    }
    
    public static interface ValueFunction<T>
    {
        double value(T t);
    }
    
    public static interface IteratorVectorFunction<T>
    {
        RealVector value(T t);
    }
    
    
    //// Sum Functions ////
    
    public static <T> double sum(ValueFunction<T> f, Iterable<T> it)
    {
        double sum = 0;
        for (T t : it)
        {
            sum += f.value(t);
        }
        return sum;
    }
    
    public static <T> RealVector sumV(IteratorVectorFunction<T> f, Iterable<T> it, int dimension)
    {
        RealVector sum = new ArrayRealVector(dimension);
        for (T t : it)
        {
            sum = sum.add(f.value(t));
        }
        return sum;
    }
    
    public static double sum(IndexedFunction f, int dimension, int exclude)
    {
        double sum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            if(i != exclude)
            {
                sum += f.value(i);
            }
        }
        
        return sum;
    }
    
    public static double sum(IndexedFunction f, int dimension)
    {
        double sum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            sum += f.value(i);
        }
        
        return sum;
    }
    
    public static double sum(double[] ds, NumberFunction f)
    {
        double sum = 0;
        for(int i = 0; i < ds.length; ++i)
        {
            sum += f.value(ds[i]);
        }
        
        return sum;
    }
    
    
    //// Old Stuff ////
    
    @Deprecated
	public static interface Adder
	{
		double add(int i);
	}
	
    @Deprecated
	public static double sum(int i_0, int N, Adder adder)
	{
		double sum = 0;
		for(int i = i_0; i <= N; ++i)
		{
			sum += adder.add(i);
		}
		return sum;
	}
}
