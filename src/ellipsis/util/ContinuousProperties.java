package ellipsis.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Properties;

/**
 * Given a key that does not exist in the properties an approximation based on a linear extrapolation
 * will be returned from {@link #get(Object)}.
 * All keys and values must be doubles.
 * After all elements are added or loaded into the properties {@link #sort()} must be called, otherwise
 * results will be unpredictable.
 * @author bmillar
 *
 */
public class ContinuousProperties extends Properties
{
    private static final long serialVersionUID = -583217139998400039L;
    
    private Object keys[];
    
    public void sort()
    {
        keys = new Object[size()];
        int i = 0;
        for (Object key : keySet())
        {
            keys[i++] = key;
        }
        Arrays.sort(keys);
    }
    
    /**
     * Key and value will be converted to Doubles.
     */
    @Override
    public synchronized Object put(Object key, Object value)
    {
        key = toDouble(key);
        value = toDouble(value);
        
        return super.put(key, value);
    }

    private Double toDouble(Object d)
    {
        if(d instanceof Double)
            return (Double) d;
        if(d instanceof Number)
            return ((Number)d).doubleValue();
        else
            return Double.valueOf(d.toString());
    }

    @Override
    public synchronized Object get(Object key)
    {
        double x = toDouble(key);
        
        // If the exact value exists then return it:
        Object value = super.get(x);
        if(value != null)
            return value;
        
        // Otherwise, find the nearest match:
        int i = Arrays.binarySearch(keys, x);
        if(i >= 0) // Theoretically impossible since an exact match should be returned above
            throw new RuntimeException("Key did not exist in properties, yet search returned non-negative index.");
        
        int insertionPoint = -(i+1); // Ref. http://docs.oracle.com/javase/1.4.2/docs/api/java/util/Arrays.html#binarySearch(byte[], byte)
        
        // Can't extrapolate for first or last elements, 
        // so return first or last respectively:
        if(insertionPoint == 0)
            return super.get(keys[0]);
        if(insertionPoint == size())
            return super.get(keys[keys.length-1]);
        
        // Get upper and lower values:
        int upperIndex = insertionPoint;
        int lowerIndex = insertionPoint - 1;
        
        Object upperKey = keys[upperIndex];
        Object lowerKey = keys[lowerIndex];
        Object upperValue = super.get(upperKey);
        Object lowerValue = super.get(lowerKey);
        
        double x1 = (Double)lowerKey;
        double x2 = (Double)upperKey;
        double y1 = (Double)lowerValue;
        double y2 = (Double)upperValue;
        
        // Extrapolate:
        double m = (y2-y1)/(x2-x1); // slope
        double y = m*(x-x1) + y1; // linear relationship
//System.out.println("Upper:keys["+upperIndex+"]="+upperKey+"=>"+upperValue+"Lower:keys["+lowerIndex+"]="+lowerKey+"=>"+lowerValue);
        return y;
    }

    public double getDouble(String key)
    {
        return (Double)get(key);
    }
    
    // Test
    public static void main(String[] args) throws IOException
    {
        FileInputStream in = new FileInputStream(new File("/svn/energynet/casestudies/src/ellipsis/energy/casestudies/ieee/wind24.properties"));
        ContinuousProperties p = new ContinuousProperties();
        p.load(in);
        p.sort();
        for (int i = 0; i < 2400;) // 15 minute intervals
        {
            String key = ""+i;
            System.out.println(key+"="+p.get(key));
            if(i%100 == 45)
                i += 55; // Go to next hour
            else
                i += 15; // Next 15 minutes
        }
    }
}
