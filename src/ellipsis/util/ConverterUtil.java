package ellipsis.util;

import org.apache.commons.math3.complex.Complex;

public class ConverterUtil
{
    public static Complex toComplex(Object object)
    {
        if(object instanceof Complex)
            return (Complex)object;
        
        if(object instanceof Number)
            return new Complex(((Number)object).doubleValue());
        
        return new Complex(Double.valueOf(object.toString()));
    }
}
