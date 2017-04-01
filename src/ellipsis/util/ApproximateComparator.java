package ellipsis.util;

import static java.lang.Math.abs;

public class ApproximateComparator
{
    /**
     * @return a - b < precision
     */
    public static boolean approxLessThan(double a, double b, double precision)
    {
        double difference = a - b;
        return difference < precision;
    }

    /**
     * @return a - b > -precision
     */
    public static boolean approxGreaterThan(double a, double b, double precision)
    {
        double difference = a - b;
        return difference > -precision;
    }


    /**
     * @return |a - b| < precision
     */
    public static boolean approxEquals(double a, double b, double precision)
    {
        double difference = a - b;
        return abs(difference) < precision;
    }
}