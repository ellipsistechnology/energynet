package ellipsis.util;

import java.util.Collection;
import java.util.HashSet;


public class ListUtil
{
    public static interface Setter<T, V>
    {
        void set(T target, V value);
    }
    
    public static interface Getter<T, R>
    {
        R get(T target);
    }
    
    public static <T, V> void setEach(Collection<T> c, Setter<T, V> setter, V value)
    {
        for (T t : c)
        {
            setter.set(t, value);
        }
    }
    
    public static <T, R> Collection<R> getEach(Collection<T> c, Getter<T, R> getter)
    {
        Collection<R> r = new HashSet<R>();
        for (T t : c)
        {
            r.add(getter.get(t));
        }
        return r;
    }
    
    public static <T> double sumEach(Collection<T> c, Getter<T, Double> getter)
    {
        double sum = 0;
        for (T t : c)
        {
            sum += getter.get(t);
        }
        return sum;
    }
}
