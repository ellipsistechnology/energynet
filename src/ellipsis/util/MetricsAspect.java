package ellipsis.util;

import java.lang.reflect.Method;
import java.util.Map;

import com.mls.util.AspectInterceptor;
import com.mls.util.AspectInterceptor.Aspected;
import com.mls.util.Timer;

/**
 * 
 * @author bmillar
 *
 */
public class MetricsAspect implements Aspected
{
	private Map<Object, Object> results;
	
	public static <T> T newWrappingInterceptor(Class<T> type, T wrapped, Map<Object, Object> results)
	{
		return AspectInterceptor.newWrappingInterceptor(type, new MetricsAspect(results), wrapped);
	}
	
	public static <T> T newInterceptor(Class<T> type, Map<Object,Object> results)
	{
		return AspectInterceptor.newInterceptor(type, new MetricsAspect(results));
	}
	
	private MetricsAspect(Map<Object, Object> results)
	{
		this.results = results;
	}

	@Override
	public void preProcess(Object obj, Method method, Object[] args)
	{
		String id = id(obj.getClass().getSuperclass(), method);
		Timer.getGlobalTimer(id).start();
		Object count = results.get(id);
		if(count != null)
			results.put(id, ((int)count)+1);
		else
			results.put(id, 1);
	}

	public static String id(Class<?> c, Method method)
	{
		return c.getName()+"."+method.getName();
	}

	@Override
	public void postProcess(Object obj, Method method, Object[] args, Object result)
	{
		Timer.getGlobalTimer(id(obj.getClass().getSuperclass(), method)).stop();
	}

	@Override
	public void handleThrowable(Object obj, Method method, Object[] args, Throwable t){}

	@Override
	public void finalProcess(Object obj, Method method, Object[] args){}
}