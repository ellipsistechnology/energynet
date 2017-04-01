package ellipsis.util;

import java.util.HashMap;
import java.util.Map;

public class TreeCache<T> 
{
	private static class Node<T>
	{
		private T value;
		private Map<Object, Node<T>> children = new HashMap<Object, Node<T>>();
	}
	
	public static interface LazyCallback<T> 
	{ 
		T value(); 
	}
	
	private Node<T> root = new Node<T>();
	
	public T get(Object... keys)
	{
		Node<T> current = root;
		for (int i = 0; i < keys.length; i++) 
		{
			if(current == null)
				return null;
			current = current.children.get(keys[i]);
		}
		if(current == null)
			return null;
		synchronized (current) 
		{
			return current.value;
		}
	}
	
	public void set(T value, Object... keys)
	{
		Node<T> current = lazyGetNode(keys);
		synchronized (current) 
		{
			current.value = value;
		}
	}
	
	public T lazyGet(LazyCallback<T> callback, Object... keys)
	{
		Node<T> current = lazyGetNode(keys);
		
		T value;
		synchronized (current) 
		{
			value = current.value;
		}
		if(value == null)
		{
			value = callback.value();
			synchronized (current) 
			{
				current.value = value;
			}
		}
		
		synchronized (current) 
		{
			return current.value;
		}
	}

	private Node<T> lazyGetNode(Object... keys) 
	{
		Node<T> current = root;
		for (int i = 0; i < keys.length; i++) 
		{
			Node<T> node = current.children.get(keys[i]);
			if(node == null)
			{
				synchronized (current) 
				{
					node = new Node<T>();
					current.children.put(keys[i], node);
				}
			}
			current = node;
		}
		return current;
	}
	
	public static void main(String[] args) 
	{
		Object[] keys1 = new String[]{"A", "B", "C", "D"};
		String value1 = "V1";
		
		Object[] keys2 = new String[]{"A", "B", "C", "E"};
		String value2 = "V2";
		
		Object[] keys3 = new String[]{"A", "B", "C", "F"};
		final String value3 = "V3";
		
		Object[] keys4 = new String[]{"X", "B", "C", "D"};
		final String value4 = "V4";
		
		Object[] keys5 = new String[]{"A", "X", "C", "D"};
		final String value5 = "V5";
		
		Object[] keys6 = new String[]{"A", "B", "X", "D"};
		String value6 = "V6";
		
		TreeCache<String> tc = new TreeCache<String>();
		
		tc.set(value1, keys1);
		tc.set(value2, keys2);
		tc.set(value6, keys6);
		
		assert tc.get(keys1) == value1;
		assert tc.get(keys2) == value2;
		assert tc.get(keys6) == value6;
		
		assert tc.lazyGet(new LazyCallback<String>() {
			
			@Override
			public String value() {
				return value3;
			}
		}, keys3) == value3;
		
		assert tc.lazyGet(new LazyCallback<String>() {
			
			@Override
			public String value() {
				assert false : "Called callback even though value should have already been cached";
				return value3;
			}
		}, keys3) == value3 : "lazyGet returned incorrect value";

		assert tc.lazyGet(new LazyCallback<String>() {
			
			@Override
			public String value() {
				return value4;
			}
		}, keys4) == value4;
		
		assert tc.lazyGet(new LazyCallback<String>() {
			
			@Override
			public String value() {
				return value5;
			}
		}, keys5) == value5;
		
		assert tc.get(keys3) == value3 : "get returned the incorrect value for value 3 as set by lazyGet";
		assert tc.get(keys5) == value5 : "get returned the incorrect value for value 5 as set by lazyGet";
	}
}