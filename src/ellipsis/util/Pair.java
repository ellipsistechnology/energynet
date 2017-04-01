package ellipsis.util;

import java.io.Serializable;

public class Pair<K,V> implements Serializable
{
	private static final long serialVersionUID = -6236948249874760563L;
	private K k;
	private V v;

	public Pair(K k, V v)
	{
		this.setKey(k);
		this.setValue(v);
	}

	@Override
	public String toString()
	{
		return "<"+getKey()+","+getValue()+">";
	}

	public K getKey()
	{
		return k;
	}

	public void setKey(K k)
	{
		this.k = k;
	}

	public V getValue()
	{
		return v;
	}

	public void setValue(V v)
	{
		this.v = v;
	}
}