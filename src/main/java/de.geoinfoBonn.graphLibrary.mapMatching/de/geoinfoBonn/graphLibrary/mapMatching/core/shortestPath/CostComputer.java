package de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath;

/**
 * Interface of a DistanceComputer which determines how the distance is computed
 * 
 * @author
 *
 * @param <V> NodeData
 * @param <E> ArcData
 */
public interface CostComputer<V, E> {
	public double computeDistance(int uId, V u, int vId, V v, E e);
}
