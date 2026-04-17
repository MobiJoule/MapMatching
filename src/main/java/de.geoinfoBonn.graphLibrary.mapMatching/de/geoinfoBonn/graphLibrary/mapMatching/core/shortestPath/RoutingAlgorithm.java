package de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath;

import java.util.List;

import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;

/**
 * This interface needs to be implemented by every routing algorithm in order
 * for DijkstraRank to work.
 * 
 * @author Axel Forsch
 *
 * @param <V>
 * @param <E>
 */
public interface RoutingAlgorithm<V, E> {

	public double run(int sourceId, int targetId);

	public double run(int sourceId, int targetId, NodeVisitor<DiGraphNode<V, E>> visitor);

	public List<Integer> getPath(int targetId);

	public void setCostComputer(CostComputer<V, E> dc);

	public interface NodeVisitor<V> {
		boolean visit(V node);
	}
}
