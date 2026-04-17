package de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath;

import java.util.List;

import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.data.WeightedArcData;

/* a node visitor implementation that returns false if a given dijkstra run has reached all the nodes in a specified list 
 */
public class MultiTargetNodeVisitor<V, E extends WeightedArcData> implements Dijkstra.NodeVisitor<DiGraphNode<V, E>> {

	public Dijkstra<V, E> myDijkstra;
	public List<DiGraphNode<V, E>> targets;

	public MultiTargetNodeVisitor(List<DiGraphNode<V, E>> targets, Dijkstra<V, E> myDijkstra) {
		this.myDijkstra = myDijkstra;
		this.targets = targets;
	}

	@Override
	public boolean visit(DiGraphNode<V, E> node) {
		while (!targets.isEmpty()
				&& myDijkstra.getCost(targets.get(targets.size() - 1)) <= myDijkstra.getCost(node)) {
			targets.remove(targets.size() - 1);
		}
		return !targets.isEmpty();
	}
}
