package de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.data.WeightedArcData;
import de.geoinfoBonn.graphLibrary.mapMatching.core.structures.MinHeap;
import de.geoinfoBonn.graphLibrary.mapMatching.core.structures.MinHeap.HeapItem;

/**
 * Simple implementation of the Dijkstra algorithm. The data structures are
 * initialized. Due to some internal time stamp mechanism, the class supports
 * multiple executions without the need to reinitialize the data structures.
 */
public class Dijkstra<V, E> implements RoutingAlgorithm<V, E> {

	protected double distToSource = 0;
	protected double cost[];
	protected int stamps[];
	protected HeapItem<DiGraphNode<V, E>> items[];
	protected DiGraphNode<V, E> pred[];
	protected DiGraph<V, E> graph;

	protected double currentDist = 0;
	protected int currentStamp = 0;

	protected CostComputer<V, E> DEFAULT_EDGE_COST = new CostComputer<V, E>() {

		@Override
		public double computeDistance(int uId, V u, int vId, V v, E e) {
			if (e instanceof WeightedArcData)
				return ((WeightedArcData) e).getValue();
			throw new IllegalArgumentException(
					"Edge weights for default-edge-distance must implement WeightedArcData!");
		}
	};

	protected CostComputer<V, E> dc = DEFAULT_EDGE_COST;

	@SuppressWarnings("unchecked")
	public Dijkstra(DiGraph<V, E> g) {
		this.cost = new double[g.n()];
		this.stamps = new int[g.n()];
		this.items = new HeapItem[g.n()];
		this.pred = new DiGraphNode[g.n()];
		this.graph = g;
	}

	@SuppressWarnings("unchecked")
	public Dijkstra(DiGraph<V, E> g, CostComputer<V, E> dc) {
		this.cost = new double[g.n()];
		this.stamps = new int[g.n()];
		this.items = new HeapItem[g.n()];
		this.pred = new DiGraphNode[g.n()];
		this.dc = dc;
		this.graph = g;
	}

	/**
	 * Runs the algorithm starting at the source node. When the target is reached,
	 * the search is aborted and the cost to the target is returned. In case
	 * that the target is not reachable from the source, the method returns
	 * Double.MAX_VALUE.
	 */
	public double run(DiGraphNode<V, E> source, DiGraphNode<V, E> target) {
		run(source, new NodeVisitor<DiGraphNode<V, E>>() {
			@Override
			public boolean visit(DiGraphNode<V, E> node) {
				return node != target;
			}
		});
		if (stamps[target.getId()] == currentStamp) {
			return cost[target.getId()];
		}
		return Double.MAX_VALUE;
	}

	public double run(DiGraphNode<V, E> source, DiGraphNode<V, E> target, NodeVisitor<DiGraphNode<V, E>> visitor) {
		run(source, new NodeVisitor<DiGraphNode<V, E>>() {
			@Override
			public boolean visit(DiGraphNode<V, E> node) {
				visitor.visit(node);
				return node != target;
			}
		});
		if (stamps[target.getId()] == currentStamp) {
			return cost[target.getId()];
		}
		return Double.MAX_VALUE;
	}

	/**
	 * Runs the algorithm starting at the source node. When the target is reached,
	 * the search is aborted and the distance to the target is returned. In case
	 * that the target is not reachable from the source, the method returns
	 * Double.MAX_VALUE.
	 */
	@Override
	public double run(int sourceId, int targetId) {
		DiGraphNode<V, E> source = graph.getNode(sourceId);
		DiGraphNode<V, E> target = graph.getNode(targetId);

		return run(source, target);
	}

	@Override
	public double run(int sourceId, int targetId, NodeVisitor<DiGraphNode<V, E>> visitor) {
		DiGraphNode<V, E> source = graph.getNode(sourceId);
		DiGraphNode<V, E> target = graph.getNode(targetId);

		return run(source, target, visitor);
	}

	public void run(DiGraphNode<V, E> source) {

		run(source, new NodeVisitor<DiGraphNode<V, E>>() {
			@Override
			public boolean visit(DiGraphNode<V, E> node) {
				return true;
			}
		});
	}

	public boolean run(DiGraphNode<V, E> start, NodeVisitor<DiGraphNode<V, E>> visitor) {
		HeapItem<DiGraphNode<V, E>> settled = new HeapItem<>();

		currentStamp++;
		cost[start.getId()] = distToSource;
		pred[start.getId()] = null;

		MinHeap<DiGraphNode<V, E>> queue = new MinHeap<DiGraphNode<V, E>>();

		items[start.getId()] = queue.insertItem(distToSource, start);
		stamps[start.getId()] = currentStamp;

		while (queue.size() > 0) {
			HeapItem<DiGraphNode<V, E>> item = queue.extractMin();
			DiGraphNode<V, E> u = item.getValue();
			currentDist = item.getKey();

			items[u.getId()] = settled;
			if (!visitor.visit(u)) {
				return false;
			}
			for (DiGraphArc<V, E> edge : u.getOutgoingArcs()) {
				DiGraphNode<V, E> v = edge.getTarget();

				double alt = cost[u.getId()]
						+ dc.computeDistance(u.getId(), u.getNodeData(), v.getId(), v.getNodeData(), edge.getArcData());
				if (items[v.getId()] != settled && (stamps[v.getId()] < currentStamp || alt < cost[v.getId()])) {
					cost[v.getId()] = alt;
					pred[v.getId()] = u;
					if (stamps[v.getId()] < currentStamp || items[v.getId()] == null) {
						items[v.getId()] = queue.insertItem(alt, v);
					} else {
						queue.decreaseKey(items[v.getId()], alt);
					}
					stamps[v.getId()] = currentStamp;
				}
			}
		}
		return true;
	}

	/**
	 * Assumes that the Dijkstra algorithm has been executed before. Returns the
	 * path of node to the target (stored in an ArrayList from start to target).
	 */
	public List<DiGraphNode<V, E>> getPath(DiGraphNode<V, E> target) {
		LinkedList<DiGraphNode<V, E>> path = new LinkedList<DiGraphNode<V, E>>();

		if (stamps[target.getId()] < currentStamp) {
			return new ArrayList<DiGraphNode<V, E>>();
		}

		DiGraphNode<V, E> current = target;
		while (current != null) {
			path.addFirst(current);
			current = pred[current.getId()];
		}

		return new ArrayList<DiGraphNode<V, E>>(path);
	}

	/**
	 * Assumes that the Dijkstra algorithm has been executed before. Returns the
	 * path of node to the target (stored in an ArrayList from start to target).
	 */
	public List<DiGraphArc<V, E>> getPathArcs(DiGraphNode<V, E> target) {
		LinkedList<DiGraphArc<V, E>> path = new LinkedList<>();

		if (stamps[target.getId()] < currentStamp) {
			return new ArrayList<>();
		}

		DiGraphNode<V, E> current = target;
		DiGraphNode<V, E> next = null;
		while (current != null) {
			if (next != null) {
				path.addFirst(current.getFirstOutgoingArcTo(next));
			}
			next = current;
			current = pred[current.getId()];
		}

		return new ArrayList<>(path);
	}

	@Override
	public List<Integer> getPath(int targetId) {
		LinkedList<Integer> path = new LinkedList<>();

		if (stamps[targetId] < currentStamp) {
			return new ArrayList<>();
		}

		DiGraphNode<V, E> current = graph.getNode(targetId);
		while (current != null) {
			path.addFirst(current.getId());
			current = pred[current.getId()];
		}

		return new ArrayList<>(path);
	}

	/**
	 * Assumes that the dijkstra has been executed before. Returns the distance of
	 * <code>node</code> to the start node of the last run.
	 * 
	 * @param node
	 * @return distance of the node to the start node of the last run. If the node
	 *         is not reachable from the start node, then it returns
	 *         Double.MAX_VALUE.
	 */
	public double getCost(DiGraphNode<V, E> node) {
		return stamps[node.getId()] < currentStamp ? Double.MAX_VALUE : cost[node.getId()];
	}

	public DiGraphNode<V, E> getPredecessor(DiGraphNode<V, E> node) {
		if (stamps[node.getId()] < currentStamp) {
			return null;
		}
		return pred[node.getId()];
	}

	@Override
	public void setCostComputer(CostComputer<V, E> dc) {
		this.dc = dc;
	}

	public CostComputer<V, E> getCostComputer() {
		return dc;
	}

	public static interface NodeIterator<V, E> {
		Iterator<DiGraphNode<V, E>> getIterator(DiGraphNode<V, E> s);

		double getWeightOfCurrentArc(DiGraphNode<V, E> s, DiGraphNode<V, E> t);
	}

	/**
	 * NodeIterator that depends mainly on the graph structure. Additional
	 * adjacencies may be added via an additional NodeIterator.
	 */
	public static class BasicAdjacentNodeIterator<V, E> implements NodeIterator<V, E> {

		protected DiGraphArc<V, E> currArc;
		protected NodeIterator<V, E> addIt;
		protected boolean start;
		protected CostComputer<V, E> dc;

		public BasicAdjacentNodeIterator(NodeIterator<V, E> additionalIterator, CostComputer<V, E> dc) {
			this.addIt = additionalIterator;
			this.start = true;
			this.currArc = null;
			this.dc = dc;
		}

		@Override
		public double getWeightOfCurrentArc(DiGraphNode<V, E> s, DiGraphNode<V, E> t) {
			if (currArc != null) {
				return dc.computeDistance(currArc.getSource().getId(), currArc.getSource().getNodeData(),
						currArc.getTarget().getId(), currArc.getTarget().getNodeData(), currArc.getArcData());
			}
			if (start) {
				return Double.NaN;
			}
			return addIt.getWeightOfCurrentArc(s, t);
		}

		@Override
		public Iterator<DiGraphNode<V, E>> getIterator(DiGraphNode<V, E> s) {
			Iterator<DiGraphArc<V, E>> it = s.getOutgoingArcs().iterator();
			Iterator<DiGraphNode<V, E>> extraIt;
			if (addIt == null) {
				extraIt = null;
			} else {
				extraIt = addIt.getIterator(s);
			}
			return new Iterator<DiGraphNode<V, E>>() {

				@Override
				public boolean hasNext() {
					if (it.hasNext()) {
						return true;
					}
					if (extraIt == null) {
						return false;
					}
					return extraIt.hasNext();
				}

				@Override
				public DiGraphNode<V, E> next() {
					start = false;
					if (it.hasNext()) {
						currArc = it.next();
						return currArc.getTarget();
					}
					if (extraIt == null) {
						return null;
					}
					currArc = null;
					return extraIt.next();
				}
			};
		}
	}

	/**
	 * NodeIterator that depends only on the graph structure.
	 */
	protected NodeIterator<V, E> DEFAULT_ADJACENT_NODE_ITERATOR = new BasicAdjacentNodeIterator<V, E>(null, this.dc);

	/**
	 * Runs the algorithm starting at the source node using the given NodeVisitor.
	 * Runs as long as the NodeVisitor's method visit returns true. Ignores the
	 * graph structure; adjacencies depend on the NodeIterator.
	 * 
	 * @throws OutOfStreetNetworkException
	 */
	public boolean run(DiGraphNode<V, E> source, NodeVisitor<DiGraphNode<V, E>> visitor, NodeIterator<V, E> nit) {
		currentStamp++;
		cost[source.getId()] = distToSource;
		pred[source.getId()] = null;
		double weightOfArc;

		MinHeap<DiGraphNode<V, E>> queue = new MinHeap<DiGraphNode<V, E>>();

		items[source.getId()] = queue.insertItem(distToSource, source);
		stamps[source.getId()] = currentStamp;

		while (queue.size() > 0) {
			HeapItem<DiGraphNode<V, E>> item = queue.extractMin();
			DiGraphNode<V, E> u = item.getValue();
			currentDist = item.getKey();

			if (!visitor.visit(u)) {
				return false;
			}

			for (Iterator<DiGraphNode<V, E>> it = nit.getIterator(u); it.hasNext();) {
				DiGraphNode<V, E> v = it.next();

				weightOfArc = nit.getWeightOfCurrentArc(u, v);
				discoverNode(u, v, queue, cost[u.getId()] + weightOfArc);
			}
		}
		return true;
	}

	private void discoverNode(DiGraphNode<V, E> curr, DiGraphNode<V, E> target, MinHeap<DiGraphNode<V, E>> queue,
			double alt) {
		if (stamps[target.getId()] < currentStamp || alt < cost[target.getId()]) {
			cost[target.getId()] = alt;
			pred[target.getId()] = curr;
			if (stamps[target.getId()] < currentStamp || items[target.getId()] == null) {
				items[target.getId()] = queue.insertItem(alt, target);
			} else {
				queue.decreaseKey(items[target.getId()], alt);
			}
			stamps[target.getId()] = currentStamp;
		}
	}

	public double getCurrentDist() {
		return currentDist;
	}

	public int[] getStamps() {
		int[] copyOfStamps = new int[stamps.length];
		for (int i = 0; i < stamps.length; ++i)
			copyOfStamps[i] = stamps[i];
		return copyOfStamps;
	}

	public int getCurrentStamp() {
		return currentStamp;
	}

	public DiGraphNode<V, E>[] getPred() {
		@SuppressWarnings("unchecked")
		DiGraphNode<V, E>[] copyOfPred = new DiGraphNode[pred.length];
		for (int i = 0; i < pred.length; ++i)
			copyOfPred[i] = pred[i];
		return copyOfPred;
	}

	public double getDistToSource() {
		return distToSource;
	}

	public void setDistToSource(double distToSource) {
		this.distToSource = distToSource;
	}

}