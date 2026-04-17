package de.geoinfoBonn.graphLibrary.mapMatching.core.generic;

import java.awt.geom.Point2D;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import de.geoinfoBonn.graphLibrary.mapMatching.core.dcel.Dcel;

/**
 * Implementation of a graph network for computing different algorithms
 * 
 * @author
 *
 * @param <V> - NodeData
 * @param <E> - ArcData
 */
public class DiGraph<V, E> implements VisitableGraph<V, E>, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Implementation of a GraphNode corresponding to DiGraph
	 * 
	 * @author
	 *
	 * @param <V> - NodeData
	 * @param <E> - ArcData
	 */
	public static class DiGraphNode<V, E> implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * list of every outgoing arc
		 */
		private List<DiGraphArc<V, E>> outgoingArcs;

		/**
		 * list of every incoming arc
		 */
		private List<DiGraphArc<V, E>> incomingArcs;
		/**
		 * NodeData of the DiGraphNode
		 * 
		 */
		private V nodeData;

		/**
		 * id increases by every addition of a node to the DiGraph
		 */
		private int id;

		/**
		 * Constructor
		 * 
		 * @param nodeData -
		 * @param id       -
		 */
		private DiGraphNode(V nodeData, int id) {
			outgoingArcs = new ArrayList<DiGraphArc<V, E>>();
			incomingArcs = new ArrayList<DiGraphArc<V, E>>();
			this.nodeData = nodeData;
			this.id = id;
		}

		/**
		 * Getter-Method for the outgoing arcs
		 * 
		 * @return list of outgoing arcs
		 */
		public List<DiGraphArc<V, E>> getOutgoingArcs() {
			return outgoingArcs;
		}

		/**
		 * Getter-Method for the incoming arcs
		 * 
		 * @return list of incoming arcs
		 */
		public List<DiGraphArc<V, E>> getIncomingArcs() {
			return incomingArcs;
		}

		/**
		 * Getter-Method for the id of a node
		 * 
		 * @return id of the node
		 */
		public int getId() {
			return id;
		}

		/**
		 * Setter-Method for the id
		 * 
		 */
		public void setId(int id) {
			this.id = id;
		}

		/**
		 * Getter-Method for the NodeData
		 * 
		 * @return NodeData V
		 */
		public V getNodeData() {
			return nodeData;
		}

		/**
		 * toString-method that returns the toString-method of the NodeData V
		 */
		@Override
		public String toString() {
			return "[" + id + "," + nodeData.toString() + "]";
		}

		/**
		 * outDegree is the same as the number of outgoing arcs which is a
		 * characteristic value for a DiGraph(Node)
		 * 
		 * @return outDegree of the node
		 */
		public int outDegree() {
			return outgoingArcs.size();
		}

		/**
		 * inDegree is the same as the number of incoming arcs which is a characteristic
		 * value for a DiGraph(Node)
		 * 
		 * @return inDegree of the node
		 */
		public int inDegree() {
			return incomingArcs.size();
		}

		/**
		 * degree gives the number of unique adjacent nodes
		 * 
		 * @return number of unique adjacent nodes
		 */
		public int degree() {
			Set<DiGraphNode<V, E>> incidentNodes = new HashSet<>();
			for (DiGraphArc<V, E> inc : incomingArcs)
				incidentNodes.add(inc.source);
			for (DiGraphArc<V, E> out : outgoingArcs)
				incidentNodes.add(out.target);
			return incidentNodes.size();
		}

		/**
		 * 
		 * @param target - DiGraphNode
		 * @return the first arc in outgoing arcs which has the targetNode as its target
		 */
		public DiGraphArc<V, E> getFirstOutgoingArcTo(DiGraphNode<V, E> target) {
			for (DiGraphArc<V, E> a : outgoingArcs) {
				if (a.getTarget() == target) {
					return a;
				}
			}
			return null;
		}

		/**
		 * 
		 * @param source - DiGraphNode
		 * @return the first arc in incoming arcs which has the sourceNode as its source
		 */
		public DiGraphArc<V, E> getFirstIncomingArcFrom(DiGraphNode<V, E> source) {
			for (DiGraphArc<V, E> a : incomingArcs) {
				if (a.getSource() == source) {
					return a;
				}
			}
			return null;
		}
	}

	/**
	 * Implementation of a GraphArc corresponding to DiGraph
	 * 
	 * @author
	 *
	 * @param <V> - NodeData
	 * @param <E> - ArcData
	 */
	public static class DiGraphArc<V, E> implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * source node of the arc
		 */
		private DiGraphNode<V, E> source;

		/**
		 * target node of the arc
		 */
		private DiGraphNode<V, E> target;

		/**
		 * for adding extra information on the arc (e.g. length)
		 */
		private E arcData;
		/**
		 * id increases by every addition of an arc to the DiGraph
		 */
		private int id;

		/**
		 * Constructor
		 * 
		 * @param source  - DiGraphNode of the DiGraph
		 * @param target  - DiGraphNode of the DiGraph
		 * @param arcData -
		 * @param id
		 */
		protected DiGraphArc(DiGraphNode<V, E> source, DiGraphNode<V, E> target, E arcData, int id) {
			this.source = source;
			this.target = target;
			this.arcData = arcData;
			this.id = id;
		}

		/**
		 * Getter-method of the source node
		 * 
		 * @return source node of the arc
		 */
		public DiGraphNode<V, E> getSource() {
			return source;
		}

		/**
		 * Getter-method of the target node
		 * 
		 * @return target node of the arc
		 */
		public DiGraphNode<V, E> getTarget() {
			return target;
		}

		/**
		 * Getter-method for the arcData
		 * 
		 * @return arcData of the arc
		 */
		public E getArcData() {
			return arcData;
		}

		/**
		 * Setter-method of the id of the arc
		 * 
		 * @param id - arc changes id to this value
		 */
		public void setId(int id) {
			this.id = id;
		}

		/**
		 * Getter-method of the id
		 * 
		 * @return id of the arc
		 */
		public int getId() {
			return id;
		}

		/**
		 * Finds the corresponding twin of the arc a twin has the arcs target as its
		 * source and the arcs source as its target
		 * 
		 * @return the twin of the arc. if there is no twin return null
		 * 
		 */
		public DiGraphArc<V, E> getTwin() {
			for (DiGraphArc<V, E> a : source.getIncomingArcs()) {
				if (a.getSource() == this.getTarget())
					return a;
			}
			return null;
		}

		public boolean hasTwin() {
			return getTwin() != null;
		}

		/**
		 * Setter-method of the arcData
		 * 
		 * @param arcData - data that will become the new arcData of the arc
		 */
		public void setArcData(E arcData) {
			this.arcData = arcData;
		}

		public double getInclination(boolean swapping) {
			double inc = 0;
			double Y = 0;
			double X = 0;

			if (swapping) {
				Y = ((Point2D) this.source.nodeData).getY() - ((Point2D) this.target.nodeData).getY();
				X = ((Point2D) this.source.nodeData).getX() - ((Point2D) this.target.nodeData).getX();
			} else {
				Y = ((Point2D) this.target.nodeData).getY() - ((Point2D) this.source.nodeData).getY();
				X = ((Point2D) this.target.nodeData).getX() - ((Point2D) this.source.nodeData).getX();
			}

			inc = Math.atan2(Y, X);

			if (inc >= Math.PI / 2 && inc <= Math.PI) {

				inc -= Math.PI / 2;
			} else {
				inc += 3 * Math.PI / 2;
			}

			return inc;
		}

		/**
		 * [source -- target] (arcData)
		 */
		@Override
		public String toString() {
			return id + " [" + source.toString() + " -- " + target.toString() + "] (" + arcData.toString() + ")";
		}
	}

	/**
	 * every added node is listed up here
	 */
	protected ArrayList<DiGraphNode<V, E>> nodeList;

	/**
	 * every added arc is listed up here
	 */
	protected ArrayList<DiGraphArc<V, E>> arcList;

	/**
	 * Default-constructor Generates empty graph
	 */
	public DiGraph() {
		nodeList = new ArrayList<DiGraphNode<V, E>>();
		arcList = new ArrayList<DiGraphArc<V, E>>();
	}

	/**
	 * Constructor Generates graph from adjacency matrix
	 * 
	 * @param nodes   - array with node data
	 * @param arcs    - adjacency matrix
	 * @param arcData - array with arc data
	 */
	public DiGraph(V[] nodes, boolean[][] arcs, E[][] arcData) {
		nodeList = new ArrayList<DiGraphNode<V, E>>();
		arcList = new ArrayList<DiGraphArc<V, E>>();
		for (V v : nodes) {
			addNode(v);
		}
		for (int i = 0; i < nodeList.size(); i++) {
			for (int j = 0; j < nodeList.size(); j++) {
				if (arcs[i][j]) {
					addArc(nodeList.get(i), nodeList.get(j), arcData[i][j]);
				}
			}
		}
	}

	/**
	 * adds a DiGraphNode to the DiGraph and updates all important lists etc.
	 * 
	 * @param nodeInfo - NodeData of the node
	 * @return added DiGraphNode
	 */
	public DiGraphNode<V, E> addNode(V nodeInfo) {
		DiGraphNode<V, E> v = new DiGraphNode<V, E>(nodeInfo, nodeList.size());
		nodeList.add(v);
		return v;
	}

	/**
	 * adds an DiGraphArc to the DiGraph and updates all important lists etc.
	 * 
	 * @param v1       - source node
	 * @param v2       - target node
	 * @param edgeData - ArcData of the arc
	 * @return added DiGraphArc
	 */
	public DiGraphArc<V, E> addArc(DiGraphNode<V, E> v1, DiGraphNode<V, E> v2, E edgeData) {
		DiGraphArc<V, E> a = new DiGraphArc<V, E>(v1, v2, edgeData, arcList.size());
		v1.outgoingArcs.add(a);
		v2.incomingArcs.add(a);
		arcList.add(a);
		return a;
	}

	public List<DiGraphArc<V, E>> addDoubleArc(DiGraphNode<V, E> v1, DiGraphNode<V, E> v2) {
		return addDoubleArc(v1, v2, null);
	}

	public List<DiGraphArc<V, E>> addDoubleArc(DiGraphNode<V, E> v1, DiGraphNode<V, E> v2, E edgeData) {
		List<DiGraphArc<V, E>> result = new ArrayList<DiGraphArc<V, E>>();
		result.add(addArc(v1, v2, edgeData));
		result.add(addArc(v2, v1, edgeData));
		return result;
	}

	/**
	 * Returns the number of nodes of the DiGraph
	 */
	@Override
	public int n() {
		return nodeList.size();
	}

	/**
	 * Returns the number of arcs of the DiGraph
	 */
	@Override
	public int m() {
		return arcList.size();
	}

	/**
	 * Returns the number of edges of the DiGraph, not double counting bidirectional
	 * arcs!
	 * 
	 * @return
	 */
	public int mUnique() {
		double m = 0;
		for (DiGraphArc<V, E> arc : arcList) {
			if (arc.hasTwin())
				m += 0.5;
			else
				m += 1;
		}
		return (int) Math.round(m);
	}

	/**
	 * Returns the node with id index in O(1)
	 * 
	 * @param index - id of the searched node
	 * @return DiGraphNode with id index
	 */
	public DiGraphNode<V, E> getNode(int index) {
		return nodeList.get(index);
	}

	/**
	 * Getter-method for the nodelist
	 * 
	 * @return ArrayList of all DiGraphNodes of the DiGraph
	 */
	public ArrayList<DiGraphNode<V, E>> getNodes() {
		return nodeList;
	}

	/**
	 * Returns the arc with id index in O(1)
	 * 
	 * @param index - id of the searched node
	 * @return DiGraphArc with id index
	 */
	public DiGraphArc<V, E> getArc(int index) {
		return arcList.get(index);
	}

	/**
	 * Getter-method for the arclist
	 * 
	 * @return ArrayList of all DiGraphArcs of the DiGraph
	 */
	public ArrayList<DiGraphArc<V, E>> getArcs() {
		return arcList;
	}

	public List<DiGraphNode<V, E>> toNodeList(List<Integer> idList) {
		List<DiGraphNode<V, E>> list = new LinkedList<>();
		for (int id : idList) {
			list.add(this.getNode(id));
		}
		return list;
	}

	public List<DiGraphArc<V, E>> toArcList(List<Integer> idList) {
		List<DiGraphArc<V, E>> list = new LinkedList<>();
		for (int id : idList) {
			list.add(this.getArc(id));
		}
		return list;
	}

	/**
	 * Setter-method for the node list that requires an ArrayList of DiGraphNodes
	 * 
	 * @param nodes list of nodes
	 */
	public void setNodes(ArrayList<DiGraphNode<V, E>> nodes) {
		this.nodeList = nodes;
	}

	/**
	 * Setter method for the arc list that requires an ArrayList of DiGraphArcs
	 * 
	 * @param arcs list of arcs
	 */
	public void setArcs(ArrayList<DiGraphArc<V, E>> arcs) {
		this.arcList = arcs;
	}

	/**
	 * method to remove a set of arcs in O(m) time (assuming constant-time hashing)
	 * 
	 * @param arcsToBeRemoved - HashSet of DiGraphArcs that will be removed from the
	 *                        DiGraph
	 */
	public void removeArcs(Set<DiGraphArc<V, E>> arcsToBeRemoved) {
		ArrayList<DiGraphArc<V, E>> arcListNew = new ArrayList<>();

		// O(m) in any case! It would be better if you go through the
		// arcsToBeRemoved-list and then remove these arcs from the arcList
		for (DiGraphArc<V, E> a : arcList)
			if (!arcsToBeRemoved.contains(a))
				arcListNew.add(a);

		for (DiGraphArc<V, E> a : arcsToBeRemoved) {
			DiGraphNode<V, E> sourceNode = a.getSource();
			DiGraphNode<V, E> targetNode = a.getTarget();
			sourceNode.getOutgoingArcs().remove(a);
			targetNode.getIncomingArcs().remove(a);
		}

		this.arcList = arcListNew;
	}

	/**
	 * method to remove an arc in O(m) time
	 * 
	 * @param arcToBeRemoved - DiGraphArc that will be removed from the DiGraph
	 */
	public void removeArc(DiGraphArc<V, E> arcToBeRemoved) {
		DiGraphNode<V, E> sourceNode = arcToBeRemoved.getSource();
		DiGraphNode<V, E> targetNode = arcToBeRemoved.getTarget();

		sourceNode.getOutgoingArcs().remove(arcToBeRemoved);
		targetNode.getIncomingArcs().remove(arcToBeRemoved);
		arcList.remove(arcToBeRemoved);
	}

	/**
	 * method to remove a set of nodes and all their incident edges in O(n + m) time
	 * 
	 * @param nodesToBeRemoved - HashSet of DiGraphNodes that will be removed from
	 *                         the DiGraph
	 */
	public void removeNodes(Set<DiGraphNode<V, E>> nodesToBeRemoved) {
		// remove incident arcs
		Set<DiGraphArc<V, E>> arcsToBeRemoved = new HashSet<>();
		for (DiGraphNode<V, E> v : nodesToBeRemoved) {
			arcsToBeRemoved.addAll(v.incomingArcs);
			arcsToBeRemoved.addAll(v.outgoingArcs);
		}
		this.removeArcs(arcsToBeRemoved);

		ArrayList<DiGraphNode<V, E>> nodeListNew = new ArrayList<>();
		for (DiGraphNode<V, E> v : nodeList)
			if (!nodesToBeRemoved.contains(v))
				nodeListNew.add(v);

		nodeList = nodeListNew;
	}

	public void sort(Comparator<DiGraphArc<V, E>> outgoingComp, Comparator<DiGraphArc<V, E>> incomingComp) {
		for (DiGraphNode<V, E> node : nodeList) {
			node.outgoingArcs.sort(outgoingComp);
			node.incomingArcs.sort(incomingComp);
		}
	}

	/**
	 * Sort this graph using the {@link LineComparator}s necessary for generating a
	 * {@link Dcel} afterwards.
	 * 
	 * @return returns this graph to allow method chaining
	 */
	public DiGraph<V, E> sort() {
		this.sort(new LineComparator<>(false), new LineComparator<>(true));
		return this;
	}

	public DiGraphArc<V, E> getOuterArc() {
		// if(this.nodeList.get(0).getNodeData() instanceof Point2D) {
		// DiGraphNode p1 = nodeList.get(0);
		// DiGraphNode p2 = nodeList.get(0).getIncomingArcs();
		//
		// }else if(this.arcList.get(0).getArcData() instanceof Line2D){
		//
		// }else {
		// throw new RuntimeException("Cannot find outer Arc. No geometric
		// information.");
		// }
		return nodeList.get(0).getIncomingArcs().get(0);
		// return nodeList.get(0).getOutgoingArcs().get(0);
	}

	/**
	 * method to remove a node and all its incident edges in O(n + m) time
	 * 
	 * @param nodeToBeRemoved - DiGraphNode that will be removed from the DiGraph
	 */
	public void removeNode(DiGraphNode<V, E> nodeToBeRemoved) {
		HashSet<DiGraphArc<V, E>> arcsToBeRemoved = new HashSet<DiGraphArc<V, E>>();
		arcsToBeRemoved.addAll(nodeToBeRemoved.incomingArcs);
		arcsToBeRemoved.addAll(nodeToBeRemoved.outgoingArcs);
		this.removeArcs(arcsToBeRemoved);
		nodeList.remove(nodeToBeRemoved);
	}

	/**
	 * method that returns arc between two DiGraphNodes in O(outDegree(source))
	 * 
	 * @param source - source node of the searched arc
	 * @param target - target node of the searched arc
	 * @return DiGraphArc between source and target if there is one
	 */
	public DiGraphArc<V, E> getArc(DiGraphNode<V, E> source, DiGraphNode<V, E> target) {
		for (DiGraphArc<V, E> arc : source.outgoingArcs)
			if (arc.getTarget() == target)
				return arc;
		return null;
	}

	/**
	 * method that returns a List of the arcs in a path between a list of
	 * DiGraphNodes
	 * 
	 * @param pathNodes - DiGraphNodes in the path
	 * @return List of DiGraphArcs visited in the path
	 */
	public LinkedList<DiGraphArc<V, E>> getPathArcs(List<DiGraphNode<V, E>> pathNodes) {
		LinkedList<DiGraphArc<V, E>> pathArcs = new LinkedList<DiGraphArc<V, E>>();
		DiGraphNode<V, E> u = null;
		for (DiGraphNode<V, E> v : pathNodes) {
			if (u != null) {
				// add arc uv to list
				DiGraphArc<V, E> uv = getArc(u, v);
				if (uv != null)
					pathArcs.add(uv);
			}
			u = v;
		}
		return pathArcs;
	}

	public ArrayList<ArrayList<Integer>> updateIDs() {
		ArrayList<ArrayList<Integer>> lookuptable = new ArrayList<>();
		ArrayList<Integer> nodeTable = new ArrayList<>();
		ArrayList<Integer> arcTable = new ArrayList<>();

		for (int i = 0; i < this.n(); i++) {
			DiGraphNode<V, E> node = this.getNode(i);
			nodeTable.add(node.getId());
			node.setId(i);
		}
		for (int i = 0; i < this.m(); i++) {
			DiGraphArc<V, E> arc = this.getArc(i);
			arcTable.add(arc.getId());
			arc.setId(i);
		}
		lookuptable.add(nodeTable);
		lookuptable.add(arcTable);

		return lookuptable;
	}

	@Override
	public void visitArcs(ArcVisitor visitor) {
		for (DiGraphArc<V, E> arc : arcList)
			visitor.visit(arc.source.id, arc.target.id, OUTGOING_ARC);
	}

	@Override
	public void visitIncidentArcs(int node, ArcVisitor visitor, byte filter) {
		if ((filter & INCOMING_ARC) != 0) {
			for (DiGraphArc<V, E> incArc : nodeList.get(node).incomingArcs) {
				visitor.visit(incArc.source.id, incArc.target.id, INCOMING_ARC);
			}
		}
		if ((filter & OUTGOING_ARC) != 0) {
			for (DiGraphArc<V, E> outArc : nodeList.get(node).outgoingArcs) {
				visitor.visit(outArc.source.id, outArc.target.id, OUTGOING_ARC);
			}
		}
	}

	@Override
	public void visitIncidentArcs(int node, ArcVisitor visitor) {
		this.visitIncidentArcs(node, visitor, (byte) (INCOMING_ARC | OUTGOING_ARC));
	}

	@Override
	public void visitOutgoingArcs(int node, ArcVisitor visitor) {
		this.visitIncidentArcs(node, visitor, OUTGOING_ARC);
	}

	@Override
	public void visitIncomingArcs(int node, ArcVisitor visitor) {
		this.visitIncidentArcs(node, visitor, INCOMING_ARC);
	}

	@Override
	public void visitArcs(ArcVisitorData<V, E> visitor) {
		for (DiGraphArc<V, E> arc : arcList) {
			visitor.visit(arc.source.id, arc.source.nodeData, arc.target.id, arc.target.nodeData, OUTGOING_ARC,
					arc.arcData);
		}
	}

	@Override
	public void visitIncidentArcs(int node, ArcVisitorData<V, E> visitor, byte filter) {
		if ((filter & INCOMING_ARC) != 0) {
			for (DiGraphArc<V, E> incArc : nodeList.get(node).incomingArcs) {
				visitor.visit(incArc.source.id, incArc.source.nodeData, incArc.target.id, incArc.target.getNodeData(),
						INCOMING_ARC, incArc.arcData);
			}
		}
		if ((filter & OUTGOING_ARC) != 0) {
			for (DiGraphArc<V, E> outArc : nodeList.get(node).outgoingArcs) {
				visitor.visit(outArc.source.id, outArc.source.nodeData, outArc.target.id, outArc.target.getNodeData(),
						OUTGOING_ARC, outArc.arcData);
			}
		}
	}

	@Override
	public void visitIncidentArcs(int node, ArcVisitorData<V, E> visitor) {
		this.visitIncidentArcs(node, visitor, (byte) (INCOMING_ARC | OUTGOING_ARC));
	}

	@Override
	public void visitOutgoingArcs(int node, ArcVisitorData<V, E> visitor) {
		this.visitIncidentArcs(node, visitor, OUTGOING_ARC);
	}

	@Override
	public void visitIncomingArcs(int node, ArcVisitorData<V, E> visitor) {
		this.visitIncidentArcs(node, visitor, INCOMING_ARC);
	}

	@Override
	public void visitNodes(NodeVisitorData<V> visitor) {
		for (DiGraphNode<V, E> node : nodeList) {
			visitor.visit(node.getId(), node.getNodeData());
		}
	}

	/**
	 * Getter-method of the NodeData
	 * 
	 * @param node - id of the DiGraphNode
	 * @return NodeData of the DiGraphNode with id node
	 */
	@Override
	public V getNodeData(int node) {
		return getNode(node).getNodeData();
	}
}
