package de.geoinfoBonn.graphLibrary.mapMatching.core.generic;

import java.util.List;

public interface VisitableGraph<V, E> {

	// Declaration of flag constants
	public static final byte INCOMING_ARC = 1;
	public static final byte OUTGOING_ARC = 2;

	/*
	 * -----------------------------------------------------------------------------
	 * ----- ----------------------- VISITOR DECLARATIONS
	 * -------------------------------------
	 * -----------------------------------------------------------------------------
	 * -----
	 */
	public static interface NodeVisitorData<V> {
		void visit(int nodeId, V nodeData);
	}

	public static interface ArcVisitor {
		void visit(int source, int target, byte flags);
	}

	public static interface ArcVisitorData<V, E> {
		void visit(int source, V sourceData, int target, V targetData, byte flags, E data);
	}

	public static interface ArcVisitorIndex {
		public void visit(int source, int target, int flags, int arcDataIndex);
	}

	public interface GraphFactory<V, E> {
		int[] getArcs(int nodeId);

		byte[] getFlags(int nodeId);

		List<E> getData(int nodeId);
	}

	/*
	 * ------- NODE VISITOR DATA -------
	 */
	public void visitNodes(NodeVisitorData<V> visitor);

	public V getNodeData(int node);

	/*
	 * ------- ARC VISITOR -------
	 */
	public void visitArcs(ArcVisitor visitor);

	public void visitIncidentArcs(int node, ArcVisitor visitor, byte filter);

	public void visitIncidentArcs(int node, ArcVisitor visitor);

	public void visitOutgoingArcs(int node, ArcVisitor visitor);

	public void visitIncomingArcs(int node, ArcVisitor visitor);

	/*
	 * ------- ARC VISITOR DATA -------
	 */
	public void visitArcs(ArcVisitorData<V, E> visitor);

	public void visitIncidentArcs(int node, ArcVisitorData<V, E> visitor, byte filter);

	public void visitIncidentArcs(int node, ArcVisitorData<V, E> visitor);

	public void visitOutgoingArcs(int node, ArcVisitorData<V, E> visitor);

	public void visitIncomingArcs(int node, ArcVisitorData<V, E> visitor);

	/*
	 * ------- GETTER & SETTER -------
	 */
	public int n(); // Number of nodes

	public int m(); // Number of edges
}
