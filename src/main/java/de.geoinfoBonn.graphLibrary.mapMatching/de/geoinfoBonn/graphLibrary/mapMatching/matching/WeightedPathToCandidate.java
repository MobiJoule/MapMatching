package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DoubleWeightDataWithInfo;

public class WeightedPathToCandidate<I> {
	private final DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> source;
	private final double d;
	private final CandidateMatch<I> cm;

	public WeightedPathToCandidate(DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> source, double d,
			CandidateMatch<I> cm) {
		this.source = source;
		this.d = d;
		this.cm = cm;
	}

	public DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> getSource() {
		return source;
	}

	public double getDistance() {
		return d;
	}

	public DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> getTarget() {
		return cm.getNode();
	}

	public CandidateMatch<I> getTargetCandidate() {
		return cm;
	}

	@Override
	public String toString() {
		return "WPTC [d=" + d + ",mapPoint=" + cm.getMapPoint() + "]";
	}
}
