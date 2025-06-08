package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;
import java.util.List;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;

public class WeightedPathToCandidate<I> {
	private List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p;
	private List<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> arcs;
	private double d;
	private CandidateMatch<I> cm;

	public WeightedPathToCandidate(List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p, double d,
			CandidateMatch<I> cm) {
		this.p = p;
		this.d = d;
		this.cm = cm;
	}

	public WeightedPathToCandidate(List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p, double d,
			CandidateMatch<I> cm, List<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> arcs) {
		this.p = p;
		this.d = d;
		this.cm = cm;
		this.arcs = arcs;
	}

	public List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> getPath() {
		return p;
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

	public List<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> getArcs() {
		return arcs;
	}

	@Override
	public String toString() {
		return "WPTC [d=" + d + ",mapPoint=" + cm.getMapPoint() + "]";
	}
}
