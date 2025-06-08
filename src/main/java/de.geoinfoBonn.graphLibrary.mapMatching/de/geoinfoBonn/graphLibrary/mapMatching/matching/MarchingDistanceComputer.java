package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;

import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.core.shortestPath.DistanceComputer;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.Typed;

public class MarchingDistanceComputer<I extends Typed>
		implements DistanceComputer<Point2D, DoubleWeightDataWithInfo<I>> {

	private double candidateCost;
	private double unmatchedCost;
	private double tessalationCost;
	private double obstacleBoundaryCost;

	public MarchingDistanceComputer(double candidateCost, double unmatchedCost, double tessalationCost,
			double obstacleBoundaryCost) {
		this.candidateCost = candidateCost;
		this.unmatchedCost = unmatchedCost;
		this.tessalationCost = tessalationCost;
		this.obstacleBoundaryCost = obstacleBoundaryCost;
	}

	@Override
	public double computeDistance(int uId, Point2D u, int vId, Point2D v, DoubleWeightDataWithInfo<I> e) {
		double typeMultiplier;
		switch (e.getInfo().getType()) {
		case TESSA:
			typeMultiplier = tessalationCost;
			break;
		case UNMATCHED:
			typeMultiplier = unmatchedCost;
			break;
		case HOLE:
			typeMultiplier = obstacleBoundaryCost;
			break;
		case BOUNDARY:
			typeMultiplier = tessalationCost;
			break;
		case CANDIDATE:
			typeMultiplier = candidateCost;
			break;
		case ROAD:
			typeMultiplier = 1;
			break;
		default:
			System.err.println("DEFAULT! " + e.getInfo());
			typeMultiplier = 1;
		}
		return typeMultiplier * e.getValue();
	}
}
