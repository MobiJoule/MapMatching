package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.core.geometry.PointComparator;

public class CandidateMatch<I> implements Comparable<CandidateMatch<I>> {

	private Point2D gpsPoint = null;
	private Point2D mapPoint = null;
	private DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> myNode = null;

	public DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> mapSegment = null;

	/**
	 * Eine (moegliche) Zuordnung zwischen einem Strassenpunkt und einem Trackpunkt
	 * 
	 * @param gpspoint der Trackpunkt
	 * @param point    der Strassenpunkt
	 * @param s        das Strassensegment, auf dem point liegt
	 */
	public CandidateMatch(Point2D gpspoint, Point2D point, DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> s) {
		this.gpsPoint = gpspoint;
		this.mapPoint = point;
		this.mapSegment = s;
	}

	public DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> getNode() {
		return myNode;
	}

	public void setNode(DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> myNode) {
		this.myNode = myNode;
	}

	@Override
	public int compareTo(CandidateMatch<I> o) {
		Point2D mapPoint1 = this.mapPoint;
		Point2D mapPoint2 = o.mapPoint;

		PointComparator myPointLexOrder = new PointComparator();
		int c = myPointLexOrder.compare(mapPoint1, mapPoint2);
		if (c != 0)
			return c;

		Point2D gpsPoint1 = this.gpsPoint;
		Point2D gpsPoint2 = o.gpsPoint;
		return myPointLexOrder.compare(gpsPoint1, gpsPoint2);

	}

	/**
	 * By default, candidate cost is the square of the euclidean distance between
	 * gps point and track point
	 * 
	 * @return square of euclidean distance between gps point and track point
	 */
	public double getCandidateCost() {
		double dx = this.gpsPoint.getX() - this.mapPoint.getX();
		double dy = this.gpsPoint.getY() - this.mapPoint.getY();
		return dx * dx + dy * dy;
	}

	public Point2D getMapPoint() {
		return mapPoint;
	}

	public Point2D getGpsPoint() {
		return gpsPoint;
	}

	public I getSegmentType() {
		if (mapSegment == null)
			return null;
		return mapSegment.getArcData().getInfo();
	}

	public DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> getSegment() {
		return mapSegment;
	}

	@Override
	public String toString() {
		return "CandidateMatch [gpsPoint=" + gpsPoint + ",mapPoint=" + mapPoint + ",cost=" + getCandidateCost() + "]";
	}
}
