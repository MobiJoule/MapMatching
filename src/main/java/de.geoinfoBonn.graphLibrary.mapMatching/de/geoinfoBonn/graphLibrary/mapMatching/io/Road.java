package de.geoinfoBonn.graphLibrary.mapMatching.io;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.TreeMap;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.core.geometry.PointComparator;

/**
 * A polyline with an additional info attribute, e.g., representing a road type
 * 
 * @author haunert
 */

public class Road<I> {

	/**
	 * the road geometry as a sequence of points
	 */
	private final ArrayList<Point2D> points;

	/**
	 * the road type, e.g., motorway
	 */
	private final I info;

	private final double weight;

	/**
	 * constructor for setting road attributes
	 * 
	 * @param points the sequence of points defining the road's geometry
	 * @param info   information about the road (e.g. type)
	 * @param weight weight/cost for traversing the road
	 */
	public Road(ArrayList<Point2D> points, I info, double weight) {
		this.points = points;
		this.info = info;
		this.weight = weight;
	}

	/**
	 * getter for accessing point sequence
	 * 
	 * @return the sequence of points defining the road's geometry
	 */
	public ArrayList<Point2D> getPoints() {
		return points;
	}

	/**
	 * getter for accessing road type
	 * 
	 * @return the road type
	 */
	public I getInfo() {
		return info;
	}

	/**
	 * method for computing a graph from a set of roads
	 * 
	 * @param roads: a list of roads
	 * @return a DiGraph representing the road network
	 */
	public static <I> DiGraph<Point2D, DoubleWeightDataWithInfo<I>> buildGraph(List<Road<I>> roads) {

		DiGraph<Point2D, DoubleWeightDataWithInfo<I>> g = new DiGraph<>();

		Comparator<Point2D> comp = new PointComparator();

		TreeMap<Point2D, DiGraph.DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> nodesByPoints = new TreeMap<>(comp);

		for (Road<I> r : roads) {

			// compute length of road
			Point2D prevPoint = null;
			double roadLength = 0.0;
			for (Point2D currentPoint : r.getPoints()) {
				if (prevPoint != null) {
					roadLength += currentPoint.distance(prevPoint);
				}
				prevPoint = currentPoint;
			}

			DiGraph.DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> prevNode = null;
			for (Point2D currentPoint : r.getPoints()) {
				DiGraph.DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> currentNode = nodesByPoints.get(currentPoint);
				if (currentNode == null) {
					currentNode = g.addNode(currentPoint);
					nodesByPoints.put(currentPoint, currentNode);
				}
				if (prevNode != null) {
					double segLength = currentPoint.distance(prevNode.getNodeData());
					DoubleWeightDataWithInfo<I> arcData = new DoubleWeightDataWithInfo<>(
							r.weight * segLength / roadLength, r.getInfo());
					g.addArc(prevNode, currentNode, arcData);
//					g.addArc(currentNode, prevNode, arcData);  // comment out to make all links bidirectional
				}
				prevNode = currentNode;
			}
		}
		return g;
	}
}