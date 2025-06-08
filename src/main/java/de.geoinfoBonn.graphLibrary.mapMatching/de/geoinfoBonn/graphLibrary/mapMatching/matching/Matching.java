package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.strtree.STRtree;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.core.geometry.Calculations;
import de.geoinfoBonn.graphLibrary.core.geometry.PointComparator;
import de.geoinfoBonn.graphLibrary.core.shortestPath.Dijkstra;
import de.geoinfoBonn.graphLibrary.core.shortestPath.MultiTargetNodeVisitor;

public class Matching<I> {

	private final LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> path;
	private final LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> pathArcs;
	private final LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> matches;
	private ArrayList<Point2D> track;
	private final HashSet<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> offRoadNodes;
	private final DiGraph<Point2D, DoubleWeightDataWithInfo<I>> g;
	private final int oldNumberOfArcs;
	private final int oldNumberOfNodes;

	public static double RADIUS = 25.0;
	public static int MAX_CAND_N = Integer.MAX_VALUE;
	public static double CANDIDATE_COST_WEIGHT = 0.01;
	public static double OFF_ROAD_WEIGHT = 1.5;
	public static boolean ADD_OFFROAD_CANDIDATE = true;
	public static boolean VERBOSE = false;

	private static final PointComparator POINT_LEX_ORDER = new PointComparator();

	/**
	 * @param input_graph a geometric graph representing the road network
	 * @param input_track a series of points representing the GPS trajectory
	 */
	public Matching(DiGraph<Point2D, DoubleWeightDataWithInfo<I>> input_graph, ArrayList<Point2D> input_track) {

		if (Matching.VERBOSE) {
			System.out.println("#################################################################");
			System.out.println("RADIUS: \t\t" + RADIUS);
			System.out.println("MAX_CAND_N: \t\t" + MAX_CAND_N);
			System.out.println("CANDIDATE_COST_WEIGHT: \t" + CANDIDATE_COST_WEIGHT);
			System.out.println("OFF_ROAD_WEIGHT: \t" + OFF_ROAD_WEIGHT);
			System.out.println("ADD_OFFROAD_CANDIDATE: \t" + ADD_OFFROAD_CANDIDATE);
			System.out.println("#################################################################");
		}

		g = input_graph;
		oldNumberOfArcs = g.getArcs().size();
		oldNumberOfNodes = g.getNodes().size();
		if (VERBOSE)
			System.out.println("Size of graph before node insertion:" + oldNumberOfNodes + " " + oldNumberOfArcs);

		// compute strtree with segments of graph (for double-edges: only arcs pointing
		// from left to right)
		STRtree segments = new STRtree();
		for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> s : g.getArcs()) {
			if (s.getTwin() == null
					|| POINT_LEX_ORDER.compare(s.getSource().getNodeData(), s.getTarget().getNodeData()) < 0) {
				Envelope env = new Envelope();
				env.expandToInclude(s.getSource().getNodeData().getX(), s.getSource().getNodeData().getY());
				env.expandToInclude(s.getTarget().getNodeData().getX(), s.getTarget().getNodeData().getY());
				segments.insert(env, s);
			}
		}

		// find for each gps point the k best candidates
		ArrayList<LinkedList<CandidateMatch<I>>> candidates = getBestKCandidatesForEachTrackPoint(segments, input_track,
				RADIUS, MAX_CAND_N);

		// add off-road candidates to candidate lists
		track = null;
		offRoadNodes = new HashSet<>();
		if (ADD_OFFROAD_CANDIDATE) {
			track = input_track;
			for (int i = 0; i < input_track.size(); i++) {
				CandidateMatch<I> offRoadCandidate = new CandidateMatch<>(input_track.get(i), input_track.get(i),
						null);
				candidates.get(i).add(offRoadCandidate);
				DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> offRoadNode = g.addNode(input_track.get(i));
				offRoadCandidate.setNode(offRoadNode);
				offRoadNodes.add(offRoadNode);
			}
		} else {
			// reduce track to gps points that have at least one candidate
			track = new ArrayList<>();
			ArrayList<LinkedList<CandidateMatch<I>>> candidates_new = new ArrayList<>();
			for (int i = 0; i < input_track.size(); i++) {
				if (!candidates.get(i).isEmpty()) {
					track.add(input_track.get(i));
					candidates_new.add(candidates.get(i));
				}
			}
			candidates = candidates_new;
		}

		// compute for each arc of the graph the candidates on it
		HashMap<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>, LinkedList<CandidateMatch<I>>> candidatesToInject = new HashMap<>();
		for (LinkedList<CandidateMatch<I>> candidatesOfSamePoint : candidates) {
			for (CandidateMatch<I> c : candidatesOfSamePoint) {
				if (c.mapSegment != null) {
                    LinkedList<CandidateMatch<I>> pointsOnSeg = candidatesToInject.computeIfAbsent(c.mapSegment, k -> new LinkedList<>());
                    pointsOnSeg.add(c);
				}
			}
		}

		// add candidate points as new graph nodes and arcs connecting them
		for (Entry<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>, LinkedList<CandidateMatch<I>>> e : candidatesToInject
				.entrySet()) {
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> arc1 = e.getKey();
			LinkedList<CandidateMatch<I>> candsOnArc1 = e.getValue();

			if (POINT_LEX_ORDER.compare(arc1.getSource().getNodeData(), arc1.getTarget().getNodeData()) < 0) {
				Collections.sort(candsOnArc1);
			} else {
				Collections.sort(candsOnArc1, Collections.reverseOrder());
			}

			ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> newNodes = new ArrayList<>();
			newNodes.add(arc1.getSource());
			for (CandidateMatch<I> p : candsOnArc1) {
				DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> newNode = g.addNode(p.getMapPoint());
				p.setNode(newNode);
				newNodes.add(newNode);
			}
			newNodes.add(arc1.getTarget());

			double oldArcLength = arc1.getSource().getNodeData().distance(arc1.getTarget().getNodeData());
			for (int i = 1; i < newNodes.size(); i++) {
				double newArcLength = newNodes.get(i).getNodeData().distance(newNodes.get(i - 1).getNodeData());
				DoubleWeightDataWithInfo<I> dw1 = new DoubleWeightDataWithInfo<>(
						arc1.getArcData().getValue() * newArcLength / oldArcLength, arc1.getArcData().getInfo());
				g.addArc(newNodes.get(i - 1), newNodes.get(i), dw1);

				if (arc1.getTwin() != null) {
					DoubleWeightDataWithInfo<I> dw2 = new DoubleWeightDataWithInfo<I>(
							arc1.getTwin().getArcData().getValue() * newArcLength / oldArcLength,
							arc1.getTwin().getArcData().getInfo());
					g.addArc(newNodes.get(i), newNodes.get(i - 1), dw2);
				}
			}
		}

		if (VERBOSE)
			System.out.println("Size of graph after node insertion: " + g.getNodes().size() + " " + g.getArcs().size());

		// add off-road edges
		if (ADD_OFFROAD_CANDIDATE) {
			for (int i = 0; i < input_track.size(); i++) {
				CandidateMatch<I> c1 = candidates.get(i).getLast(); // fetch off-road candidate for i-th track point
				Point2D p1 = c1.getGpsPoint();
				if (i < input_track.size() - 1) {
					for (CandidateMatch<I> c2 : candidates.get(i + 1)) {
						Point2D p2 = c2.getMapPoint();
						g.addArc(c1.getNode(), c2.getNode(),
								new DoubleWeightDataWithInfo<>(OFF_ROAD_WEIGHT * p1.distance(p2), null));
					}
				}
				if (i > 0) {
					for (CandidateMatch<I> c0 : candidates.get(i - 1)) {
						if (c0.mapSegment != null) {
							Point2D p0 = c0.getMapPoint();
							g.addArc(c0.getNode(), c1.getNode(),
									new DoubleWeightDataWithInfo<>(OFF_ROAD_WEIGHT * p0.distance(p1), null));
						}
					}
				}
			}
			if (VERBOSE)
				System.out.println(
						"Size of graph with off-road edges: " + g.getNodes().size() + " " + g.getArcs().size());

		}

		// add dummy nodes and arcs
		ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> dummies = new ArrayList<>();
		for (int i = 0; i < candidates.size() - 1; i++) {
			LinkedList<CandidateMatch<I>> candidatesForCurrentPoint = candidates.get(i);
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> dummy = g.addNode(null);
			dummies.add(dummy);
			for (CandidateMatch<I> cm : candidatesForCurrentPoint) {
				if(i == 0) {
					g.addArc(dummy, cm.getNode(),
							new DoubleWeightDataWithInfo<>(CANDIDATE_COST_WEIGHT * 10 * cm.getCandidateCost(), null));
				} else {
					g.addArc(dummy, cm.getNode(),
							new DoubleWeightDataWithInfo<>(CANDIDATE_COST_WEIGHT * cm.getCandidateCost(), null));

				}
			}
		}

		// initialize dijkstra instance
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);

		// search paths
		ArrayList<LinkedList<WeightedPathToCandidate<I>>> allWPs = new ArrayList<LinkedList<WeightedPathToCandidate<I>>>();
		for (int i = 0; i < candidates.size() - 1; i++) {

			// define source
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> dummy = dummies.get(i);

			// update weights of dummy arcs
			if (i > 0) {
				for (WeightedPathToCandidate<I> myPath : allWPs.get(i - 1)) {
					DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> dummyArc = g.getArc(dummy, myPath.getTarget());
					if (dummyArc != null) {
						dummyArc.setArcData(new DoubleWeightDataWithInfo<>(
								dummyArc.getArcData().getValue() + myPath.getDistance(), null));
					}
				}
			}

			// define targets
			LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> targets = new LinkedList<>();
			for (CandidateMatch<I> cm : candidates.get(i + 1)) {
				targets.add(cm.getNode());
			}

			// let dijkstra run from source to targets
			dijkstra.run(dummy, new MultiTargetNodeVisitor<>(targets, dijkstra));

			// memorize solutions of dijkstra
			LinkedList<WeightedPathToCandidate<I>> myPathList = new LinkedList<>();
			for (CandidateMatch<I> cm : candidates.get(i + 1)) {
				List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p = g
						.toNodeList(dijkstra.getPath(cm.getNode().getId()));
				double d = dijkstra.getDistance(cm.getNode());
				WeightedPathToCandidate<I> wp = new WeightedPathToCandidate<>(p, d, cm);
				myPathList.add(wp);
			}
			allWPs.add(myPathList);
		}

		// find last candidate of optimal solution
		double minTotalCost = Double.POSITIVE_INFINITY;
		CandidateMatch<I> bestCandidate = null;
		LinkedList<WeightedPathToCandidate<I>> myPathList = allWPs.get(candidates.size() - 2);
		for (WeightedPathToCandidate<I> myPath : myPathList) {
			double totalCost = myPath.getDistance() + myPath.getTargetCandidate().getCandidateCost();
			if (totalCost < minTotalCost) {
				minTotalCost = totalCost;
				bestCandidate = myPath.getTargetCandidate();
			}
		}
		if (VERBOSE)
			System.out.println("total cost: " + minTotalCost);

		// reconstruct solution
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> lastNode = bestCandidate.getNode();
		path = new LinkedList<>();
		pathArcs = new LinkedList<>();
		matches = new LinkedList<>();
		matches.add(lastNode);

		for (int i = candidates.size() - 2; i >= 0; i--) {
			for (WeightedPathToCandidate<I> myWP : allWPs.get(i)) {
				if (myWP.getTarget() == lastNode) {
					List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> subPath = myWP.getPath();
					subPath.remove(0); // remove dummy node
					pathArcs.addAll(0, g.getPathArcs(subPath));
					path.addAll(0, subPath);
					lastNode = subPath.get(0);
					matches.add(0, lastNode);
					break;
				}
			}
		}

		// print some statistics
//		if (VERBOSE) {
//			System.out.println("Path: ");
//			for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> v : path) {
//				System.out.println(v.getId() + " " + v.getNodeData());
//			}
//		}

	}

	/**
	 * Method to find for each track point the k nearest segments of the graph and
	 * the closest point on each of those segments
	 * 
	 * @param segments  a quadtree with the segments of the network
	 * @param gps_track a list with the track points
	 * @param r         a search radius
	 * @param k         the number of candidates per track point
	 * @return for each track point a list of candidate matches
	 */
	private ArrayList<LinkedList<CandidateMatch<I>>> getBestKCandidatesForEachTrackPoint(STRtree segments,
			ArrayList<Point2D> gps_track, double r, int k) {
		ArrayList<LinkedList<CandidateMatch<I>>> allCandidates = new ArrayList<>();
		for (Point2D gps_point : gps_track) {
			LinkedList<CandidateMatch<I>> candidates = getBestKCandidatesForTrackPoint(segments, gps_point, r, k);
			allCandidates.add(candidates);
		}
		return allCandidates;
	}

	/**
	 * 
	 * @param segments
	 * @param gps_point
	 * @param r
	 * @param k
	 * @return
	 */
	private LinkedList<CandidateMatch<I>> getBestKCandidatesForTrackPoint(STRtree segments, Point2D gps_point, double r,
			int k) {
		LinkedList<CandidateMatch<I>> candidates = new LinkedList<>();

		// square-shaped search window around gps_point
		Envelope searchEnvelope = new Envelope(gps_point.getX() - r, gps_point.getX() + r, gps_point.getY() - r,
				gps_point.getY() + r);

		// test each segment in search window
		for (Object o : segments.query(searchEnvelope)) {
			@SuppressWarnings("unchecked")
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> seg = (DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>) o;

			// find point on seg closest to gps_point
			Line2D line = new Line2D.Double(seg.getSource().getNodeData(), seg.getTarget().getNodeData());
			Point2D nearestPoint = Calculations.nearestPointOnSegment(gps_point, line);

			// if point on seg is within search radius, add it to candidate list
			if (Math.hypot(gps_point.getX() - nearestPoint.getX(), gps_point.getY() - nearestPoint.getY()) <= r) {
				CandidateMatch<I> cm = new CandidateMatch<>(gps_point, nearestPoint, seg);
				candidates.add(cm);
			}
		}

		// sort candidates by distance
		Comparator<CandidateMatch<I>> myCandidateComp = (a, b) -> {
            if (a.getMapPoint().distance(a.getGpsPoint()) < b.getMapPoint().distance(b.getGpsPoint()))
                return -1;
            if (a.getMapPoint().distance(a.getGpsPoint()) > b.getMapPoint().distance(b.getGpsPoint()))
                return 1;
            return POINT_LEX_ORDER.compare(a.getMapPoint(), b.getMapPoint());
        };
		candidates.sort(myCandidateComp);

		// remove all but the best k candidates
		while (candidates.size() > MAX_CAND_N)
			candidates.removeLast();
		return candidates;
	}

	/**
	 * Method to restore the graph by deleting all new nodes and their incident
	 * edges that were added for map matching
	 */
	public void restoreGraph() {
		for (int n = g.getArcs().size(); n > oldNumberOfArcs; n = g.getArcs().size()) {
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> a = g.getArcs().get(n - 1);
			a.getSource().getOutgoingArcs().remove(a);
			a.getTarget().getIncomingArcs().remove(a);
			g.getArcs().remove(n - 1);
		}

		while (g.getNodes().size() > oldNumberOfNodes) {
			g.getNodes().remove(g.getNodes().size() - 1);
		}
	}

	/**
	 * Method to retrieve the matching solution
	 * 
	 * @return the path in the graph corresponding to the trajectory
	 */
	public LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> getPath() {
		return path;
	}

	public ArrayList<ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> getChunks() {
		ArrayList<ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> chunks = new ArrayList<>();
		ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> currentChunk = new ArrayList<>();
		for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> node : getPath()) {
			if (!offRoadNodes.contains(node)) {
				currentChunk.add(node);
			} else {
				if (!currentChunk.isEmpty()) {
					if (currentChunk.size() > 1) {
						chunks.add(currentChunk);
					}
					currentChunk = new ArrayList<>();
				}
			}
		}
		if (currentChunk.size() > 1) {
			chunks.add(currentChunk);
		}
		return chunks;
	}

	public List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> getShortestPathForWholeTrajectory() {
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> s = path.get(0);
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> t = path.get(path.size() - 1);
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);
		dijkstra.run(s, t);
		return g.toNodeList(dijkstra.getPath(t.getId()));
	}

	public ArrayList<List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> getShortestPathsForChunks() {
		ArrayList<List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> shortestPaths = new ArrayList<>();
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);
		for (ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> chunk : this.getChunks()) {
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> s = chunk.get(0);
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> t = chunk.get(chunk.size() - 1);
			dijkstra.run(s.getId(), t.getId());
			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> shortestpath = g
					.toNodeList(dijkstra.getPath(t.getId()));
			if (shortestpath.size() > 1) {
				shortestPaths.add(shortestpath);
			}
		}
		return shortestPaths;
	}

	public LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> getMatches() {
		return matches;
	}

	public ArrayList<Point2D> getTrack() {
		return track;
	}

	public LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> getPathArcs() {
		return pathArcs;
	}

}
