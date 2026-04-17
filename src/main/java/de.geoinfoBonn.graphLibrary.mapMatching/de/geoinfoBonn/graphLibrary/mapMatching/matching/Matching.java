package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;
import java.util.List;
import java.util.Map.Entry;

import de.geoinfoBonn.graphLibrary.mapMatching.io.Road.RoadInfo;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.strtree.STRtree;

import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.mapMatching.core.geometry.Calculations;
import de.geoinfoBonn.graphLibrary.mapMatching.core.geometry.PointComparator;
import de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath.Dijkstra;
import de.geoinfoBonn.graphLibrary.mapMatching.core.shortestPath.MultiTargetNodeVisitor;
import org.tinylog.Logger;

public class Matching<I> {

	private final LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> path;
	private final LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> pathArcs;
	private final LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> matches;
	private final LinkedList<Integer> matchesArcCount;
	private ArrayList<Point2D> track;
	private final HashSet<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> offRoadNodes;
	private final DiGraph<Point2D, DoubleWeightDataWithInfo<I>> g;
	private final int oldNumberOfArcs;
	private final int oldNumberOfNodes;

	public static double RADIUS = 100.0;
	public static int MAX_CAND_N = Integer.MAX_VALUE; // has been 50 so far...
	public static double CANDIDATE_COST_WEIGHT = 0.01;
	public static double OFF_ROAD_WEIGHT = 15;
	public static boolean ADD_OFFROAD_CANDIDATE = true;
	public static boolean VERBOSE = false;

	public static double DEVIATION_PENALTY_FACTOR = 1.4; // must be >= 1
	public static double DISTANCE_PENALTY_FACTOR = 0.6; // must be <= 1

	private static final PointComparator POINT_LEX_ORDER = new PointComparator();

	/**
	 * @param input_graph a geometric graph representing the road network
	 * @param input_track input track
	 */
	public Matching(DiGraph<Point2D, DoubleWeightDataWithInfo<I>> input_graph,
					Track input_track,
					I offroadLinkInfo) {

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

		// compute strtree with segments of graph (for double-edges: only arcs pointing from left to right)
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
		ArrayList<Point2D> input_track_points = input_track.getTrackPoints();
		ArrayList<LinkedList<CandidateMatch<I>>> candidates = getBestKCandidatesForEachTrackPoint(segments, input_track_points);

		// add off-road candidates to candidate lists
		track = null;
		offRoadNodes = new HashSet<>();
		if (ADD_OFFROAD_CANDIDATE) {
			track = input_track_points;
			for (int i = 0; i < input_track_points.size(); i++) {
				CandidateMatch<I> offRoadCandidate = new CandidateMatch<>(input_track_points.get(i), input_track_points.get(i),
						null);
				candidates.get(i).add(offRoadCandidate);
				DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> offRoadNode = g.addNode(input_track_points.get(i));
				offRoadCandidate.setNode(offRoadNode);
				offRoadNodes.add(offRoadNode);
			}
		} else {
			// reduce track to gps points that have at least one candidate
			track = new ArrayList<>();
			ArrayList<LinkedList<CandidateMatch<I>>> candidates_new = new ArrayList<>();
			for (int i = 0; i < input_track_points.size(); i++) {
				if (!candidates.get(i).isEmpty()) {
					track.add(input_track_points.get(i));
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
				candsOnArc1.sort(Collections.reverseOrder());
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
					DoubleWeightDataWithInfo<I> dw2 = new DoubleWeightDataWithInfo<>(
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
			for (int i = 0; i < input_track_points.size(); i++) {
				CandidateMatch<I> c1 = candidates.get(i).getLast(); // fetch off-road candidate for i-th track point
				Point2D p1 = c1.getGpsPoint();
				if (i < input_track_points.size() - 1) {
					for (CandidateMatch<I> c2 : candidates.get(i + 1)) {
						Point2D p2 = c2.getMapPoint();
						g.addArc(c1.getNode(), c2.getNode(),
								new DoubleWeightDataWithInfo<>(OFF_ROAD_WEIGHT * p1.distance(p2), offroadLinkInfo));
					}
				}
				if (i > 0) {
					for (CandidateMatch<I> c0 : candidates.get(i - 1)) {
						if (c0.mapSegment != null) {
							Point2D p0 = c0.getMapPoint();
							g.addArc(c0.getNode(), c1.getNode(),
									new DoubleWeightDataWithInfo<>(OFF_ROAD_WEIGHT * p0.distance(p1), offroadLinkInfo));
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
				double dummyWeight = CANDIDATE_COST_WEIGHT * cm.getCandidateCost();
				if(i == 0) {
					dummyWeight *= 10;
				}
				g.addArc(dummy, cm.getNode(), new DoubleWeightDataWithInfo<>(dummyWeight, null));
			}
		}

		// initialize dijkstra instance
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);

		// search paths
		ArrayList<LinkedList<WeightedPathToCandidate<I>>> allWPs = new ArrayList<>();
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

			// Current and next trajectory point
			Point2D currPoint = input_track_points.get(i);
			Point2D nextPoint = input_track_points.get(i + 1);

			// Compute expected distance if speeds are available
			Double dt = null;
			Double modelled_speed = input_track.getSpeed(i+1);
			if (modelled_speed != null) {
				dt = input_track.getDiffTime(i+1).doubleValue();
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
				List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p = g.toNodeList(dijkstra.getPath(cm.getNode().getId()));
				double cost = dijkstra.getCost(cm.getNode());

				// Get current and next match
				Point2D currMatch = p.get(1).getNodeData();
				Point2D nextMatch = p.getLast().getNodeData();

				//Get path arcs
				LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> pathArcs = getPathArcs(p);

				// Test cost
				double costDummy = dijkstra.getCost(p.get(1));
				double costNet = cost - costDummy;
				double costNetTest = pathArcs.stream().mapToDouble(
						a -> a.getArcData().getValue()).sum();
				if (Math.abs(costNet - costNetTest) > 0.001d) {
					throw new RuntimeException("Cost test failed! costMap="+costNet+". costMapTest="+costNetTest+
							". This is probably happening because costs are running out of control...");
				}

				// Initialise penalties
				double deviation_penalty = 0.;
				double distance_penalty = 0.;

				// DEVIATION PENALTY (FOR SMOOTHED TRAJECTORIES ONLY)
				if (dt != null) {

					// X and Y components of next match vector
					double bx = nextMatch.getX() - nextPoint.getX();
					double by = nextMatch.getY() - nextPoint.getY();

					// Network distance distance
					double networkDistance = pathArcs.stream().mapToDouble(
							a -> a.getArcData().getValue() /
									((RoadInfo) a.getArcData().getInfo()).getWeightAdjustment()).sum();
					double networkDistanceSq = networkDistance * networkDistance;

					// Current match distance
					double diff_parallel;
					double diff_perpendicular;
					double diff_network;
					if (currMatch.equals(currPoint)) {
						diff_parallel = bx * bx + by * by;
						diff_perpendicular = 0.;
						diff_network = networkDistanceSq;
					} else {

						// Euclidean distance (squared)
						double euclideanDistanceSq = currPoint.distanceSq(nextPoint);

						// X and Y components of current match vector
						double ax = currMatch.getX() - currPoint.getX();
						double ay = currMatch.getY() - currPoint.getY();

						// Current match vector distance (squared)
						double a2 = ax * ax + ay * ay;

						// Parallel and orthogonal distance (squared) of next match vector, projected to current vector
						double dot = ax * bx + ay * by;
						double det = ay * bx - ax * by;

						double bp2 = dot * dot / a2;
						double br2 = det * det / a2;

						if (Math.abs(Math.atan2(det,dot)) > Math.PI / 2) {
							bp2 *= -1;
						}

						// Perpendicular and parallel components
						diff_parallel = Math.abs(a2 - bp2);
						diff_perpendicular = Math.abs(br2);

						// Network vs euclidean distance
						diff_network = Math.abs(networkDistanceSq - euclideanDistanceSq);

//						if (networkDistanceSq > euclideanDistanceSq) {
//							boolean hasOffroad = pathArcs.stream().anyMatch(a -> a.getArcData().getInfo().equals(offroadLinkInfo));
//
//							if (!hasOffroad) {
//								double excess = networkDistance - Math.sqrt(euclideanDistanceSq);
//
//								DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> firstArc = pathArcs.getFirst();
//								DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> lastArc = pathArcs.getLast();
//
//								// Match distance
//								double d_curr = Math.min(15., currPoint.distance(currMatch));
//								double d_next = Math.min(15., nextPoint.distance(nextMatch));
//
//								// Angle between vectors
//								double cx = firstArc.getTarget().getNodeData().getX() - firstArc.getSource().getNodeData().getX();
//								double cy = firstArc.getTarget().getNodeData().getY() - firstArc.getSource().getNodeData().getY();
//								double dx = lastArc.getTarget().getNodeData().getX() - lastArc.getSource().getNodeData().getX();
//								double dy = lastArc.getTarget().getNodeData().getY() - lastArc.getSource().getNodeData().getY();
//								double angle = Math.atan2(cx * dy - dx * cy, cx * dx + cy * dy);
//
//								// Turning allowance for right turns only
//								if (angle < 0) {
//									angle = Math.max(angle, Math.PI / -2.);
//									networkDistance += Math.max(-1 * excess, (d_curr + d_next) * Math.tan(angle / 2));
//									diff_network = networkDistance*networkDistance - euclideanDistanceSq;
//								}
//							}
//						}

					}

					deviation_penalty = (diff_parallel + diff_perpendicular) * CANDIDATE_COST_WEIGHT * DEVIATION_PENALTY_FACTOR / dt;
					distance_penalty = diff_network * CANDIDATE_COST_WEIGHT * DISTANCE_PENALTY_FACTOR;
				}

				// NETWORK DISTANCE PENALTY (ONLY WHEN EXPECTED DISTANCE AVAILABLE)
				// This doesn't really do anything. Square before taking difference instead...
//				if (expected_distance != null) {
//					double network_distance = pathArcs.stream().mapToDouble(
//							a -> a.getArcData().getValue() /
//									((RoadInfo) a.getArcData().getInfo()).getWeightAdjustment()).sum();
//					double excess = network_distance - expected_distance;
//					distance_penalty = excess * excess * CANDIDATE_COST_WEIGHT * DISTANCE_PENALTY_FACTOR;
//				}

				// Update cost
				cost += deviation_penalty + distance_penalty;

//				// For debugging only
//				Logger.info("i=" + i + "/" + candidates.size() + " p0: " + p.get(1).getId() + " dummy: " +
//						costDummy + " costNet: " + costNet + " distNet: " + network_distance + " deviation_penalty: " +
//						deviation_penalty + " distance_penalty: " + distance_penalty);

				// Store path
				WeightedPathToCandidate<I> wp = new WeightedPathToCandidate<>(p, cost, cm);
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
		matchesArcCount = new LinkedList<>();
		matchesArcCount.add(0);

		for (int i = candidates.size() - 2; i >= 0; i--) {
			for (WeightedPathToCandidate<I> myWP : allWPs.get(i)) {
				if (myWP.getTarget() == lastNode) {
					List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> subPath = myWP.getPath();
					subPath.removeFirst(); // remove dummy node
					LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> subPathArcs = getPathArcs(subPath);
					pathArcs.addAll(0, subPathArcs);
					path.addAll(0, subPath);
					lastNode = subPath.getFirst();
					matches.addFirst(lastNode);
					matchesArcCount.addFirst(subPathArcs.size());
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

	public LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> getPathArcs(List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> pathNodes) {
		LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> pathArcs = new LinkedList<>();
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> u = null;
		for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> v : pathNodes) {
			if (u != null) {
				// add arc uv to list
				DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> uv = g.getArc(u, v);
				if (uv != null && u.getNodeData() != null && !u.getNodeData().equals(v.getNodeData()))
					pathArcs.add(uv);
			}
			u = v;
		}
		return pathArcs;
	}

	/**
	 * Method to find for each track point the k nearest segments of the graph and
	 * the closest point on each of those segments
	 * 
	 * @param segments  a quadtree with the segments of the network
	 * @param gps_track a list with the track points
	 * @return for each track point a list of candidate matches
	 */
	private ArrayList<LinkedList<CandidateMatch<I>>> getBestKCandidatesForEachTrackPoint(STRtree segments, ArrayList<Point2D> gps_track) {

		// Track size
		int track_size = gps_track.size();

		// Create candidates
		ArrayList<LinkedList<CandidateMatch<I>>> allCandidates = new ArrayList<>(track_size);
		for (Point2D gps_point : gps_track) {
			LinkedList<CandidateMatch<I>> candidates = getBestKCandidatesForTrackPoint(segments, gps_point);
			allCandidates.add(candidates);
		}
//		Logger.debug("Preliminary candidates: " + allCandidates.stream().map(LinkedList::size).reduce(0, Integer::sum));

		// Radius2
		double radius2 = RADIUS * RADIUS;


		// Setup sorter
		Comparator<CandidateMatch<I>> myCandidateComp = (a, b) -> {
			if (a.getMapPoint().distance(a.getGpsPoint()) < b.getMapPoint().distance(b.getGpsPoint()))
				return -1;
			if (a.getMapPoint().distance(a.getGpsPoint()) > b.getMapPoint().distance(b.getGpsPoint()))
				return 1;
			return POINT_LEX_ORDER.compare(a.getMapPoint(), b.getMapPoint());
		};

		// Create new adopted candidates
		Map<Integer, LinkedList<CandidateMatch<I>>> allAdoptedCandidates = new HashMap<>();
		for(int i = 0 ; i < track_size ; i++) {
//			Logger.debug("Running adoption for track point " + i + " / " + track_size);
			for (CandidateMatch<I> candidate : allCandidates.get(i)) {
				Point2D map_point = candidate.getMapPoint();
				if (!map_point.equals(candidate.getSegment().getSource().getNodeData()) &&
						!map_point.equals(candidate.getSegment().getTarget().getNodeData())) {
					I segmentInfo = candidate.getSegmentInfo();
					for (int j = 0; j < track_size; j++) {
						if (i != j) {
							if (allCandidates.get(j).stream().anyMatch(c -> c.getSegmentInfo().equals(segmentInfo) &&
									!c.getMapPoint().equals(c.getSegment().getSource().getNodeData()) &&
									!c.getMapPoint().equals(c.getSegment().getTarget().getNodeData()))) {
								Point2D gps_point = gps_track.get(j);
								double vx = gps_point.getX() - map_point.getX();
								double vy = gps_point.getY() - map_point.getY();
								if (vx*vx + vy*vy <= radius2) {
									CandidateMatch<I> cm = new CandidateMatch<>(gps_point, map_point, candidate.getSegment());
									allAdoptedCandidates.computeIfAbsent(j, m -> new LinkedList<>()).add(cm);
								}
							}
						}
					}
				}
			}
		}

		// Sort and remove candidates above limit
		for(int i = 0 ; i < track_size ; i++) {
			LinkedList<CandidateMatch<I>> candidates = allCandidates.get(i);
			LinkedList<CandidateMatch<I>> adoptedCandidates = allAdoptedCandidates.get(i);

			if (adoptedCandidates != null) {

				// Check how many adopted candidates can be added
				int original_candidates = candidates.size();
				int allowed_insertions = MAX_CAND_N - original_candidates;

				// Reduce size of adopted candidates list to be within limit
				if (allowed_insertions <= 0) {
					Logger.warn("Candidates for track point " + i + " exceed allowed amount (" + original_candidates + ")\n" +
							"You have a dense network. Consider increasing candidate allowance with input option '-k' (currently " + MAX_CAND_N + ").");
				} else if (adoptedCandidates.size() > allowed_insertions) {

					// Sort adopted candidates
					adoptedCandidates.sort(myCandidateComp);

					// Remove adopted candidates that exceed limit
					while (adoptedCandidates.size() > allowed_insertions) {
						adoptedCandidates.removeLast();
					}
				}

				// Add adopted candidates
				candidates.addAll(adoptedCandidates);
			}

			// Sort all candidates
			candidates.sort(myCandidateComp);

			// Remove candidates exceeding limit
			while (candidates.size() > MAX_CAND_N) {
				candidates.removeLast();
			}
		}

		// add new candidates
//		Logger.debug("Adopted candidates: " + allAdoptedCandidates.values().stream().map(LinkedList::size).reduce(0, Integer::sum));
//		Logger.debug("Total candidates: " + allCandidates.stream().map(LinkedList::size).reduce(0, Integer::sum));

		return allCandidates;
	}

	/**
	 * 
	 * @param segments
	 * @param gps_point
	 * @return
	 */
	private LinkedList<CandidateMatch<I>> getBestKCandidatesForTrackPoint(STRtree segments, Point2D gps_point) {
		LinkedList<CandidateMatch<I>> candidates = new LinkedList<>();

		// square-shaped search window around gps_point
		Envelope searchEnvelope = new Envelope(gps_point.getX() - RADIUS, gps_point.getX() + RADIUS,
				gps_point.getY() - RADIUS, gps_point.getY() + RADIUS);

		// test each segment in search window
		for (Object o : segments.query(searchEnvelope)) {
			@SuppressWarnings("unchecked")
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> seg = (DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>) o;

			// find point on seg closest to gps_point
			Line2D line = new Line2D.Double(seg.getSource().getNodeData(), seg.getTarget().getNodeData());
			Point2D nearestPoint = Calculations.nearestPointOnSegment(gps_point, line);

			// if point on seg is within search radius, add it to candidate list
			if (Math.hypot(gps_point.getX() - nearestPoint.getX(), gps_point.getY() - nearestPoint.getY()) <= RADIUS) {
				CandidateMatch<I> cm = new CandidateMatch<>(gps_point, nearestPoint, seg);
				candidates.add(cm);
			}
		}

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
			g.getNodes().removeLast();
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
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> s = path.getFirst();
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> t = path.getLast();
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);
		dijkstra.run(s, t);
		return g.toNodeList(dijkstra.getPath(t.getId()));
	}

	public ArrayList<List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> getShortestPathsForChunks() {
		ArrayList<List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> shortestPaths = new ArrayList<>();
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);
		for (ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> chunk : this.getChunks()) {
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> s = chunk.getFirst();
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> t = chunk.getLast();
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

	public LinkedList<Integer> getMatchesArcCount() {
		return matchesArcCount;
	}

	public ArrayList<Point2D> getTrack() {
		return track;
	}

	public LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> getPathArcs() {
		return pathArcs;
	}

}
