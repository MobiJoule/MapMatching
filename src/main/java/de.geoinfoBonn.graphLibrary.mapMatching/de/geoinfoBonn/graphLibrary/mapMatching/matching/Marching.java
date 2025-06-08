package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.function.Function;
import java.util.logging.Logger;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.index.strtree.STRtree;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.core.geometry.Calculations;
import de.geoinfoBonn.graphLibrary.core.geometry.PointComparator;
import de.geoinfoBonn.graphLibrary.core.shortestPath.Dijkstra;
import de.geoinfoBonn.graphLibrary.core.shortestPath.MultiTargetNodeVisitor;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.DummyNode;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.EdgeType;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.InterDummyNode;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.Typed;

public class Marching<I extends Typed> {
	public static final Logger logger = Logger.getLogger(MethodHandles.lookup().lookupClass().getName());

	private DiGraph<Point2D, DoubleWeightDataWithInfo<I>> g;
	private STRtree segments;

	private LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> path;
	private LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> pathArcs;
	private LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> matches;
	private ArrayList<LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>>> matchSegments;
	private ArrayList<CandidateMatch<I>> matchedCandidates;
	private ArrayList<Point2D> track;
	private HashSet<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> unmatchedNodes;

	private final int originalNumberOfArcs;
	private final int originalNumberOfNodes;

	public static boolean VERBOSE = false;
	private static final PointComparator POINT_LEX_ORDER = new PointComparator();

	private Variant variant = Variant.UNMATCHED;

	private double radius = 25.0;
	private int maxCandN = Integer.MAX_VALUE;
	private double candidateCostWeight = 0.01;
	private double unmatchedCostWeight = 2.5;
	private boolean addUnmatchedCandidates = true;
	private double tessalationCostWeight = 1.2;
	private double obstacleBoundaryCostWeight = 4.0;

	private Function<EdgeType, I> infoGenerator;

	/**
	 * 
	 * @param inputGraph    graph representation of the network to match to
	 * @param infoGenerator function to map between EdgeType and the generic type
	 *                      parameter to add temporary arcs to the graph
	 */
	public Marching(DiGraph<Point2D, DoubleWeightDataWithInfo<I>> inputGraph, Function<EdgeType, I> infoGenerator) {
		this.infoGenerator = infoGenerator;
		this.initializeInputGraph(inputGraph);
		this.originalNumberOfArcs = g.getArcs().size();
		this.originalNumberOfNodes = g.getNodes().size();
		logger.fine("Size of original input graph :" + originalNumberOfNodes + " " + originalNumberOfArcs);
	}

	public void match(ArrayList<Point2D> inputTrack) {
		long startTime = System.currentTimeMillis();

		// find for each gps point the k best candidates
		ArrayList<LinkedList<CandidateMatch<I>>> candidates = getBestKCandidatesForEachTrackPoint(segments, inputTrack,
				radius, maxCandN);

		// add off-road candidates to candidate lists
		track = null;
		unmatchedNodes = new HashSet<>();
		if (variant == Variant.UNMATCHED || variant == Variant.TESSELLATION) {
			track = inputTrack;
			for (int i = 0; i < inputTrack.size(); i++) {
				CandidateMatch<I> unmatchedCandidate = new CandidateMatch<>(inputTrack.get(i), inputTrack.get(i), null);
				candidates.get(i).add(unmatchedCandidate);
				DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> unmatchedNode = g.addNode(inputTrack.get(i));
				unmatchedCandidate.setNode(unmatchedNode);
				unmatchedNodes.add(unmatchedNode);
			}
		} else {
			// reduce track to gps points that have at least one candidate
			track = new ArrayList<Point2D>();
			ArrayList<LinkedList<CandidateMatch<I>>> candidates_new = new ArrayList<>();
			for (int i = 0; i < inputTrack.size(); i++) {
				if (candidates.get(i).size() > 0) {
					track.add(inputTrack.get(i));
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
					LinkedList<CandidateMatch<I>> pointsOnSeg = candidatesToInject.get(c.mapSegment);
					if (pointsOnSeg == null) {
						pointsOnSeg = new LinkedList<CandidateMatch<I>>();
						candidatesToInject.put(c.mapSegment, pointsOnSeg);
					}
					pointsOnSeg.add(c);
				}
			}
		}

		// add candidate points as new graph nodes and arcs connecting them
		for (Entry<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>, LinkedList<CandidateMatch<I>>> e : candidatesToInject
				.entrySet()) {
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> arc1 = e.getKey();
			LinkedList<CandidateMatch<I>> candsOnArc1 = e.getValue();

			if (POINT_LEX_ORDER.compare(arc1.getSource().getNodeData(), arc1.getTarget().getNodeData()) == -1) {
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
				double scaledWeight = arc1.getArcData().getValue() * newArcLength / oldArcLength;
				DoubleWeightDataWithInfo<I> dw1 = new DoubleWeightDataWithInfo<>(scaledWeight,
						arc1.getArcData().getInfo());
				g.addArc(newNodes.get(i - 1), newNodes.get(i), dw1);

				if (arc1.getTwin() != null) {
					scaledWeight = arc1.getTwin().getArcData().getValue() * newArcLength / oldArcLength;
					DoubleWeightDataWithInfo<I> dw2 = new DoubleWeightDataWithInfo<>(scaledWeight,
							arc1.getTwin().getArcData().getInfo());
					g.addArc(newNodes.get(i), newNodes.get(i - 1), dw2);
				}
			}
		}

		logger.fine("Size of graph after node insertion: " + g.getNodes().size() + " " + g.getArcs().size());

		// add off-road edges
		if (addUnmatchedCandidates) {
			for (int i = 0; i < inputTrack.size(); i++) {
				CandidateMatch<I> c1 = candidates.get(i).getLast(); // fetch off-road candidate for track point i

				if (i < inputTrack.size() - 1) {
					for (CandidateMatch<I> c2 : candidates.get(i + 1)) {
						addUnmatchedSegment(c1, c2);
					}
				}
				if (i > 0) {
					for (CandidateMatch<I> c0 : candidates.get(i - 1)) {
						if (c0.mapSegment != null) {
							addUnmatchedSegment(c0, c1);
						}
					}
				}
			}
			logger.fine("Size of graph with off-road edges: " + g.getNodes().size() + " " + g.getArcs().size());
		}

		// add dummy nodes and arcs
		ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> dummies = new ArrayList<>();
		ArrayList<HashMap<CandidateMatch<I>, DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>>> intermediateDummies = new ArrayList<>();
		for (int i = 0; i < candidates.size() - 1; i++) {
			LinkedList<CandidateMatch<I>> candidatesForCurrentPoint = candidates.get(i);
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> dummy = g.addNode(new DummyNode());

			HashMap<CandidateMatch<I>, DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> currIntermediateDummies = new HashMap<>();
			for (CandidateMatch<I> cm : candidatesForCurrentPoint) {
				DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> intermediateDummy = g.addNode(new InterDummyNode());
				currIntermediateDummies.put(cm, intermediateDummy);
				// is later filled with distance of path until here
				g.addArc(dummy, intermediateDummy,
						new DoubleWeightDataWithInfo<>(0, infoGenerator.apply(EdgeType.ROAD)));
				g.addArc(intermediateDummy, cm.getNode(),
						new DoubleWeightDataWithInfo<>(cm.getCandidateCost(), infoGenerator.apply(EdgeType.CANDIDATE)));
			}

			dummies.add(dummy);
			intermediateDummies.add(currIntermediateDummies);
		}

		// initialize dijkstra instance
		Dijkstra<Point2D, DoubleWeightDataWithInfo<I>> dijkstra = new Dijkstra<>(g);
		dijkstra.setDistanceComputer(new MarchingDistanceComputer<>(candidateCostWeight, unmatchedCostWeight,
				tessalationCostWeight, obstacleBoundaryCostWeight));

		// search paths
		ArrayList<LinkedList<WeightedPathToCandidate<I>>> allWPs = new ArrayList<>();
		for (int i = 0; i < candidates.size() - 1; i++) {

			// define source
			DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> dummy = dummies.get(i);

			// update weights of dummy arcs
			if (i > 0) {
				for (WeightedPathToCandidate<I> myPath : allWPs.get(i - 1)) {
					DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> intermediateDummy = intermediateDummies.get(i)
							.get(myPath.getTargetCandidate());
					DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> dummyArc = g.getArc(dummy, intermediateDummy);
					if (dummyArc != null) {
						dummyArc.setArcData(new DoubleWeightDataWithInfo<>(myPath.getDistance(),
								infoGenerator.apply(EdgeType.ROAD)));
					} else {
						logger.severe("Dummy arc is null!");
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
				List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> p = dijkstra.getPath(cm.getNode());
				List<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> arcs = dijkstra.getPathArcs(cm.getNode());

				double d = dijkstra.getDistance(cm.getNode());
				WeightedPathToCandidate<I> wp = new WeightedPathToCandidate<>(p, d, cm, arcs);
				myPathList.add(wp);
			}
			allWPs.add(myPathList);
		}

		// find last candidate of optimal solution
		double minTotalCost = Double.POSITIVE_INFINITY;
		CandidateMatch<I> bestCandidate = null;
		LinkedList<WeightedPathToCandidate<I>> myPathList = allWPs.get(candidates.size() - 2);
		for (WeightedPathToCandidate<I> myPath : myPathList) {
			double totalCost = myPath.getDistance()
					+ candidateCostWeight * myPath.getTargetCandidate().getCandidateCost();
			if (totalCost < minTotalCost) {
				minTotalCost = totalCost;
				bestCandidate = myPath.getTargetCandidate();
			}
		}

		// reconstruct solution
		DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> lastNode = bestCandidate.getNode();
		path = new LinkedList<>();
		pathArcs = new LinkedList<>();
		matches = new LinkedList<>();
		matches.add(lastNode);
		matchSegments = new ArrayList<>();
		matchedCandidates = new ArrayList<>();

		matchedCandidates.add(bestCandidate);

		for (int i = candidates.size() - 2; i >= 0; i--) {
			for (WeightedPathToCandidate<I> myWP : allWPs.get(i)) {
				if (myWP.getTarget() == lastNode) {
					List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> subPath = myWP.getPath();
					subPath.remove(0); // remove dummy node
					subPath.remove(0); // remove intermediateDummy
					var currPathArcs = g.getPathArcs(subPath);
					pathArcs.addAll(0, currPathArcs);
					path.addAll(0, subPath);
					matchSegments.add(0, currPathArcs);
					lastNode = subPath.get(0);
					matches.add(0, lastNode);
					matchedCandidates.add(0, myWP.getTargetCandidate());
					break;
				}
			}
		}

		logger.info("total cost: " + minTotalCost + ", running time: "
				+ (System.currentTimeMillis() - startTime) / 1000.0 + "s");
	}

	private void addUnmatchedSegment(CandidateMatch<I> cFrom, CandidateMatch<I> cTo) {
		Point2D pFrom = cFrom.getMapPoint();
		Point2D pTo = cTo.getMapPoint();

		g.addArc(cFrom.getNode(), cTo.getNode(),
				new DoubleWeightDataWithInfo<>(pFrom.distance(pTo), infoGenerator.apply(EdgeType.UNMATCHED)));
	}

	public LineString getSegment(Point2D p, Point2D q, GeometryFactory gf) {
		return gf.createLineString(
				new Coordinate[] { new Coordinate(p.getX(), p.getY()), new Coordinate(q.getX(), q.getY()) });
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
		Collections.sort(candidates, new DistanceComperator<>());

		// remove all but the best k candidates
		while (candidates.size() > maxCandN)
			candidates.removeLast();
		return candidates;
	}

	private void initializeInputGraph(DiGraph<Point2D, DoubleWeightDataWithInfo<I>> input_graph) {
		g = input_graph;

		// compute strtree with segments of graph (for double-edges: only arcs pointing
		// from left to right)
		segments = new STRtree();
		for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> s : g.getArcs()) {
			if (s.getTwin() == null
					|| POINT_LEX_ORDER.compare(s.getSource().getNodeData(), s.getTarget().getNodeData()) == -1) {
				Envelope env = new Envelope();
				env.expandToInclude(s.getSource().getNodeData().getX(), s.getSource().getNodeData().getY());
				env.expandToInclude(s.getTarget().getNodeData().getX(), s.getTarget().getNodeData().getY());
				segments.insert(env, s);
			}
		}
	}

	public void setRadius(double radius) {
		this.radius = radius;
	}

	public void setMaxCandN(int maxCandN) {
		this.maxCandN = maxCandN;
	}

	public void setCandidateCostWeight(double candidateCostWeight) {
		this.candidateCostWeight = candidateCostWeight;
	}

	public void setVariant(Variant variant) {
		this.variant = variant;
	}

	public void setUnmatchedCostWeight(double unmatchedCostWeight) {
		this.unmatchedCostWeight = unmatchedCostWeight;
	}

	public void setObstacleBoundaryCostWeight(double obstacleBoundaryCostWeight) {
		this.obstacleBoundaryCostWeight = obstacleBoundaryCostWeight;
	}

	public void setTessalationCostWeight(double tessalationCostWeight) {
		this.tessalationCostWeight = tessalationCostWeight;
	}

	/**
	 * Method to restore the graph by deleting all new nodes and their incident
	 * edges that were added for map matching
	 */
	public void restoreGraph() {
		for (int n = g.getArcs().size(); n > originalNumberOfArcs; n = g.getArcs().size()) {
			DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> a = g.getArcs().get(n - 1);
			a.getSource().getOutgoingArcs().remove(a);
			a.getTarget().getIncomingArcs().remove(a);
			g.getArcs().remove(n - 1);
		}

		while (g.getNodes().size() > originalNumberOfNodes) {
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
			if (!unmatchedNodes.contains(node)) {
				currentChunk.add(node);
			} else {
				if (currentChunk.size() > 0) {
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

	public LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> idsToPath(List<Integer> ids) {
		LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> subpath = new LinkedList<>();
		for (Integer i : ids)
			for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>> arc : matchSegments.get(i))
				subpath.add(arc.getSource());
		return subpath;
	}

	public ArrayList<Point2D> idsToTracks(List<Integer> ids) {
		ArrayList<Point2D> subtrack = new ArrayList<>();
		for (Integer i : ids)
			subtrack.add(track.get(i));
		return subtrack;
	}

	public LinkedList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> getMatches() {
		return matches;
	}

	public ArrayList<CandidateMatch<I>> getMatchedCandidates() {
		return matchedCandidates;
	}

	public ArrayList<Point2D> getTrack() {
		return track;
	}

	public LinkedList<DiGraphArc<Point2D, DoubleWeightDataWithInfo<I>>> getPathArcs() {
		return pathArcs;
	}

	public enum Variant {
		STRICT, // no unmatched candidates
		UNMATCHED, // add unmatched candidates
		TESSELLATION // add unmatched candidates + better tesselation costs
	}
}
