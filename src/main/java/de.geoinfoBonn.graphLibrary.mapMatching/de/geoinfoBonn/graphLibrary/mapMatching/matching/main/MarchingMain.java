package de.geoinfoBonn.graphLibrary.mapMatching.matching.main;

import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Optional;
import java.util.logging.Logger;

import org.geotools.api.referencing.crs.CoordinateReferenceSystem;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.mapMatching.io.Road;
import de.geoinfoBonn.graphLibrary.mapMatching.io.RoadReader;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Marching;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Marching.Variant;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Track;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.types.EdgeType;

public class MarchingMain extends AbstractMain {

	private static final Logger logger = Logger.getLogger(MarchingMain.class.getName());

	public static void main(String[] args) {
		if (containsOptionalArg(args, "-h")) {
			printHelpText();
			System.exit(0);
		}

		boolean assertsEnabled = false;
		assert assertsEnabled = true;
		if (assertsEnabled)
			logger.info("Assertions enabled.");

		Variant variant = Variant.UNMATCHED;
		if (containsOptionalArg(args, "-variant"))
			variant = Variant.valueOf(getOptionalArg(args, "-variant").toUpperCase());

		boolean exportSingleFiles = false;
		if (containsOptionalArg(args, "-single"))
			exportSingleFiles = true;

		double radius = 25.0;
		int maxCandN = Integer.MAX_VALUE;
		double candidateCostWeight = 0.02;
		double unmatchedCostWeight = 3.0;
		double obstacleBoundaryCostWeight = 5;
		double tessalationCostWeight = 1.5;

		if (containsOptionalArg(args, "-r"))
			radius = Double.parseDouble(getOptionalArg(args, "-r"));
		if (containsOptionalArg(args, "-k"))
			maxCandN = Integer.parseInt(getOptionalArg(args, "-k"));
		if (containsOptionalArg(args, "-c"))
			candidateCostWeight = Double.parseDouble(getOptionalArg(args, "-c"));
		if (containsOptionalArg(args, "-w"))
			unmatchedCostWeight = Double.parseDouble(getOptionalArg(args, "-w"));
		if (containsOptionalArg(args, "-ob"))
			obstacleBoundaryCostWeight = Double.parseDouble(getOptionalArg(args, "-ob"));
		if (containsOptionalArg(args, "-t"))
			tessalationCostWeight = Double.parseDouble(getOptionalArg(args, "-t"));

		File outputDir = containsOptionalArg(args, "-o") ? new File(getOptionalArg(args, "-o")) : new File("out");
		outputDir.mkdirs();

		logger.fine("#################################################################");
		logger.fine("RADIUS: \t\t\t" + radius);
		logger.fine("MAX_CAND_N: \t\t" + maxCandN);
		logger.fine("CANDIDATE_COST_WEIGHT: \t" + candidateCostWeight);
		logger.fine("OFF_ROAD_WEIGHT: \t\t" + unmatchedCostWeight);
		logger.fine("OBSTACLE_BOUNDARY_COST: \t" + obstacleBoundaryCostWeight);
		logger.fine("TESSALATION_COST: \t" + tessalationCostWeight);
		if (variant == Variant.UNMATCHED)
			logger.fine("Variant: \t\t\twithout tessalation");
		else
			logger.fine("Variant: \t\t\twith tessalation");
		logger.fine("EXPORT SINGLE: \t\t" + exportSingleFiles);
		logger.fine("#################################################################");

		String weightColName = containsOptionalArg(args, "-l") ? getOptionalArg(args, "-l") : null;

		// load network and trajectory data
		LinkedList<Road<EdgeType>> roads = RoadReader.importFromGpkg(args[0], f -> EdgeType.ROAD, weightColName,null);
		DiGraph<Point2D, DoubleWeightDataWithInfo<EdgeType>> g = Road.buildGraph(roads);
		ArrayList<Track> trajectories = Track.importFromShapefile(args[1], "gpsies_id");
		logger.info("Number of trajectories in input file: " + trajectories.size());

		CoordinateReferenceSystem crs = null;
		crs = RoadReader.readCRS(args[1].split("\\.")[0] + ".prj");


		ArrayList<Track> paths = new ArrayList<>();
		ArrayList<Track> matches = new ArrayList<>();
		ArrayList<Track> chunks = new ArrayList<>();
		ArrayList<Track> chunkPaths = new ArrayList<>();
		ArrayList<Track> globalPaths = new ArrayList<>();
		ArrayList<Track> segments = new ArrayList<>();

		Marching<EdgeType> m = new Marching<>(g, e -> e);
		m.setVariant(variant);
		m.setRadius(radius);
		m.setMaxCandN(maxCandN);
		m.setCandidateCostWeight(candidateCostWeight);
		m.setUnmatchedCostWeight(unmatchedCostWeight);
		m.setObstacleBoundaryCostWeight(obstacleBoundaryCostWeight);
		m.setTessalationCostWeight(tessalationCostWeight);

		int counter = 0;
		for (Track track : trajectories) {
			logger.info("Starting track #" + counter++ + "/" + trajectories.size() + " (id=" + track.getId()
					+ ", subtrack=" + track.getSubtrack() + ")");

			m.match(track.getTrackPoints());

			Track matched = new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(m.getPath()));
			paths.add(matched);

			int segmentCounter = 1;
			for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<EdgeType>> arc : m.getPathArcs()) {
				ArrayList<Point2D> segment = new ArrayList<Point2D>();
				segment.add(arc.getSource().getNodeData());
				segment.add(arc.getTarget().getNodeData());
				segments.add(new Track(track.getId(), track.getSubtrack(),segmentCounter++, segment));
			}

			ArrayList<Point2D> track_matchedPoints = extractPointsFromPath(m.getMatches());
			for (int i = 0; i < track.getTrackPoints().size(); i++) {
				ArrayList<Point2D> match = new ArrayList<Point2D>();
				match.add(track.getTrackPoints().get(i));
				match.add(track_matchedPoints.get(i));
				matches.add(new Track(track.getId(), track.getSubtrack(), match));
			}

			int chunkCounter = 1;
			for (ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<EdgeType>>> chunk : m.getChunks()) {
				chunks.add(new Track(track.getId(), track.getSubtrack(),chunkCounter++, extractPointsFromPath(chunk)));
			}
			for (List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<EdgeType>>> chunkPath : m
					.getShortestPathsForChunks()) {
				chunkPaths.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(chunkPath)));
			}

			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<EdgeType>>> p = m.getShortestPathForWholeTrajectory();
			if (p != null && p.size() > 1) {
				ArrayList<Point2D> p_points = extractPointsFromPath(p);
				globalPaths.add(new Track(track.getId(), track.getSubtrack(), p_points));
			}

			m.restoreGraph();
			logger.info("Track done!");

			if (!exportSingleFiles) {
				exportResults(args, outputDir, Optional.of(track.getId() + "_"), Optional.empty(), crs, paths, matches,
						chunks, chunkPaths, globalPaths, segments);

				paths = new ArrayList<>();
				matches = new ArrayList<>();
				chunks = new ArrayList<>();
				chunkPaths = new ArrayList<>();
				globalPaths = new ArrayList<>();
				segments = new ArrayList<>();
			}

		}

		logger.info("Writing results to shapefiles...");

		if (exportSingleFiles)
			exportResults(args, outputDir, Optional.empty(), Optional.empty(), crs, paths, matches, chunks,
					chunkPaths, globalPaths, segments);

		logger.info("Finished.");
	}

	public static void exportResults(String[] args, File outputDir, Optional<String> prefix, Optional<String> suffix,
									 CoordinateReferenceSystem crs, ArrayList<Track> paths, ArrayList<Track> matches,
									 ArrayList<Track> chunks, ArrayList<Track> chunkPaths,
									 ArrayList<Track> globalPaths, ArrayList<Track> segments) {
//		Track.writeToGpkg(outputDir, prefix.orElse("") + "paths" + suffix.orElse(""), paths, false, false, crs);
//
//		// Optional outputs
//		if (containsOptionalArg(args, "-ma"))
//			Track.writeToGpkg(outputDir, prefix.orElse("") + "matches" + suffix.orElse(""), matches, false, false,
//					crs);
//
//		if (containsOptionalArg(args, "-ch")) // each section of the input that could be matched onto roads
//			Track.writeToGpkg(outputDir, prefix.orElse("") + "chunks" + suffix.orElse(""), chunks, true, false,
//					crs);
//
//		if (containsOptionalArg(args, "-chp"))
//			Track.writeToGpkg(outputDir, prefix.orElse("") + "chunkPaths" + suffix.orElse(""), chunkPaths, true, false, crs);
//
//		if (containsOptionalArg(args, "-gp"))
//			Track.writeToGpkg(outputDir, prefix.orElse("") + "globalPaths" + suffix.orElse(""), globalPaths, false, false, crs);
//
//		if (containsOptionalArg(args, "-s"))
//			Track.writeToGpkg(outputDir, prefix.orElse("") + "segments" + suffix.orElse(""), segments, false,  false, crs);
	}

	public static ArrayList<Point2D> extractPointsFromPath(
			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<EdgeType>>> path) {
		ArrayList<Point2D> pointList = new ArrayList<Point2D>();
		for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<EdgeType>> node : path) {
			pointList.add(node.getNodeData());
		}
		return pointList;
	}

	private static void printHelpText() {
		System.out.println("Executable for map marching.");
		System.out.println();
		System.out.println("Mandatory arguments:");
		System.out.println("1) path to road shapefile (must be program first argument) [string] ");
		System.out.println("2) path to trajectory shapefile (must be second program argument) [string] ");
		System.out.println();
		System.out.println("Optional arguments:");
		System.out.println("-h  , flag to (just) print help");
		System.out.println("-r  , parameter for radius (default: 25.0) [double]");
		System.out.println("-k  , parameter for maximum number of candidates (default: INF) [int]");
		System.out.println("-c  , parameter for candidate cost weight (default: 0.01) [double]");
		System.out.println("-w  , parameter for offroad weight (default: 1.5) [double]");
		System.out.println("-v  , flag to print verbose output");
		System.out.println("-l  , name of weight column in road data [string]");
		System.out.println("-t  , name of type column in road data [string]");
		System.out.println("-ma , output path for matches [string]");
		System.out.println("-ch , output path for chunks [string]");
		System.out.println("-chp, output path for chunk paths [string]");
		System.out.println("-gp , output path for global paths [string]");
		System.out.println("-s  , output path for segments [string]");
		System.out.println("-single, flag to write all outputs to the same file, otherwise one file per track");
		System.out.println();
		System.out.println("Developed by Prof. Dr.-Ing. Jan-Henrik Haunert and Axel Forsch");
		System.out.println("Contact: forsch@igg.uni-bonn.de");
	}
}
