package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.data.shapefile.ShapefileDumper;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

public class Track {

	private static final Logger logger = Logger.getLogger(Track.class.getName());

	private String id;
	private int subtrack;
	private ArrayList<Point2D> trackPoints;

	public Track(String id, int subtrack, ArrayList<Point2D> trackPoints) {
		this.id = id;
		this.subtrack = subtrack;
		this.trackPoints = trackPoints;
	}

	public ArrayList<Point2D> getTrackPoints() {
		return trackPoints;
	}

	public String getId() {
		return id;
	}

	public int getSubtrack() {
		return subtrack;
	}

	public static ArrayList<Track> importFromShapefile(String filename) {
		return importFromShapefile(filename, null);
	}

	public static ArrayList<Track> importFromShapefile(String filename, String nameColumn) {
		ArrayList<Track> trajectories = new ArrayList<Track>();
		try {
			if (filename.endsWith(".shp")) {
				File shpfile = new File(filename);
				Map<String, Object> map = new HashMap<>();
				map.put("url", shpfile.toURI().toURL());

				DataStore dataStore = DataStoreFinder.getDataStore(map);
				String typeName = dataStore.getTypeNames()[0];

				FeatureSource<SimpleFeatureType, SimpleFeature> source = dataStore.getFeatureSource(typeName);
				Filter filter = Filter.INCLUDE; // ECQL.toFilter("BBOX(THE_GEOM, 10,20,30,40)")

				FeatureCollection<SimpleFeatureType, SimpleFeature> myFeatureCollection = source.getFeatures(filter);
				try (FeatureIterator<SimpleFeature> features = myFeatureCollection.features()) {
					while (features.hasNext()) {
						SimpleFeature myFeature = features.next();
						String id = nameColumn != null ? (String) myFeature.getAttribute(nameColumn)
								: myFeature.getID();
						Geometry myGeometry = (Geometry) myFeature.getDefaultGeometry();
						if (myGeometry.getGeometryType().equals("LineString")) {
							ArrayList<Point2D> points = new ArrayList<Point2D>();
							Coordinate[] xyz = myGeometry.getCoordinates();
							for (int j = 0; j < xyz.length; j++) {
								Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
								points.add(p);
							}
							trajectories.add(new Track(id, 0, points));
						} else if (myGeometry.getGeometryType().equals("MultiLineString")) {
							for (int geoIndex = 0; geoIndex < myGeometry.getNumGeometries(); geoIndex++) {
								Geometry myGeometryPart = myGeometry.getGeometryN(geoIndex);
								ArrayList<Point2D> points = new ArrayList<Point2D>();
								Coordinate[] xyz = myGeometryPart.getCoordinates();
								for (int j = 0; j < xyz.length; j++) {
									Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
									points.add(p);
								}
								trajectories.add(new Track(id, geoIndex, points));
							}
						} else {
							logger.warning(
									"The shapefile does not contain a polyline but a " + myGeometry.getGeometryType());
						}
					}
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				} finally {
					dataStore.dispose();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return trajectories;

	}

	public static void writeToShapeFile(File outputFile, ArrayList<Track> paths, ArrayList<Integer> subtrackIds,
			ArrayList<String> types, CoordinateReferenceSystem crs) {
		Track.writeToShapeFile(outputFile.getParentFile(), outputFile.getName(), paths, subtrackIds, types, crs);
	}

	public static void writeToShapeFile(File outputDir, String filename, ArrayList<Track> paths,
			ArrayList<Integer> subtrackIds, ArrayList<String> types, CoordinateReferenceSystem crs) {

		long starttime = System.currentTimeMillis();
		SimpleFeatureType featureType = Track.createFeatureType(subtrackIds != null, types != null, crs);
		DefaultFeatureCollection collection = Track.createFeatureCollection(featureType, paths, subtrackIds, types);
		if (collection.isEmpty()) {
			logger.warning("Empty collection! Skipping " + filename + "!");
			return;
		}
		Track.writeToFile(outputDir, filename, featureType, collection);

		logger.info(
				"Tracks written to " + filename + " in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
	}

	private static void writeToFile(File outputDir, String filename, SimpleFeatureType featureType,
			DefaultFeatureCollection collection) {
		try {
			ShapefileDumper dumper = new ShapefileDumper(outputDir);
			// optiona, set a target charset
			dumper.setCharset(Charset.forName("UTF-8"));
			// split when shp or dbf reaches 100MB
			int maxSize = 100 * 1024 * 1024;
			dumper.setMaxDbfSize(maxSize);
			// actually dump data
			dumper.dump(filename, collection);
		} catch (IOException e) {
			logger.severe("Error writing shape.");
			e.printStackTrace();
		}
	}

	private static DefaultFeatureCollection createFeatureCollection(SimpleFeatureType featureType,
			ArrayList<Track> paths, ArrayList<Integer> subtrackIds, ArrayList<String> types) {
		long starttime = System.currentTimeMillis();

		DefaultFeatureCollection collection = new DefaultFeatureCollection();
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory(null);
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		int pathIndex = 0;
		for (Track path : paths) {
			if (path.length() < 100) {
				logger.warning("track " + path + " is shorter than 100m, skipping export");
				continue;
			}
			CoordinateArraySequence seq = new CoordinateArraySequence(path.trackPoints.size());
			int i = 0;
			for (Point2D p : path.trackPoints) {
				seq.setOrdinate(i, 0, p.getX());
				seq.setOrdinate(i, 1, p.getY());
				i++;
			}
			LineString ls = new LineString(seq, gf);

			featureBuilder.add(ls);
			featureBuilder.add(path.getId());
			featureBuilder.add(path.getSubtrack());
			if (subtrackIds != null) {
				featureBuilder.add(subtrackIds.get(pathIndex));
			}
			if (types != null) {
				featureBuilder.add(types.get(pathIndex));
			}

			SimpleFeature feature = featureBuilder.buildFeature(null);
			collection.add(feature);
			pathIndex++;
		}

		logger.finest("FeatureCollection created in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
		return collection;
	}

	private static SimpleFeatureType createFeatureType(boolean withSubtrackId, boolean withType,
			CoordinateReferenceSystem crs) {
		long starttime = System.currentTimeMillis();

		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName("Location");
		if (crs == null)
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		else
			builder.setCRS(crs); // <- Coordinate reference system

		// add attributes in order
		builder.add("the_geom", LineString.class);
		builder.add("track", String.class);
		builder.add("sid", Integer.class);
		if (withSubtrackId)
			builder.add("subtrack", Integer.class);
		if (withType)
			builder.length(64).add("type", String.class); // <- 15 chars width for name field

		// build the type
		final SimpleFeatureType LOCATION = builder.buildFeatureType();

		logger.finest("FeatureType created in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
		return LOCATION;
	}

	public LineString getSegment(Point p, Point q, GeometryFactory gf) {
		return gf.createLineString(
				new Coordinate[] { new Coordinate(p.getX(), p.getY()), new Coordinate(q.getX(), q.getY()) });
	}

	public LineString getSegment(Point2D p, Point2D q, GeometryFactory gf) {
		return gf.createLineString(
				new Coordinate[] { new Coordinate(p.getX(), p.getY()), new Coordinate(q.getX(), q.getY()) });
	}

	public int size() {
		return trackPoints.size();
	}

	public double length() {
		return IntStream.range(0, trackPoints.size() - 1)//
				.mapToDouble(i -> trackPoints.get(i).distance(trackPoints.get(i + 1)))//
				.sum();
	}
}
