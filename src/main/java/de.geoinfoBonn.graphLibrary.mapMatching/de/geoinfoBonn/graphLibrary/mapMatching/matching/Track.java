package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.tinylog.Logger;
import java.util.stream.IntStream;

import org.geotools.api.data.*;
import org.geotools.api.feature.type.PropertyDescriptor;
import org.geotools.api.filter.Filter;
import org.geotools.api.referencing.ReferenceIdentifier;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.store.ReprojectingFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.geopkg.FeatureEntry;
import org.geotools.geopkg.GeoPackage;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.impl.CoordinateArraySequence;
import org.geotools.api.referencing.crs.CoordinateReferenceSystem;
import org.geotools.api.feature.simple.SimpleFeature;
import org.geotools.api.feature.simple.SimpleFeatureType;

import javax.sql.DataSource;

public class Track {

	private final long id;
	private final int subtrack;
	private final Integer section;
	private final Long type;
	private final ArrayList<Point2D> trackPoints;

	public Track(long id, int subtrack, ArrayList<Point2D> trackPoints) {
		this.id = id;
		this.subtrack = subtrack;
		this.section = null;
		this.type = null;
		this.trackPoints = trackPoints;
	}

	public Track(long id, int subtrack, int section, ArrayList<Point2D> trackPoints) {
		this.id = id;
		this.subtrack = subtrack;
		this.section = section;
		this.type = null;
		this.trackPoints = trackPoints;
	}

	public Track(long id, int subtrack, int section, Long type, ArrayList<Point2D> trackPoints) {
		this.id = id;
		this.subtrack = subtrack;
		this.section = section;
		this.type = type;
		this.trackPoints = trackPoints;
	}

	public ArrayList<Point2D> getTrackPoints() {
		return trackPoints;
	}

	public long getId() {
		return id;
	}

	public int getSubtrack() {
		return subtrack;
	}

	public static ArrayList<Track> importFromShapefile(String filename) {
		return importFromShapefile(filename, null);
	}

	public static ArrayList<Track> importFromShapefile(String filename, String nameColumn) {
		ArrayList<Track> trajectories = new ArrayList<>();
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

				long i = 0;
				try (FeatureIterator<SimpleFeature> features = myFeatureCollection.features()) {
					while (features.hasNext()) {
						SimpleFeature myFeature = features.next();
						long id = nameColumn != null ? (long) myFeature.getAttribute(nameColumn) : i;
						Geometry myGeometry = (Geometry) myFeature.getDefaultGeometry();
						if (myGeometry.getGeometryType().equals("LineString")) {
							ArrayList<Point2D> points = new ArrayList<>();
							Coordinate[] xyz = myGeometry.getCoordinates();
							for (int j = 0; j < xyz.length; j++) {
								Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
								points.add(p);
							}
							trajectories.add(new Track(id, 0, points));
						} else if (myGeometry.getGeometryType().equals("MultiLineString")) {
							for (int geoIndex = 0; geoIndex < myGeometry.getNumGeometries(); geoIndex++) {
								Geometry myGeometryPart = myGeometry.getGeometryN(geoIndex);
								ArrayList<Point2D> points = new ArrayList<>();
								Coordinate[] xyz = myGeometryPart.getCoordinates();
								for (int j = 0; j < xyz.length; j++) {
									Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
									points.add(p);
								}
								trajectories.add(new Track(id, geoIndex, points));
							}
						} else {
							Logger.warn(
									"The shapefile does not contain a polyline but a " + myGeometry.getGeometryType());
						}
						i++;
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

	public static void writeToGpkg(GeoPackage out, String description, ArrayList<Track> paths,
								   boolean withSectionIds, boolean withTypes, CoordinateReferenceSystem crs) {

		long starttime = System.currentTimeMillis();
		SimpleFeatureType featureType = Track.createFeatureType(description, withSectionIds, withTypes, crs);
		SimpleFeatureCollection collection = Track.createFeatureCollection(featureType, paths, withSectionIds, withTypes);
		if (collection.isEmpty()) {
			Logger.warn("Empty collection! Skipping " + description + "!");
			return;
		}
		collection = forceXY(collection);

		try {

			FeatureEntry entry = out.feature(description);
			if(entry == null) {
				entry = new FeatureEntry();
				entry.setDescription(description);
				entry.setBounds(collection.getBounds());
				out.create(entry, collection.getSchema());
			} else {
				entry.getBounds().expandToInclude(collection.getBounds());
			}

			// Write (new) features
			Transaction tx = new DefaultTransaction();
			try {
				try (SimpleFeatureWriter w = out.writer(entry, true, null, tx);
					 SimpleFeatureIterator it = collection.features()) {
					while (it.hasNext()) {
						SimpleFeature f = it.next();
						SimpleFeature g = w.next();
						g.setAttributes(f.getAttributes());
						for (PropertyDescriptor pd : collection.getSchema().getDescriptors()) {
							/* geopkg spec requires booleans to be stored as SQLite integers this fixes
							 * bug reported by GEOT-5904 */
							String name = pd.getName().getLocalPart();
							if (pd.getType().getBinding() == Boolean.class) {
								int bool = 0;
								if (f.getAttribute(name) != null) {
									bool = (Boolean) (f.getAttribute(name)) ? 1 : 0;
								}
								g.setAttribute(name, bool);
							}
						}
						w.write();
					}
				}
				tx.commit();
			} catch (Exception ex) {
				tx.rollback();
				throw new IOException(ex);
			} finally {
				tx.close();
			}

		} catch (IOException e) {
			out.close();
			throw new RuntimeException(e);
		}

		Logger.info(
				description + " written in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
	}

	static SimpleFeatureCollection forceXY(SimpleFeatureCollection fc) {
		CoordinateReferenceSystem sourceCRS = fc.getSchema().getCoordinateReferenceSystem();
		if ((CRS.getAxisOrder(sourceCRS) == CRS.AxisOrder.EAST_NORTH)
				|| (CRS.getAxisOrder(sourceCRS) == CRS.AxisOrder.INAPPLICABLE)) {
			return fc;
		}

		for (ReferenceIdentifier identifier : sourceCRS.getIdentifiers()) {
			try {
				String _identifier = identifier.toString();
				CoordinateReferenceSystem flippedCRS = CRS.decode(_identifier, true);
				if (CRS.getAxisOrder(flippedCRS) == CRS.AxisOrder.EAST_NORTH) {
					ReprojectingFeatureCollection result = new ReprojectingFeatureCollection(fc, flippedCRS);
					return result;
				}
			} catch (Exception e) {
				// couldn't flip - try again
			}
		}
		return fc;
	}



	private static ListFeatureCollection createFeatureCollection(SimpleFeatureType featureType,
																 ArrayList<Track> paths, boolean withSubtrackIds, boolean withTypes) {
		long starttime = System.currentTimeMillis();


		List<SimpleFeature> featureList = new ArrayList<>();
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory(null);
		SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(featureType);

		for (Track path : paths) {
//			if (path.length() < 100) {
//				//logger.warning("track " + path + " is shorter than 100m, skipping export");
//				continue;
//			}
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
			if (withSubtrackIds) {
				featureBuilder.add(path.getSection());
			}
			if (withTypes) {
				featureBuilder.add(path.getType());
			}

			featureList.add(featureBuilder.buildFeature(null));
		}

		Logger.debug("FeatureCollection created in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
		return new ListFeatureCollection(featureType, featureList);
	}

	private static SimpleFeatureType createFeatureType(String description, boolean withSubtrackId, boolean withType,
													   CoordinateReferenceSystem crs) {

		SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();
		builder.setName(description);
		if (crs == null) {
			Logger.warn("Printing tracks with default coordinate reference system WGS84");
			builder.setCRS(DefaultGeographicCRS.WGS84); // <- Default Coordinate reference system
		} else {
			builder.setCRS(crs); // <- Coordinate reference system
		}

		// add attributes in order
		builder.add("the_geom", LineString.class);
		builder.add("track", Long.class);
		builder.add("sid", Integer.class);
		if (withSubtrackId)
			builder.add("subtrack", Integer.class);
		if (withType)
			builder.add("type", Long.class);

		// build the type
		return builder.buildFeatureType();
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

	public int getSection() {
		assert section != null;
		return section;
	}

	public Long getType() {
		return type;
	}
}
