package de.geoinfoBonn.graphLibrary.mapMatching.io;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.function.Function;
import java.util.logging.Logger;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

public class RoadReader<I> {

	private static final Logger LOGGER = Logger.getLogger(RoadReader.class.getName());

	/**
	 * method for reading a shapefile with roads
	 * 
	 * @param filename:      the name of the shapefile
	 * @param typeColName:   the name of the column containing the road type
	 * @param weightColName: the name of the column containing the weight data
	 * @return a list of roads
	 * @throws IOException
	 */
	public static <I> LinkedList<Road<I>> importFromShapefile(String filename, Function<SimpleFeature, I> infoGenerator,
			String weightColName) {
		LinkedList<Road<I>> roads = new LinkedList<Road<I>>();
		try {
			if (filename.endsWith(".shp")) {
				File shpfile = new File(filename);
				Map<String, Object> map = new HashMap<>();
				map.put("url", shpfile.toURI().toURL());

				DataStore dataStore = DataStoreFinder.getDataStore(map);
				String typeName = dataStore.getTypeNames()[0];

				FeatureSource<SimpleFeatureType, SimpleFeature> source = dataStore.getFeatureSource(typeName);
				Filter filter = Filter.INCLUDE; // ECQL.toFilter("BBOX(THE_GEOM, 10,20,30,40)")

				FeatureCollection<SimpleFeatureType, SimpleFeature> collection = source.getFeatures(filter);

				try (FeatureIterator<SimpleFeature> features = collection.features()) {
					while (features.hasNext()) {
						SimpleFeature feature = features.next();

						I type = null;
						if (infoGenerator != null) {
							type = infoGenerator.apply(feature);
//							type = (I) feature.getAttribute(typeColName);
						}
						Geometry myGeometry = (Geometry) feature.getDefaultGeometry();
						double totalLength = myGeometry.getLength();
						double weight = totalLength;
						if (weightColName != null) {
							weight = (double) feature.getAttribute(weightColName);
						}

						if (myGeometry.getGeometryType().equals("LineString")) {
							ArrayList<Point2D> points = new ArrayList<Point2D>();
							Coordinate[] xyz = myGeometry.getCoordinates();
							for (int j = 0; j < xyz.length; j++) {
								Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
								points.add(p);
							}
							Road<I> r = new Road<I>(points, type, weight);
							roads.add(r);
						} else if (myGeometry.getGeometryType().equals("MultiLineString")) {
							for (int geoIndex = 0; geoIndex < myGeometry.getNumGeometries(); geoIndex++) {
								Geometry myGeometryPart = myGeometry.getGeometryN(geoIndex);
								double partLength = myGeometryPart.getLength();
								ArrayList<Point2D> points = new ArrayList<Point2D>();
								Coordinate[] xyz = myGeometryPart.getCoordinates();
								for (int j = 0; j < xyz.length; j++) {
									Point2D p = new Point2D.Double(xyz[j].x, xyz[j].y);
									points.add(p);
								}
								Road<I> r = new Road<I>(points, type, weight * partLength / totalLength);
								roads.add(r);
							}
						} else {
							LOGGER.warning("Feature is not a polyline but a " + myGeometry.getGeometryType());
						}
					}
				} finally {
					dataStore.dispose();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return roads;
	}

	public static CoordinateReferenceSystem readCRS(String filename) throws FileNotFoundException {
		File prj = new File(filename);
		if (!prj.exists() || !filename.endsWith(".prj")) {
			throw new FileNotFoundException("Please provide filename of .prj file! (tried to read: " + filename + ")");
		}

		CoordinateReferenceSystem crs = null;
		try {
			crs = CRS.parseWKT(Files.readString(prj.toPath()));
			return crs;
		} catch (FactoryException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
}
