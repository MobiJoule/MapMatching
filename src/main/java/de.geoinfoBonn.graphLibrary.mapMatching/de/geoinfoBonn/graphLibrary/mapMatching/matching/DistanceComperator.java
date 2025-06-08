package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.util.Comparator;

import de.geoinfoBonn.graphLibrary.core.geometry.PointComparator;

public class DistanceComperator<I> implements Comparator<CandidateMatch<I>> {

	private static final PointComparator POINT_LEX_ORDER = new PointComparator();

	@Override
	public int compare(CandidateMatch<I> a, CandidateMatch<I> b) {
		if (a.getMapPoint().distance(a.getGpsPoint()) < b.getMapPoint().distance(b.getGpsPoint()))
			return -1;
		if (a.getMapPoint().distance(a.getGpsPoint()) > b.getMapPoint().distance(b.getGpsPoint()))
			return 1;
		return POINT_LEX_ORDER.compare(a.getMapPoint(), b.getMapPoint());
	}
}
