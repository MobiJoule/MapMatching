package de.geoinfoBonn.graphLibrary.mapMatching.matching.types;

public enum EdgeType implements Typed {
	ROAD, TESSA, HOLE, CANDIDATE, BOUNDARY, UNMATCHED;

	@Override
	public EdgeType getType() {
		return this;
	}
}
