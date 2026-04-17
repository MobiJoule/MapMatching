package de.geoinfoBonn.graphLibrary.mapMatching.core.generic.data;

import java.util.ArrayList;
import java.util.List;

/**
 * Models arc data which can have an arbitrary amount of weights.
 * 
 * @author Axel Forsch
 *
 */
public class MultiWeightData {
	protected ArrayList<Long> w;

	/**
	 * Creates a MultiWeightData object with empty weight list.
	 */
	public MultiWeightData() {
		this.w = new ArrayList<>();
	}

	/**
	 * Creates a MultiWeightData object with given weights.
	 * 
	 * @param weights list of weights
	 */
	public MultiWeightData(ArrayList<Long> weights) {
		this.w = weights;
	}

	/**
	 * Creates a (deep) copy of the given MultiWeightData
	 * 
	 * @param lwd data to copy
	 */
	public MultiWeightData(MultiWeightData lwd) {
		this.w = new ArrayList<>();
		this.w.addAll(lwd.w);
	}

	public long getWeight(int index) {
		return w.get(index);
	}

	public List<Long> getWeights() {
		return w;
	}

	public void removeWeight(int index) {
		w.remove(index);
	}

	public int addWeight(long value) {
		w.add(value);
		return w.size() - 1;
	}

	public void setWeight(int index, long value) {
		if (w.size() > index) {
			w.set(index, value);
		} else {
			throw new IllegalArgumentException("Index (" + index + ") exceeds number of weights (" + w.size() + ").");
		}
	}

	public int numberOfWeights() {
		return w.size();
	}

	public ArrayList<MultiWeightData> splitData() {

		ArrayList<Long> w1 = new ArrayList<>(w.size());
		ArrayList<Long> w2 = new ArrayList<>(w.size());
		for (int i = 0; i < w.size(); i++) {
			long currW = w.get(i);
			if (currW % 2 != 0) {
				w1.add((currW + 1) / 2);
				w2.add((currW - 1) / 2);
			} else {
				w1.add(currW / 2);
				w2.add(currW / 2);
			}
		}

		ArrayList<MultiWeightData> mwdList = new ArrayList<>(2);
		MultiWeightData mwd1 = new MultiWeightData(w1);
		MultiWeightData mwd2 = new MultiWeightData(w2);
		mwdList.add(mwd1);
		mwdList.add(mwd2);
		return mwdList;

	}

	@Override
	public String toString() {
		return w.toString();
	}

	public MultiWeightData add(MultiWeightData mwd) {
		if (mwd.numberOfWeights() != this.numberOfWeights())
			System.err.println("Error, MultiWeightData need to have same amount of weights to add them");

		for (int i = 0; i < this.numberOfWeights(); i++) {
			w.set(i, w.get(i) + mwd.w.get(i));
		}
		return this;
	}

	public MultiWeightData subtract(MultiWeightData mwd) {
		if (mwd.numberOfWeights() != this.numberOfWeights())
			System.err.println("Error, MultiWeightData need to have same amount of weights to subtract them");

		for (int i = 0; i < this.numberOfWeights(); i++) {
			w.set(i, w.get(i) - mwd.w.get(i));
		}
		return this;
	}

	public MultiWeightData makeNegative() {
		for (int i = 0; i < this.numberOfWeights(); i++) {
			w.set(i, -Math.abs(w.get(i)));
		}
		return this;
	}

	public MultiWeightData makeZero() {
		for (int i = 0; i < this.numberOfWeights(); i++) {
			w.set(i, 0l);
		}
		return this;
	}
}
