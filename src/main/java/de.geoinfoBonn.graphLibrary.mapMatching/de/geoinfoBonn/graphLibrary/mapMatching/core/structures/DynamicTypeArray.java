package de.geoinfoBonn.graphLibrary.mapMatching.core.structures;

import java.util.Arrays;

/**
 * An implementation for dynamically sized arrays. Alternative to ArrayLists
 * which allows to use primitive types. For usage, a python script needs to be
 * used to transform template "Type" to wanted datatype.
 * 
 * @author Axel Forsch
 * @version 1.2
 * @since 02. Sep. 2020
 */
public class DynamicTypeArray {
	Type[] content;
	int idx;

	/**
	 * Initializes the dynamic array with one empty array element.
	 */
	public DynamicTypeArray() {
		this(1);
	}

	/**
	 * Initializes the dynamic array with the given initial size
	 * 
	 * @param initialSize number of elements in the initial array
	 */
	public DynamicTypeArray(int initialSize) {
		this.content = new Type[initialSize];
		this.idx = 0;
	}

	/**
	 * Copy constructor doing a deep copy of the element.
	 * 
	 * @param copy
	 */
	public DynamicTypeArray(DynamicTypeArray copy) {
		this.content = Arrays.copyOf(copy.content, copy.content.length);
		this.idx = copy.idx;
	}

	/**
	 * Initializes the dynamic array with the given data. Capacity is Initialized
	 * exactly fitting to the content.
	 * 
	 * @param content Initial data of the array
	 */
	public DynamicTypeArray(Type[] content) {
		this.content = content;
		idx = content.length;
	}

	/**
	 * Sets the content at specified index to the new element.
	 * 
	 * @param index   index of content array that should be altered
	 * @param element new content
	 * @return true if set was successful, else false
	 */
	public boolean set(int index, Type element) {
		if (index < this.idx) {
			content[index] = element;
			return true;
		} else if (index == this.idx) {
			this.add(element);
			return true;
		} else {
			System.err.println("Index to big, use add to enlarge DynamicArray or add parameter for fill element.");
			return false;
		}
	}

	/**
	 * Sets the content at specified index to the new element.
	 * 
	 * @param index   index of content array that should be altered
	 * @param element new content
	 * @return true if set was successful, else false
	 */
	public boolean set(int index, Type element, Type fill) {
		if (index < this.idx) {
			content[index] = element;
			return true;
		} else if (index == this.idx) {
			this.add(element);
			return true;
		} else {
			while (index > this.idx)
				this.add(fill);
			this.add(element);
			return true;
		}
	}

	/**
	 * Method adds an element to the end of the array, if array is full, the array's
	 * size is doubled before the element is added.
	 * 
	 * @param element new element to be added
	 * @since 17. April 2018
	 */
	public void add(Type element) {
		// if content is full, double content size
		if (idx == this.capacity()) {
			Type[] doubledContent = new Type[this.capacity() * 2];
			System.arraycopy(content, 0, doubledContent, 0, this.capacity());
			content = doubledContent;
		}

		// Add content at current index position
		content[idx++] = element;
	}

	/**
	 * Removes the specified element from the vector, halves capacity if possible
	 * 
	 * @param index element to be removed
	 */
	public void remove(int index) {
		System.arraycopy(content, index + 1, content, index, this.capacity() - 1 - index);
		idx--;

		// if size needed is smaller than half the capacity, halve it
		if (this.size() < this.capacity() / 2) {
			Type[] halvedContent = new Type[this.capacity() / 2];
			System.arraycopy(content, 0, halvedContent, 0, this.capacity() / 2);
			content = halvedContent;
		}
	}

	/**
	 * Changes capacity to the given value.
	 * 
	 * @param capacity new capacity
	 */
	public void setCapacity(int capacity) {
		if (capacity < this.size()) {
			// System.out.println("Warning! Capacity smaller than size, content is cropped!
			// (in tests.DynamicTypeArray.setCapacity(int))");
		}
		Type[] newContent = new Type[capacity];
		System.arraycopy(content, 0, newContent, 0, Math.min(capacity, this.size()));
		content = newContent;
	}

	/**
	 * Trims the array, so that capacity equals size.
	 */
	public void trim() {
		this.setCapacity(this.size());
	}

	/**
	 * Returns the specified element
	 * 
	 * @param index index of returned element
	 * @return the searched element
	 */
	public Type get(int index) {
		return content[index];
	}

	/**
	 * Returns last element of the array
	 * 
	 * @return last element
	 */
	public Type getLast() {
		return get(idx - 1);
	}

	/**
	 * Returns first element of the array
	 * 
	 * @return first element
	 */
	public Type getFirst() {
		return get(0);
	}

	/**
	 * Returns the number of filled slots in the array.
	 * 
	 * @return number of non-null elements in content
	 */
	public int size() {
		return idx;
	}

	/**
	 * Returns size of the data structure, including empty slots.
	 * 
	 * @return size of content array
	 */
	public int capacity() {
		return content.length;
	}

	/**
	 * Checks if elements are present in the array
	 * 
	 * @return true if no elements are in array
	 */
	public boolean isEmpty() {
		if (this.size() == 0)
			return true;
		return false;
	}
}
