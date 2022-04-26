# Neighborhood Table

## Using `std::vector`

implemented using `std::vector`, code in neighbor_vector.cpp

- Add `vector<int> *neighbor` as the array of vectors to hold the neighborhood table.
- function `createNeighborhood` is called in every iteration of `calcVerlet` (and in initialize phase), to update the neighborhood table.  
- `calcAccel` is updated to fit the operation of neighborhood table. 
- `sign` is replaced with in-line if statement
- the method of detection of periodic boundary has been updated in line 194. 

## Using Array 

implemented using paired array, code in neighbor_array.cpp

- Add `int *First = new int[np+1]` as the array of index, whose content is the first index of the its neighbor in `Edge` array
- Add `int *Edge = new int[np*np]` as the array of neighbor, here we use the largest possible number of edges as the space to store edges.
- `calcAccel` is updated to fit the operation of neighborhood table. 
- `sign` is replaced with in-line if statement