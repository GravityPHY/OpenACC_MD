# Neighborhood Table

## Using `std::vector`

implemented using `std::vector`

- Add `vector<int> *neighbor` as the array of vectors to hold the neighborhood table.
- function `createNeighborhood` is called in every iteration of `calcVerlet` (and in initialize phase), to update the neighborhood table.  
- `calcAccel` is updated to fit the operation of neighborhood table. 
- `sign` is replaced with in-line if statement
- the method of detection of periodic boundary has been updated in line 194. 