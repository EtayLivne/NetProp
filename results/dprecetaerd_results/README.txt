l1_l2: first dict is when distances in the liquids in each node between knockout and original networks are measured with l1 distances, and seconds dict is for measure in l2 distance. 
	   the subheaders refer to how the ranking of each knockout is performed. if l1, just rank by sum of flows (equivalent to not using flows from different sources). If l2, take root of sum of squares.
	   the second method will emphasize spikes over a high average, i.e a 10 and a 0 will be prefered over two fives.

single_knockout - each test only knocked out a single gene