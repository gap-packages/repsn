gap> START_TEST("manual.tst");

# one quick test, based on a manual example
gap> G := SymmetricGroup( 6 );;
gap> H := AlternatingGroup( 6 );;
gap> chi := Irr( H )[ 2 ];;
gap> rep := IrreducibleAffordingRepresentation( chi );;
gap> A := Image(rep);;
gap> IsIntegerMatrixGroup(A);
true
gap> DimensionOfMatrixGroup(A);
5

#
gap> sub := InducedSubgroupRepresentation( G, rep );;
gap> B := Image(sub);;
gap> IsIntegerMatrixGroup(B);
true
gap> DimensionOfMatrixGroup(B);
10

#
gap> STOP_TEST("manual.tst", 0);
