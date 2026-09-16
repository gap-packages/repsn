#@local calls, G, H, chi, rep, ext
gap> START_TEST("preimages.tst");

# ModuleBasis knows preimages of its basis elements, so extending a
# representation must not compute any. Since GAP 4.17 PreImagesRepresentative
# tests membership in Range and Image, which for a matrix group over a
# cyclotomic field costs more than the extension itself.
gap> calls:= 0;;
gap> InstallMethod( PreImagesRepresentative, "count calls (repsn test)",
>        FamRangeEqFamElm, [ IsGeneralMapping, IsObject ], 100000,
>        function( map, elm ) calls:= calls + 1; TryNextMethod(); end );
gap> G:= AlternatingGroup( 6 );;
gap> H:= Group( [ (1,2,3,4,6), (1,4)(5,6) ] );;
gap> chi:= Irr( G )[2];;
gap> rep:= IrreducibleAffordingRepresentation(
>              RestrictedClassFunction( chi, H ) );;
gap> calls:= 0;;
gap> ext:= ExtendedRepresentation( chi, rep );;
#I  Need to extend a representation of degree 5. This may take a while.
gap> IsAffordingRepresentation( chi, ext );
true
gap> calls;
0

#
gap> STOP_TEST("preimages.tst");
