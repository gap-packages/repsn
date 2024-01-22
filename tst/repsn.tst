#@local level, G, chi, rep
gap> START_TEST("repsn.tst");

# switch off warnings
gap> level:= InfoLevel( InfoWarning );;
gap> SetInfoLevel( InfoWarning, 0 );

# IrreducibleAffordingRepresentation:
# linear character
gap> G:= AlternatingGroup( 5 );;
gap> chi:= TrivialCharacter( G );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true

# perfect group, nonlinear character with nontrivial kernel
gap> G:= SL( 2, 5 );;
gap> chi:= First( Irr( G ), x -> x[1] = 3 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
false

# perfect group, nonlinear character with trivial kernel
gap> G:= SL( 2, 5 );;
gap> chi:= First( Irr( G ), x -> x[1] = 2 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
true

# perfect group, abelian normal subgroup N not < Z(G), nonlinear character with nontrivial kernel
gap> G:= PerfectGroup(960, 1);;
gap> chi:= First( Irr( G ), x -> x[1] = 3 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
false

# perfect group, abelian normal subgroup N not < Z(G), nonlinear character with trivial kernel
gap> G:= PerfectGroup(960, 1);;
gap> chi:= First( Irr( G ), x -> x[1] = 15 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
true

# not perfect group, nontrivial kernel, induced
gap> G:= DihedralGroup( 12 );;
gap> chi:= First( Irr( G ), x -> x[1] = 2 and Size( Kernel( x ) ) > 1 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
false

# not perfect group, nontrivial kernel, extended
gap> G:= SmallGroup( [ 48, 32 ] );; # 2xSL(2,3);;
gap> chi:= First( Irr( G ), x -> x[1] = 2 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
false

# not perfect group, trivial kernel, induced
gap> G:= DihedralGroup( 6 );;
gap> chi:= First( Irr( G ), x -> x[1] = 2 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
true

# not perfect group, trivial kernel, extended
gap> G:= SL( 2, 3 );;
gap> chi:= First( Irr( G ), x -> x[1] = 2 );;
gap> rep:= IrreducibleAffordingRepresentation( chi );;
gap> IsOne( Image( rep, One( G ) ) );
true
gap> IsTrivial( Kernel( chi ) );
true

# not perfect group, nontrivial kernel, from perfect subgroup: example?
# not perfect group, trivial kernel, from perfect subgroup: example?

# switch on warnings if applicable
gap> SetInfoLevel( InfoWarning, level );;

#
gap> STOP_TEST("repsn.tst");
