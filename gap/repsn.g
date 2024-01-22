############################################################################################
############################################################################################
############################################################################################
##                                                                                         #
##                                         REPSN                                           #
##                                                                                         #
##                                     A GAP4 Package                                      #
##                      for Constructing Representations of Finite Groups                  #
##                                                                                         #
##                Developed by: Vahid Dabbaghian                                           #
##                              Department of Mathematics,                                 #
##                              Simon Fraser University,                                   #
##                              Burnaby, BC, Canada.                                       #
##                              E-mail: vdabbagh@sfu.ca                                    #
##                              Home Page: www.sfu.ca/~vdabbagh                            #
##                                                                                         #
############################################################################################



############################################################################################
############################################################################################
################################### A Main Function ########################################
## This program returns a representation of a group $G$ affording $\chi$.
## In this program the following programs and subroutines are called.
## 1. CharacterSubgroupRepresentation
## 2. PerfectRepresentation
## 3. ReducedSubgroupCharacter
## 4. InducedSubgroupRepresentation
## 5. ExtendedRepresentationNormal
## 6. IrreducibleAffordingRepresentation
############################################################################################

 InstallGlobalFunction( IrreducibleAffordingRepresentation,
 function( chi )
 local n, k, g, c, rep, G;
 G := UnderlyingGroup( chi );
 if chi[1] = 1 then
    return CharacterSubgroupRepresentation( chi );
 fi;
 if Size(Kernel(chi)) > 1 then
    k := NaturalHomomorphismByNormalSubgroupNC( G, Kernel( chi ) );
    g := Image(k);
    c := Filtered(Irr(Image(k)), x-> chi=RestrictedClassFunction( x, k ) );
 else
    g := G;
    c := [chi];
 fi;
 if IsPerfect( g ) then
    rep := PerfectRepresentation( g, c[1], SpecialCoversData );
 else
    n := ReducedSubgroupCharacter( c[1] );
    if n[2][1] < c[1][1] then
      rep := InducedSubgroupRepresentation( g, IrreducibleAffordingRepresentation( n[2] ) );
    fi;
    if n[2][1] = c[1][1] then
      if IsPerfect(n[1]) and PerfectRepresentation( n[1], n[2], 1 ) = true then
        return CharacterSubgroupRepresentation( chi );
      fi;
      rep := ExtendedRepresentationNormal( c[1], IrreducibleAffordingRepresentation( n[2] ) );
    fi;
 fi;
 if IsBound( k ) then
   return CompositionMapping( rep, k );
 else
   return rep;
 fi;
 end );


############################################################################################
############################################################################################
####################################### A main Function ####################################
## This is a test function. If $R$ is a representation of $G$ and $\chi$ is an irreducible
## character of $G$ then this program returns "true" if the trace of $R(x)$ is equal to
## $\chi(x)$ for all representatives $x$ of the conjugacy classes of $G$. The inputs of this
## program are a representation rep and a character chi of a group.
############################################################################################

 InstallGlobalFunction( IsAffordingRepresentation,
 function( chi, rep )
 local tbl, ccl, i;
     tbl:= UnderlyingCharacterTable( chi );
     ccl:= ConjugacyClasses( tbl );
     for i in [ 1 .. Length( ccl ) ] do
       if Trace( Representative( ccl[i] )^rep ) <> chi[i] then
         return false;
       fi;
     od;
     return true;
 end );


############################################################################################
############################################################################################
####################################### A main Function ####################################
## This program finds all $p$-subgroups of a group $G$ which are $\chi$-subgroups.
## In this program the following programs and subroutines are called,
## 1. LinearConstituentWithMultiplicityOne
## 2. PSubgroups
############################################################################################

 InstallGlobalFunction( AllCharacterPSubgroups,
 function( G, chi )
 local i, lin, P, psub, L;
 psub := [ ];
 L := PSubgroups( G );
 for P in L do
    lin := LinearConstituentWithMultiplicityOne( chi, P );
    if lin <> fail then
       Add( psub, P );
    fi;
 od;
 return psub;
 end );


############################################################################################
############################################################################################
####################################### A main Function ####################################
## This program finds standard subgroups (Sylow subgroups and their normalizers, centralizers
## and derived subgroups) of a group $G$ which are $\chi$-subgroups.
## In this program the following programs and subroutines are called,
## 1. LinearConstituentWithMultiplicityOne
## 2. StandardSubgroups
############################################################################################

 InstallGlobalFunction( AllCharacterStandardSubgroups,
 function( G, chi )
 local i, lin, P, psub, L;
 psub := [ ];
 L := StandardSubgroups( G );
 for P in L do
    lin := LinearConstituentWithMultiplicityOne( chi, P );
    if lin <> fail then
       Add( psub, P );
    fi;
 od;
 return psub;
 end );


############################################################################################
############################################################################################
####################################### A main Function ####################################
## This program finds all $\chi$-subgroups of $G$ among the lattice of subgroups.
## In this program the following programs and subroutines are called,
## 1. LinearConstituentWithMultiplicityOne
############################################################################################

 InstallGlobalFunction( AllCharacterSubgroups,
 function( G, chi )
 local H, i, lin, l, S;
 S := [ ];
 l := List( ConjugacyClassesSubgroups( G ), Representative );
 for H in l do
   lin := LinearConstituentWithMultiplicityOne( chi, H );
   if lin <> fail then
      Add( S, H );
   fi;
 od;
 return S;
 end );


############################################################################################
############################################################################################
##################################### A Main Function ######################################
## This is a test function. If $H$ is a $\chi$-subgroup then it returns "true", otherwise
## it returns "false".
## In this program the following programs and subroutines are called,
## 1. LinearConstituentWithMultiplicityOne
############################################################################################

 InstallGlobalFunction( IsCharacterSubgroup,
 function( chi, H )
 local l;
 l := LinearConstituentWithMultiplicityOne( chi, H );
 if  l = fail then
    return false;
 else
    return true;
 fi;
 end );


############################################################################################
############################################################################################
#################################### A Main Function #######################################
## This program computes an equivalent representation to the representation rep by
## transforming rep to a new basis. If $x, y,...$ are the generators of the image
## of rep, this program finds this new basis among the elements
## $x, y, xy, x^2, y^2, x^2y,...$.
############################################################################################

 InstallGlobalFunction( EquivalentRepresentation,
 function( rep )
 local l, l1, L, LL, V, f, G, GG, U, N, u, s, i, d, row, g, S, B, Mb;
 d := 0;
 L := [ ];
 LL := [ ];
 l := [ ];
 GG := [ ];
 U := GeneratorsOfGroup( Image( rep ) );
 g := PreImagesRange( rep );
 S := GeneratorsOfGroup( g );
 f := FieldOfMatrixGroup( Image(rep) );
 V := FullRowSpace( f, DimensionsMat( U[1] )[1] );
 G := [ Identity (Image( rep ) ) ];
 Append( G, List( [ 1..Length( U ) ] ,i -> U[ i ] ) );
 while d = 0 do
   for u in U do
     for s in G do
       Add( GG , s*u );
       GG := Set( GG );
     od;
   od;
   G := GG;
   Mb := MutableBasis( f, List( [ 1..Length(GG) ], i -> GG[i][1] ));
   LL := BasisVectors( Mb );
   l1 := TransposedMat( LL );
   N := LinearIndependentColumns( l1 );
   if Length(N) = Dimension( V ) then
   for i in N do
     Add( l, LL [i] );
   od;
   B := Basis( V, l );
   for u in U do
     row := List( [ 1..Dimension(V) ] , i -> Coefficients (  B, l[i]*u ) );
     Add( L, row );
   od;
     return GroupHomomorphismByImagesNC( g, Group( L ), S, L );
   fi;
 od;
 end );


############################################################################################
############################################################################################
################################### A Main Function ########################################
## This is a program for constructing representation of finite groups according
## to Dixon's method. The input of this program are a group $G$ and an irreducible
## character $\chi$ of $G$, or $G$, $\chi$ and a $\chi$-subgroup.
## In this program the following programs and subroutines are called.
## 1. CharacterSubgroup
## 2. CharacterSubgroupLinearConstituent
## 3. InducedSubgroupRepresentation
## 4. RepresentationMatrixGenerators
## 5. CharacterSubgroupRepresentation
############################################################################################

InstallGlobalFunction( CharacterSubgroupRepresentation,
 function( arg )
    local  G, U, r, chi, s, S, H, f, M, m1, Y, rep;
    chi := arg[1];
    G := UnderlyingGroup( chi );
    U := GeneratorsOfGroup( G );
    if chi[1] = 1 then
        if Length( U ) = 0 then
           if IsPermGroup(G) then
              return GroupHomomorphismByImagesNC( Group(()), Group([[1]]), [()], [[[1]]] );
           else
              return GroupHomomorphismByImagesNC( G, G, U, U );
           fi;
        fi;
        r := List( [ 1..Length( U ) ],i -> [ [ U[ i ] ^ arg[ 1 ] ] ] );
        return GroupHomomorphismByImagesNC( G, Group( r ), U, r );
    fi;
    if Length( arg ) = 1  then
       if IsClassFunction( chi ) and IsCharacter( chi ) then
          s := CharacterSubgroup( G, chi );
          H := s[1];
          f := s[2];
       else
          Error("first argument must be an ordinary character");
       fi;
    fi;
    if Length( arg ) = 2  then
        if IsSubgroup ( G, arg[2] ) and IsClassFunction( chi ) and IsCharacter( chi ) then
        H := arg[2];
        s := CharacterSubgroupLinearConstituent( chi, [ H ] );
        if s = fail then
           Error( "the given group is not a chi-subgroup");
        fi;
        H := s[1];
        f := s[2];
        else
           Error("first argument is not a character or second argument is not a subgroup of the given group");
    fi; fi;
    if Index( G, H ) = chi[1] then
    rep := CharacterSubgroupRepresentation( f );
    return InducedSubgroupRepresentation ( G, rep );
    fi;
    M := MatrixRepsn( chi, G, H, f);
    m1 := M[1];
    Y := M[2];
    r := RepresentationMatrixGenerators( U, Y, m1, f, chi, H );
    return GroupHomomorphismByImagesNC( G, Group( r ), U, r );
 end );


############################################################################################
############################################################################################
#################################### A Main Function #######################################
## If rep is a representation of a subgroup $H$ of a group $G$, this programs
## returns an induced representation of rep of group $G$.
############################################################################################

InstallGlobalFunction( InducedSubgroupRepresentation,
 function( G, rep )
 local M, M1, mat, u, U, rt, i, j, k, LL, L, m, Ls, Lr, l, d, h, r, t, nmat, lrt;
 mat := [ ];
 L := [ ];
 h := PreImagesRange( rep );
 if Size(h) = 1 then
   d := 1;
 else
   d := DimensionsMat( GeneratorsOfGroup( Image(rep) )[1] )[1];
 fi;
 U := GeneratorsOfGroup( G );
 rt := RightTransversal( G, h );
 lrt := Length( rt );
 nmat := NullMat(d,d);
 for u in U do
  M := [ ];
  for i in [1..lrt] do
   t := rt[i]*u;
   M1 := [ ];
    for j in [1..lrt] do
      r := t*rt[j]^(-1);
      if r in h then
        Add( M1, r^rep );
      else
        Add( M1, nmat );
      fi;
    od;
   Add( M, M1 );
  od;
  Add( mat, M );
 od;
 for m in mat do
 LL := [ ];
 Lr := [ ];
 Ls := [ ];
  for i in [1..lrt] do
    for j in [1..d] do
      for k in [1..lrt] do
          Add( LL, m[i][k][j] );
      od;
    Add( Ls, Concatenation(LL) );
    LL := [ ];
    od;
    Add( Lr, Ls );
  od;
  Add( L, Lr );
 od;
 l := List( L, i-> i[1] );
 return GroupHomomorphismByImagesNC( G, Group( l ), U, l );
 end );


############################################################################################
############################################################################################
#################################### A Main Function #######################################
## The inputs of this program are an irreducible character $\chi$ of a group $G$ and
## a representation rep of a subgroup $H$ of $G$ affording the irreducible character
## $\chi_H$.
## In this program the following programs and subroutines are called,
## 1. ModuleBasis
############################################################################################

 InstallGlobalFunction( ExtendedRepresentation,
 function( chi, rep )
 local S, C, mat, B, i, j, s, u, U, L, n, M, l, h, G, Mat, dim;
 S := [ ];
 C := [ ];
 mat := [ ];
 L := [ ];
 n := chi[1];
 h := PreImagesRange( rep );
 if Size(h) = 1 then
    B := [ [ [ [ 1 ] ] ], [ Identity( h ) ] ];
 else
    Mat := GeneratorsOfGroup(Image(rep));
    dim := DimensionsMat( Mat[1] )[1];
    Info( InfoWarning, 1, "Need to extend a representation of degree ", dim, ". This may take a while." );
    B := ModuleBasis( h, rep );
 fi;
 G := UnderlyingGroup( chi );
 U := GeneratorsOfGroup( G );
 for u in U do
   for i in [1..n^2] do
     for j in [1..n] do
       Append( L, TransposedMat( B[1][i] )[j] );
     od;
     Add( mat, L );
     L := [ ];
     Add( C, (B[2][i]*u)^chi );
   od;
   s := SolutionMat( TransposedMat(mat), C );
   M:= List( [ 1 .. n ], i -> s{ [ (i-1)*n+1 .. i*n ] } );
   Add( S, M );
   mat := [ ];
   C := [ ];
 od;
 return  GroupHomomorphismByImagesNC( G, Group( S ), U, S );
 end );


############################################################################################
############################################################################################
#################################### A Main Function #######################################
## The inputs of this program are an irreducible character $\chi$ of a group $G$ and
## a representation rep of a normal subgroup $H$ of $G$ affording the irreducible
## character $\chi_H$.
############################################################################################

 InstallGlobalFunction( ExtendedRepresentationNormal,
 function( chi, rep )
 local R, S, C, G, x, M1, M2, M, Mat, i, s, u, U, Uh, h, b, w, o, n, dim;
 C := [ ];
 n := chi[1];
 h := PreImage( rep );
 Uh:= GeneratorsOfGroup( h );
 U := ShallowCopy( Uh );
 Mat := List( [1..Length(U)], i -> U[i]^rep );
 Mat := ShallowCopy( Mat );
 if Size( h ) = 1 then
   dim := 1;
 else
   dim := DimensionsMat( Mat[1] )[1];
   Info( InfoWarning, 1, "Need to extend a representation of degree ", dim, ". This may take a while." );
 fi;
 G := UnderlyingGroup( chi );
 repeat
   repeat
     x:= PseudoRandom( G );
   until not x in h and x^chi <> 0;
   for u in Uh do
     R := u^rep;
     S := ( x^(-1)*u*x )^rep;
     M1 := KroneckerProduct( R, IdentityMat(n) );
     M2 := KroneckerProduct( IdentityMat(n), TransposedMat(S) );
     Append( C, M1 - M2);
   od;
     b := IdentityMat(n);
     b := Concatenation(b);
     Add( C, b );
     w := List( [1..Length(Uh)*n^2] , i->0 );
     Add( w, x^chi );
     s := SolutionMat( TransposedMat(C), w );
     M:= List( [ 1 .. n ], i -> s{ [ (i-1)*n+1 .. i*n ] } );
     Add( Mat, M );
     C := [ ];
     Add( U, x );
 until Size ( Group(U) ) = Size( G );
 return GroupHomomorphismByImagesNC( Group(U), Group( Mat ), U, Mat );
 end );


############################################################################################
############################################################################################
################################### A Main Function ########################################
## This program returns constituents and multiplicities a representation $rep$.
## In this program the following program is called.
## 1. IrreducibleAffordingRepresentation
############################################################################################

   InstallGlobalFunction( ConstituentsOfRepresentation,
 function( rep )
   local D, G, ctb, chi, con, R, C, m, i, r, UG, l;
   D := [ ];
   G := PreImagesRange( rep );
   UG := GeneratorsOfGroup( G );
   ctb := CharacterTable( G );
   con := ConjugacyClasses( ctb );
   chi := ClassFunction( ctb, List( con , i -> TraceMat( Image( rep, Representative( i )))));
   R := RestrictedClassFunction( chi, G );
   C := ConstituentsOfCharacter( R );
   m := MatScalarProducts( C, [ R ] );
   for i in [1..Length( C )] do
      r := IrreducibleAffordingRepresentation( C[i] );
      l := List( UG, i -> i^r );
      Add( D, [ m[1][i], GroupHomomorphismByImagesNC( G, Group( l ), UG, l ) ] );
   od;
   return D;
 end );


############################################################################################
############################################################################################
################################### A Main Function ########################################
## This program computes an equivalent block diagonal representation to a given
## representation, or computes a block diagonal representation of the given list of
## irreducible representations.
## In this program the following program is called.
## 1. DiagonalBlockMatrix
############################################################################################

   InstallGlobalFunction( EquivalentBlockRepresentation,
 function( reps )
   local G, UG, U, i, j, k, matlist, mat, Mat;
   Mat := [ ]; U := [ ];
   if not IsList( reps ) then
      reps := ConstituentsOfRepresentation( reps ) ;
   fi;
   if Length( reps ) = 1 and reps[1][1] = 1 then
      return reps[1][2];
   fi;
   G := PreImagesRange( reps[1][2] );
   UG := GeneratorsOfGroup( G );
   for k in [1..Length( reps )] do
       for j in [1..reps[k][1]] do
            Add( U, GeneratorsOfGroup( Image( reps[k][2] ) ) );
       od;
   od;
   for i in [1..Length( U[1] )] do
      matlist := List( U, j-> j[i] );
      mat := DiagonalBlockMatrix( matlist );
      Add( Mat, mat );
   od;
   return GroupHomomorphismByImagesNC( G, Group( Mat ), UG, Mat );
 end );


############################################################################################
############################################################################################
################################### A Main Function ########################################
## This program returns true if the representation $rep$ is reducible.
############################################################################################

 InstallGlobalFunction( IsReducibleRepresentation,
 function( rep )
   local G, ctb, chi, con;
   G := PreImagesRange( rep );
   ctb := CharacterTable( G );
   con := ConjugacyClasses( ctb );
   chi := ClassFunction( ctb, List( con , i -> TraceMat( Image( rep, Representative( i )))));
   if IsIrreducibleCharacter( chi ) = false then
      return true;
   else
      return true;
   fi;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
##
##

 InstallGlobalFunction( DiagonalBlockMatrix,
 function( matrices )
   local dimension,  partialSum,  totalLength,  result,  m,  i,  range;
   dimension := List( matrices, Length );
   partialSum := List([0..Length( matrices )], x -> Sum( dimension{[1..x]} ));
   range := List( [1..Length( matrices )], x -> [ partialSum[x]+1..partialSum[x+1] ] );
   totalLength := Sum( dimension );
   result := NullMat( totalLength, totalLength, Field(matrices[1][1][1]) );
   for m in [1..Length( matrices )] do
     for i in [1..dimension[m]] do
       result[ range[m][i] ]{range[m]} := matrices[m][i];
     od;
   od;
   return result;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
##
##

 InstallGlobalFunction( AlphaRepsn,
    function( x, f, chi, H )
    local  z, sum, C, l, a, i, t;
    sum := 0; t:=0;
    C := ConjugacyClasses(H);
    l := List(C, Representative);
    for i in [1..Length(l)] do
        a := ( l[i] ^ -1 ) ^ f;
        for z in C[i] do
          t := t +  (z * x ) ^ chi;
        od;
        t := t * a;
        sum := sum + t ;
        t := 0;
    od;
    return sum;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program finds a linear constituent with multiplicity one of $\chi_H$.
## If such a constituent does not exist, it returns "false".
##

 InstallGlobalFunction( LinearConstituentWithMultiplicityOne,
 function( chi, H )
 local rest;
 rest:= RestrictedClassFunction( chi, H );
 return First( LinearCharacters( H ), lambda -> ScalarProduct( rest, lambda ) = 1 );
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program finds a $\chi$-subgroup from the list S. If such a subgroup
## does not exist it returns "fail".
##

 InstallGlobalFunction( CharacterSubgroupLinearConstituent,
 function( chi, S )
 local  j, lin;
 for j  in [1..Length(S)]  do
    lin := LinearConstituentWithMultiplicityOne( chi, S[j]);
    if lin <> fail then
       return [ S[j] ,lin ];
    fi;
 od;
 return fail;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This Program constructs the matrix A described in Dixon's paper.
##

 InstallGlobalFunction( MatrixRepsn,
    function( chi, G, H, f )
    local  mat, matn, X, Xinv, j, n, x, d, h, Y;
    h := Size(H);
    mat := [ [ h ] ];
    X := [ Identity( G ) ];
    Xinv:= [ Identity( G ) ];
    d := DegreeOfCharacter( chi );
    while Length( X ) < d  do
        matn := [  ];
        n := Length( X );
        x := PseudoRandom( G );
        for j  in [ 1 .. n ]  do
            Y := x * Xinv[j];
            matn[j] := AlphaRepsn( Y, f, chi, H );
            mat[j][n + 1] := ComplexConjugate( matn[j] );
        od;
        matn[n + 1] := h;
        mat[n + 1] := matn;
        if RankMat( mat ) = Length( mat) then
            X[n + 1] := x;
            Xinv[n + 1] := x^-1;
        fi;
    od;
    return [ mat, X ];
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program finds a $\chi$-subgroup among the Sylow subgroups of $G$.
##

 InstallGlobalFunction( CharacterSylowSubgroup,
 function( G, chi )
 local   S, x, pair, coll, sizes, lin;
 coll:= Collected( Factors( Size( G ) ) );
 sizes:= List( coll, x -> x[1]^x[2] );
 SortParallel( sizes, coll );
 for pair in coll do
 S:= SylowSubgroup( G, pair[1] );
    lin := LinearConstituentWithMultiplicityOne( chi, S );
    if lin <> fail then
       return S;
    fi;
 od;
 return fail;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program constructs the matrices of the representation, corresponding
## to set of generators U in the main Lemma of Dixon's paper.
##

 InstallGlobalFunction( RepresentationMatrixGenerators,
    function( U, Y, m1, f, chi, H )
    local  u, i, j, m2, L, X, inv, Yinv, lenY;
    L := [ ];
    lenY := Length( Y );
    Yinv := List( Y, i -> i^-1 );
    inv := Inverse( m1 );
    for u  in U  do
      m2 := [ ];
      for i in [ 1 .. lenY ] do
        m2[i] := [ ];
        for j in [ 1 .. lenY ] do
            X := ( Y[i] * u ) * Yinv[j];
            m2[i][j] := AlphaRepsn( X, f, chi, H );
        od;
     od;
     Add( L, m2 * inv );
    od;
    return L;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program finds the smallest possible $\chi$-subgroup of $G$. At the first
## it checks among the standard subgroups of $G$ and then the lattice of subgroups of the
## Sylow subgroups. If there is no subgroup as a $\chi$-subgroup then it searches among
## the lattice of subgroups of $G$.
##

 InstallGlobalFunction( CharacterSubgroup,
 function( G, chi )
 local N, s, S, pi, i, L, reps, sp, l;
      reps := [  ];
      N := StandardSubgroups( G );
      s := CharacterSubgroupLinearConstituent( chi, N );
  #   if s <> fail the return s; fi;
      if s = fail then
         pi := Set(FactorsInt( Size(G) ) );
         for i in [1..Length(pi)] do
             S := SylowSubgroup( G, pi[i] );
             sp := List(ConjugacyClassesSubgroups(S),Representative);
             L := sp[1];
             repeat
                 Add( reps, L );
                 L := First( sp, g -> not ForAny( reps, rep -> IsConjugate( G, g, rep ) ) );
                 until L = fail;
             reps := Set( reps );
             l := List( reps, Size );
             SortParallel( l , reps );
             s := CharacterSubgroupLinearConstituent( chi, reps );
        od;
      fi;
      if s = fail  then
             if not ( HasConjugacyClassesSubgroups( G ) or HasTableOfMarks( G ) ) then
                  Info( InfoWarning, 1, "need to compute subgroup lattice of a group of order ", Size( G ) );
             fi;
            S := List( ConjugacyClassesSubgroups( G ), Representative );
            s := CharacterSubgroupLinearConstituent( chi, S );
            if s = fail  then
               Error ( "Dixon's method does not work for the given character ");
            fi;
      fi;
 return s;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## The inputs of this program are a group $G$ and a representation rep of a subgroup
## $H$ of $G$ which is the restriction of a representation of $G$ to $H$.
## Using the Burnside's Theorem, this program finds a basis for the matrix algebra
## of the image of rep.
##

 InstallGlobalFunction( ModuleBasis,
 function( G, rep )
 local V, f, B, u, x, L, l, d, i, N, P;
 l := [ ];
 L := [ ];
 B := [ ];
 d := DimensionsMat( GeneratorsOfGroup( Image(rep) )[1] )[1];
 for i in [1..2*d^2] do
   x := PseudoRandom( G );
   Add( l, x^rep );
 od;
 f := FieldOfMatrixGroup( Image(rep) );
 V := MutableBasis( f, l );
 while NrBasisVectors(V) < d^2 do
   x := PseudoRandom(G);
   u := x^rep;
   if IsContainedInSpan( V, u ) = false then
              CloseMutableBasis( V, u );
              Add( l, u );
   fi;
 od;
 for i in [1..Length(l)] do
   Add( L, Concatenation( l[i] ) );
 od;
 N := LinearIndependentColumns( TransposedMat( L ) );
 for i in N do
   Add( B, l[i] );
 od;
 P := List( [1..d^2], i->PreImagesRepresentative( rep, B[i] ) );
 return [B,P];
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This is a subprogram for the function ``MinimalNormalSubgps".
##

 InstallGlobalFunction( AddToSubgpList,
 function( nt, U )
 local i, N;
 for N in nt do
   if IsSubgroup( U, N ) then
      return nt;
   fi;
 od;
 for i in [1..Length(nt)] do
   if IsSubgroup(nt[i], U) then
      nt[i] := false;
   fi;
 od;
 nt := Filtered( nt, x -> not IsBool(x) );
 Add( nt, U );
 return nt;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This function finds the minimal normal subgroups of a group $G$.
##

 InstallGlobalFunction( MinimalNormalSubgps,
 function( G )
 local cl, nt, c, U;
 cl := ConjugacyClasses( G );
 nt := [];
 for c in cl do
   if IsPrime( Order( Representative(c) ) ) then
      U  := SubgroupNC( G, [Representative(c)] );
      U  := NormalClosure( G, U );
      nt := AddToSubgpList( nt, U );
   fi;
 od;
 return nt;
 end );


############################################################################################
############################################################################################
###################################### A Utility Function ##################################
## If $\chi$ is an irreducible character of $G$, this program reduces the problem of
## constructing a representation of $G$ affording $\chi$ to a perfect subgroup of $G$.
## It finds either a subgroup $H$ and an irreducible character $\phi$ of $H$ such that
## $\chi=\phi^G$ or a subgroup $H$ such that $\chi_H$ is irreducible. The input of this
## program is $\chi$ and the outputs are either $H$ and $\phi$ or $H$ and $\chi_H$.
##

 InstallGlobalFunction( ReducedSubgroupCharacter,
 function( chi )
 local k, i, N, I, R, Ri, Rc, C, Ci, Cc, c, L, n, G;
 G := UnderlyingGroup( chi );
 N := DerivedSeries( G );
 k := 0;
 for i in [2..Length( N )] do
   R := RestrictedClassFunction( chi, N[i] );
   C := ConstituentsOfCharacter(R);
   if Length(C) > 1 then
    I := InertiaSubgroup( G, C[1] );
    Ri := RestrictedClassFunction( chi, I );
    Ci := ConstituentsOfCharacter( Ri );
    for c in Ci do
      Rc := RestrictedClassFunction( c, N[i] );
      Cc := ConstituentsOfCharacter( Rc );
      if C[1] in Cc then
        return [I,c];
      fi;
    od;
   elif IsIrreducibleCharacter( R ) then
      k := i; c := R;
   fi;
 od;
 if k > 0 then
   return [ N[k], c ];
 fi;
 L := ChiefSeriesThrough( G, [N[2]] );
 n := Length(L);
 for i in [2..n] do
   R := RestrictedClassFunction( chi, L[i] );
   C := ConstituentsOfCharacter( R );
   if Length(C) > 1 then
       I := InertiaSubgroup( G, C[1] );
       Ri := RestrictedClassFunction( chi, I );
       Ci := ConstituentsOfCharacter( Ri );
         for c in Ci do
          Rc := RestrictedClassFunction( c, L[i] );
          Cc := ConstituentsOfCharacter( Rc );
          if C[1] in Cc then
            return [I,c];
          fi;
         od;
   fi;
 od;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program lists all $p$-subgroups of $G$ which are not $G$-conjugate.
##

 InstallGlobalFunction( PSubgroups,
 function( G )
 local P, S, reps, L, rep, g, sp, l;
 reps := [ ];
 S := List( Set(FactorsInt( Size(G) ) ), i->SylowSubgroup( G, i ) );
 for P in S do
   sp := List(ConjugacyClassesSubgroups(P),Representative);
   L := sp[1];
   repeat
   Add( reps, L );
   L := First( sp, g -> not ForAny( reps, rep -> IsConjugate( G, g, rep ) ) );
   until L = fail;
 od;
 reps := Set( reps );
 l := List( reps, Size );
 SortParallel( l , reps );
 return reps;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## For all Sylow $p$-subgroups $P$ of $G$ of order $p^a$, this program lists $P$, the
## normalizer of $P$ in $G$ if $a <= 2$, and the derived subgroup $P'$ of $P$ and the
## centralizer of $P'$ in $G$ if $a > 2$.
##

 InstallGlobalFunction( StandardSubgroups,
 function( G )
 local reps, S, P, D, l, i, col;
 col := Collected( FactorsInt( Size ( G ) ) );
 S := List(col, i-> SylowSubgroup( G, i[1] ) );
 reps := S;
 for i in [1..Length(col)] do
   if col[i][2] <= 2 then
     Add( reps, Normalizer( G, S[i] ) );
   elif not IsAbelian( S[i] ) then
     D := DerivedSubgroup( S[i] );
     Add( reps, D );
     Add( reps, Centralizer( G, D ) );
   fi;
 od;
 reps := Set( reps );
 l := List( reps, Size );
 SortParallel( l , reps );
 return reps;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program is for constructing representations of groups $6.A_6$, $2.A_7$,
## $3.A_7$, $6.A_7$ and $2.A_8$, affording some special characters. These are the
## characters, say $\chi$, such that there exists a $\chi$-subgroup which is
## not a $p$-subgroup.
##

 InstallGlobalFunction( CoversOfAlternatingGroupsRepresentation,
 function( G, chi, size, hom )
 # See section 4.4 of thesis.
 local L, g, c, C, H;
 L := [ ];
 g := Image( hom );
 C := List( ConjugacyClassesSubgroups( g ), Representative );
 for c in C do
     if Size(c) = size then
       Add( L, PreImage( hom, c ) );
     fi;
 od;
 H := CharacterSubgroupLinearConstituent( chi, L )[1];
 return CharacterSubgroupRepresentation( chi, H );
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program is for constructing representations of covering groups, say $G$, of
## group $U_4(3)$, affording an irreducible character $\chi$ of degree less than 32.
## This program finds a maximal subgroup $M$ of $G$ such that $\chi_M$ is irreducible.
##

 InstallGlobalFunction( CoversOfU43Representation,
 function( chi, hom )
 local g, k, d, c, n, r, max, R;
 g := Image( hom );
    if chi[1] in [6,15] then
      k := SylowSubgroup( g, 3 );
      d := DerivedSubgroup( k );
      c := Centralizer( k, d );
      n := Normalizer( g, c );
    fi;
    if chi[1] = 21 then
      r := IsomorphicSubgroups( g, PSL(3,4) );
      n := Image( r[1] );
    fi;
 max := PreImage( hom, n );
 R := RestrictedClassFunction( chi, max );
 return ExtendedRepresentation( chi, IrreducibleAffordingRepresentation( R ) );
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program is for constructing representations of covering group $G:=3.O_7(3)$
## of group $O_7(3)$ affording an irreducible character $\chi$ of degree less than 32.
## This program finds a maximal subgroup $M$ of $G$ such that $\chi_M$ is irreducible.
##

 InstallGlobalFunction( CoversOfO73Representation,
 function( chi, hom )
 local g, P, n, C, f, l, F, t, M, max, R;
 g := Image( hom );
 P := SylowSubgroup( g, 3 );
 n := Normalizer( g, P );
 C := List( ConjugacyClasses( n ), Representative );
 f := Filtered( C, i-> i in P );
 l := List( f, i-> NormalClosure( n, Group( i ) ) );
 F := Filtered( l, i-> Size( i ) = 243 );
 for t in F do
   M := Normalizer( g, t );
   if Index( g, M ) = 364 then
     max := PreImage( hom, M );
     R := RestrictedClassFunction( chi, max );
     return ExtendedRepresentation( chi, IrreducibleAffordingRepresentation( R ) );
   fi;
 od;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program is for constructing representations of covering group $G:=3.U_6(2)$
## of group $U_6(2)$ affording an irreducible character $\chi$ of degree less than 32.
## This program finds a maximal subgroup $M$ of $G$ such that $\chi_M$ is irreducible.
##

 InstallGlobalFunction( CoversOfU62Representation,
 function( chi, hom )
 local g, k, u, c, n, m, R;
 g := Image( hom );
 k := SylowSubgroup( g, 2 );
 u := UpperCentralSeries( k );
 c := Centralizer( k, u[3] );
 n := Normalizer( g, c );
 m := PreImage( hom, n );
 R := RestrictedClassFunction( chi, m );
 return ExtendedRepresentation( chi, IrreducibleAffordingRepresentation( R ) );
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## For computing projective representations based on Minkwitz's method
##

 InstallGlobalFunction( ProjectiveReps,
 function( chi, rep, g )
 local d, N, mat, y, G, e, T, ZG;
 G := UnderlyingGroup( chi );
 ZG := Center(G);
 N := PreImages( rep );
 d := DimensionsMat( GeneratorsOfGroup( Image(rep) )[1] )[1];
 e := chi[1]/d;
 T := List(RightTransversal(N,ZG),i->CanonicalRightCosetElement(ZG,i));
 mat := (d/(e*Size(T))) * Sum(List(T, y-> (g*(y^(-1)))^chi * ImagesRepresentative( rep,y ) ) );
 return mat;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## Group generators for Minkwitz's method
##

 InstallGlobalFunction( GeneratorsProjective,
 function( chi, rep )
 local z, N, U, g, e, ZG, G;
 G := UnderlyingGroup( chi );
 ZG := Center(G);
 N := FittingSubgroup( G );
 U := [];
 z := Size(G);
 repeat g:=Random( Difference( G, N ) );
    if g^chi<>0 then
       Add(U,g);
    fi;
 until U<>[] and Size(Group(U))=z;
 return U;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## A program for finding suitable ordinary representation of the cover of G/F to use in
## Schur's theorem
##

 InstallGlobalFunction( FilteringCharacters,
 function( chi, rep, U, Phi, h, EpiGf, IsoGf, SGf )
 local LU, d, e, l, le, ls, conj, G, Uproj, Iso, phi, Iphi, lUSGf, mu, A, S3, WS, s, s1, s2, s3, g, cfp, a, WsGf, j, TrA, tst;
 LU := Length(U);
 d := DimensionsMat( GeneratorsOfGroup( Image(rep) )[1] )[1];
 e := chi[1]/d;
 G := Group(U);                                                                  ###### G=<x1,...,xn> and U={x1...,xn} ######
 Uproj := List( U, i-> ProjectiveReps( chi, rep, i) );                           ###### [A1,...An] where Ai=V(xi) for a Proj. Reps V of G ######
 l  := List( [1..LU], i->ImagesRepresentative( h, U[i] ) );
 le := List( [1..LU], i->PreImagesRepresentative( EpiGf, l[i] ) );
 ls := List( [1..LU], i->ImagesRepresentative( IsoGf, le[i] ) );                  ###### \bar{x1},...,\bar{xn} in the cover of G/F ######
 conj := List( ConjugacyClasses(SGf), Representative );
 Iso := IsomorphismFpGroupByGenerators(G,U);
 for phi in Phi do
   Iphi := IrreducibleAffordingRepresentation( phi );
   lUSGf := List( [1..LU], i->ImagesRepresentative( Iphi, ls[i] ) );
   mu := List(lUSGf, t -> Trace(t));
   if ForAll(mu, t -> (t <> 0)) then
     A:=[]; S3:=[]; WS:=[];
     for s in conj do
        s1:= PreImagesRepresentative( IsoGf, s );
        s2:= ImagesRepresentative( EpiGf, s1 );
        s3:= PreImagesRepresentative( h, s2 );                                   ####### \bar{s} in the cover ######
        g := ImagesRepresentative( Iso, s3 );                                    ####### \bar{s}=w(x1,...xn)  ######
        cfp := ExtRepOfObj(g);
        Add(S3,s3);
        a := IdentityMat( d );
        for j in [1..Length(cfp)/2] do
            a := a * (e/(mu[cfp[2*j-1]])*Uproj[cfp[2*j-1]])^cfp[2*j];            ####### A=w(A1,...,An) #######
        od;
        Add(A,a);
        WsGf := Identity( SGf );
        for j in [1..Length(cfp)/2] do
            WsGf := WsGf * (ls[cfp[2*j-1]])^cfp[2*j];                            ####### WsGf=w(\bar{x1},...,\bar{xn})  #######
        od;
        Add(WS, WsGf);
     od;
     TrA := List( [1..Length(conj)], i->Trace(A[i]) );
     tst:= Filtered( [1..Length(conj)], i->TrA[i]<>0 );
     if List( tst, i-> WS[i]^phi) = List( tst, i-> S3[i]^chi/TrA[i] ) then
        return [lUSGf, mu];
     fi;
   fi;
 od;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## A program for computing representation of a perfect group $G$ with abelian
## $Soc(G/Z(G))$ affording \chi, in the case that \chi_F=e*\theta and e not prime.
##

 InstallGlobalFunction( IrrRepsNotPrime,
 function(chi,theta,f)
 local G, h, e, Gf, EpiGf, SchurGf, IsoGf, SGf, chiSGf, Vf, U, lVG, SGfReps, r, i;
    G := UnderlyingGroup( chi );
    h := NaturalHomomorphismByNormalSubgroup( G, f );
    e:= chi[1]/theta[1];
    Gf := Image(h);
    EpiGf := EpimorphismSchurCover( Gf, [2,3,5,7,11] );
    SchurGf := PreImage(EpiGf);
    IsoGf := IsomorphismPermGroup(SchurGf);
    SGf := Image(IsoGf);                                                   ####### SGf is a perm. reps of a cover of G/F #######
    chiSGf:= Filtered( Irr(SGf), i->i[1] = e);
    Vf := IrreducibleAffordingRepresentation( theta );
## The following 2 lines is for computing a projective representation of G of degree \theta(1)
## Here it computes the matrices V(x) for generators x in U
   U := GeneratorsProjective( chi, Vf );
   lVG := List( U, x->ProjectiveReps( chi, Vf, x ) );
## The program "FilteringCharacters" checks to find the suitable representation of
## the central cover of G/F
   SGfReps := FilteringCharacters( chi, Vf, U, chiSGf, h, EpiGf, IsoGf, SGf );
   r := List( [1..Length(U)], i->(e/SGfReps[2][i])*KroneckerProduct( SGfReps[1][i], lVG[i] ) );
   return GroupHomomorphismByImagesNC( Group(U), Group( r ), U, r );
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program returns a representation of a perfect group $G$ with abelian $Soc(G/Z(G))$
## using Gallagher's Theorem
##

 InstallGlobalFunction( GallagherRepsn,
 function( G, chi, h, f, C, m, U )
 local IrrG, i, R, psi, e1, r1, im, Irrim, e2, r2, r, l;
 IrrG := Irr(G);
 psi:=0;
 for i in [1..Length(IrrG)] do
    R := RestrictedClassFunction( IrrG[i], f );
    if R = C[1] then
       psi := IrrG[i];
       e1 := IrreducibleAffordingRepresentation( psi );
       r1 := List( U, x->ImagesRepresentative( e1, x ) );
       break;
    fi;
 od;
 if psi <> 0 then
   im := Image( h );
   Irrim := Irr(im);
   for i in [1..Length( Irrim )] do
      if m[1][1] = Irrim[i][1] then
        R := RestrictedClassFunction( Irrim[i], h );
        if chi = psi*R then
           e2 := IrreducibleAffordingRepresentation( Irrim[i] );
           l := List( U, x->ImagesRepresentative( h, x ) );
           r2 := List( l, x->ImagesRepresentative( e2, x ) );
           r := List( [1..Length(r1)], i->KroneckerProduct( r1[i], r2[i] ) );
           return GroupHomomorphismByImagesNC( G, Group( r ), U, r );
        fi;
      fi;
   od;
 fi;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program returns a representation of a perfect group $G$ with abelian
## $Soc(G/Z(G))$. This function calls the following three functions:
## 1. ExtendedRepresentationNormal
## 2. InducedSubgroupRepresentation
## 3. GallagherRepsn
##

 InstallGlobalFunction( PerfectAbelianSocleRepresentation,
 function( G, chi )
 local Ug, U, f, R, C, m, I, t, Rc, Cc, h, P;
 Ug := GeneratorsOfGroup( G );
 U := ShallowCopy( Ug );
 f := FittingSubgroup( G );
 R := RestrictedClassFunction( chi , f );
 C := ConstituentsOfCharacter( R );
 m := MatScalarProducts( C, [ R ] );
 ######### if \chi_F is irreducible #########
 if IsIrreducibleCharacter( R ) then
    return ExtendedRepresentationNormal( chi, IrreducibleAffordingRepresentation( R ) );
 fi;
 ######### if \chi is imprimitive #########
 if Length( C ) > 1 then
   I := InertiaSubgroup( G, C[1] );
   R := RestrictedClassFunction( chi, I );
   C := ConstituentsOfCharacter( R );
   for t in C do
     Rc := RestrictedClassFunction( t, f );
     Cc := ConstituentsOfCharacter( Rc );
     if C[1] in Cc then
        return InducedSubgroupRepresentation( G, IrreducibleAffordingRepresentation( t ));
     fi;
   od;
 fi;
 h := NaturalHomomorphismByNormalSubgroup( G, f );
 ######### if \chi_F=e*\theta for e prime #########
 if Length( C ) = 1  and IsPrime( m[1][1] ) then
    P := PreImage( h , SylowSubgroup( Image(h), m[1][1] ) );
    R := RestrictedClassFunction( chi , P );
    if IsIrreducibleCharacter( R ) then
       return ExtendedRepresentation( chi, IrreducibleAffordingRepresentation( R ) );
    else
       return GallagherRepsn( G, chi, h, f, C, m, U );
    fi;
 fi;
 ######### if \chi_F=e*\theta for e not prime #########
 if Length( C ) = 1 and IsPrime( m[1][1] ) = false then
    return IrrRepsNotPrime( chi, C[1], f );
 fi;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program returns a representation of a perfect group $G$ with abelian
## $Soc(G/Z(G))$. This function calls the following three functions:
## 1. ExtendedRepresentationNormal3
## 2. InducedSubgroupRepresentation3
## 3. GallagherRepsn
##

 InstallGlobalFunction( SocleMoreComponents,
 function( G, chi, S )
 local Ug, U, R, C, m, I, t, Rc, Cc, h, P;
 Ug := GeneratorsOfGroup( G );
 U := ShallowCopy( Ug );
 R := RestrictedClassFunction( chi , S );
 C := ConstituentsOfCharacter( R );
 m := MatScalarProducts( C, [ R ] );
 ######### if \chi_S is irreducible #########
 if IsIrreducibleCharacter( R ) then
    return ExtendedRepresentationNormal( chi, IrreducibleAffordingRepresentation( R ) );
 fi;
 ######### if \chi is imprimitive #########
 if Length( C ) > 1 then
   I := InertiaSubgroup( G, C[1] );
   R := RestrictedClassFunction( chi, I );
   C := ConstituentsOfCharacter( R );
   for t in C do
     Rc := RestrictedClassFunction( t, S );
     Cc := ConstituentsOfCharacter( Rc );
     if C[1] in Cc then
        return InducedSubgroupRepresentation( G, IrreducibleAffordingRepresentation( t ));
     fi;
   od;
 fi;
 h := NaturalHomomorphismByNormalSubgroup( G, S );
 ######### if \chi_F=e*\theta for e in [ 2, 3] #########
 if Length( C ) = 1  and m[1][1] in [ 2, 3 ] then
    P := PreImage( h , SylowSubgroup( Image(h), m[1][1] ) );
    R := RestrictedClassFunction( chi , P );
    if IsIrreducibleCharacter( R ) then
       return ExtendedRepresentation( chi, IrreducibleAffordingRepresentation( R ) );
    fi;
 else
    return GallagherRepsn( G, chi, h, S, C, m, U );
 fi;
 end );


############################################################################################
############################################################################################
##################################### A Utility Function ###################################
## This program returns a representation of a perfect group $G$.
##

 InstallGlobalFunction( PerfectRepresentation,
 function( arg )
 local  G, chi, hom, n, a, r, s, S, c, C, i, R, l, L, M, Co, cl, so, N, x, I, t, Ri, Ci, Rc, Cc, run, ls, j, k, A;
 G := arg[1];
 chi := arg[2];
 hom := NaturalHomomorphismByNormalSubgroup( G, Centre( G ) );
 L := [ ]; M := [ ];
 n := MinimalNormalSubgroups( Image(hom) );
 a := Filtered( n, IsAbelian );
 ## First program checks to find an abelian normal subgroup N not < Z(G).
 if a=n then
   cl := List( ConjugacyClasses( G ), Representative );
   so := Difference( cl, Centre( G) );
   for x in so do
      N := NormalClosure( G, Group( x ) );
      if IsAbelian(N) then
         R := RestrictedClassFunction( chi , N );
         C := ConstituentsOfCharacter( R );
         I := InertiaSubgroup( G, C[1] );
         Ri := RestrictedClassFunction( chi, I );
         Ci := ConstituentsOfCharacter( Ri );
         for t in Ci do
           Rc := RestrictedClassFunction( t, N );
           Cc := ConstituentsOfCharacter( Rc );
           if C[1] in Cc then
              Info( InfoWarning, 1, " the given character of degree ", chi[1], " is imprimitive" );
              return InducedSubgroupRepresentation( G, IrreducibleAffordingRepresentation( t ));
           fi;
         od;
      fi;
   od;
   return PerfectAbelianSocleRepresentation( G, chi );
 else
   S := Difference( n, a );
 ## S is the non-abelian part of the soc(G/Z) which is the direct product of simple groups
 ## if Length( S ) = 1 then G is a central cover of a simple group and we use Dixon's method
   if Length( arg ) = 3 and Length( S ) = 1 and arg[3] = 1 then
       return true;
   fi;
 ## if Length( S ) >= 5 then we apply Gallaghe's method.
 ## For checking the components of S we only use the size of these components.
   ls := List(S,Size);
   A := 0;
   Sort(ls);
   for k in [1..Length(ls)] do
      if ls[k] > A then
         A:=ls[k];
         j := 0;
      fi;
      j:=j+1;
      if j > 4 then
         return  SocleMoreComponents( G, chi, S );
      fi;
   od;
 ## if Length( S ) > 1 then, if number of isomorphic components of S is <= 5 then we use
 ## the Tensor Product method.
   for s in S do
      Append( M, GeneratorsOfGroup( s ) );
   od;
   c := Centralizer( Image( hom ), Group( M ) );
   if Size(c) > 1 then
      Add( S, c );
   fi;
   S := List( S, i -> PreImage( hom, i ) );
   Co := [ One( Source( hom ) ) ];
   for i in [ 1..Length( S ) ] do
      R := RestrictedClassFunction( chi, S[i] );
      C := ConstituentsOfCharacter( R );
      r := CharacterSubgroupRepresentation( C[1] );
      Add( L, GeneratorsOfGroup( Image( r ) ) );
      Co := ListX(Co, GeneratorsOfGroup( S[i] ), \*);
   od;
   if Length(L) = 1 then
      return r;
   fi;
   r := L[1];
   for i in [ 2..Length( L ) ] do
     r := ListX( r, L[i], KroneckerProduct );
   od;
 fi;
 if arg[3] <> 1 then
 M := arg[3];
 for i in [1..Length(M)] do
    if Size(G) = M[i][1]  and Size(Centre(G)) = M[i][2] and chi[1] = M[i][3] and IsSimple( Image(hom) ) then
       if Length[ M[i] ] = 5 then
          run := M[i][4]( G, chi, M[i][5], hom );
       fi;
       if Length[ M[i] ] = 4 then
          run := M[i][4]( chi, hom );
       fi;
       return run;
    fi;
 od;
 fi;
 return GroupHomomorphismByImagesNC( G, Group( r ), Co, r );
 end );
