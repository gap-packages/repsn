#######################################################################################################
#######################################################################################################
## The covering groups $6.A_6$, $2.A_7$, $3.A_7$, $6.A_7$ and $2.A_8$ have some irreducible characters
## $\chi$ such that there is a $\chi$-subgroup, say $H$ , such that $H$ is not a $p$-group. 
## Also the groups $3.O_7(3)$, $3.U_6(2)$ and the covering groups of $U_4(3)$ have some irreducible 
## characters such that their restriction to a maximal subgroup is irreducible.
##                                                                                                   
##                          |G|     ,|Z(G)|,chi[1],               functions                 , |H|    
##                     ---------------------------------------------------------------------------   
 SpecialCoversData := [[6*360       ,  6   ,  12  ,  CoversOfAlternatingGroupsRepresentation, 10 ],
                       [2*2520      ,  2   ,  20  ,  CoversOfAlternatingGroupsRepresentation, 20 ],
                       [3*2520      ,  3   ,  21  ,  CoversOfAlternatingGroupsRepresentation, 12 ],
                       [3*2520      ,  3   ,  24  ,  CoversOfAlternatingGroupsRepresentation, 20 ],
                       [6*2520      ,  6   ,  20  ,  CoversOfAlternatingGroupsRepresentation, 20 ],
                       [6*2520      ,  6   ,  21  ,  CoversOfAlternatingGroupsRepresentation, 12 ],
                       [6*2520      ,  6   ,  24  ,  CoversOfAlternatingGroupsRepresentation, 20 ],
                       [2*20160     ,  2   ,  24  ,  CoversOfAlternatingGroupsRepresentation, 15 ],
                       [2*3265920   ,  2   ,  21  ,  CoversOfU43Representation],
                       [2*3265920   ,  2   ,  20  ,  CoversOfU43Representation],
                       [3*3265920   ,  3   ,  21  ,  CoversOfU43Representation],
                       [3*3265920   ,  3   ,  15  ,  CoversOfU43Representation],
                       [4*3265920   ,  4   ,  21  ,  CoversOfU43Representation],
                       [4*3265920   ,  4   ,  20  ,  CoversOfU43Representation],
                       [6*3265920   ,  6   ,  21  ,  CoversOfU43Representation],
                       [6*3265920   ,  6   ,  6   ,  CoversOfU43Representation],
                       [3*4585351680,  3   ,  27  ,  CoversOfO73Representation],
                       [3*9196830720,  3   ,  22  ,  CoversOfU62Representation],
                       [3*9196830720,  3   ,  21  ,  CoversOfU62Representation]];;



################################################################################################
################################################################################################
## Suppose $G$ is a group with abelian $Soc(G/Z(G))$ and $F$ is the Fitting subgroup of $G$.
## If $\chi$ is an irreducible character of $G$ of degree 16, 20, 24, 28 and 30 such that
## $\chi_F$ has a single constituent with multiplicity 4, 4, 6, 4, and 6, respectively, then
## there exists a $\chi$-subgroup $H$ containing the centre of $G$.
##                                                                                            
##                          chi[1], m , |G/F|,|H/Z(G)|,          functions                    
##                         --------------------------------------------------------------------
 SpecialAbelianSocleData:=[[  16  , 4 ,  60  ,   12   ,RepresentationSpecialCharactersSubgroup],
                           [  16  , 4 ,  60  ,   16   ,RepresentationSpecialCharactersSubgroup],
                           [  16  , 4 ,  60  ,   24   ,RepresentationSpecialCharactersSubgroup],
                           [  16  , 4 ,  360 ,   36   ,RepresentationSpecialCharactersMaximal ],
                           [  20  , 4 ,  120 ,   10   ,RepresentationSpecialCharactersSubgroup],
                           [  24  , 6 ,  60  ,   24   ,RepresentationSpecialCharactersSubgroup],
                           [  24  , 6 ,  360 ,   12   ,RepresentationSpecialCharactersSubgroup],
                           [  24  , 6 ,  360 ,   16   ,RepresentationSpecialCharactersSubgroup],
                           [  28  , 4 ,  336 ,   42   ,RepresentationSpecialCharactersMaximal ],
                           [  28  , 4 ,  336 ,   48   ,RepresentationSpecialCharactersMaximal ],
                           [  30  , 6 ,  120 ,   25   ,RepresentationSpecialCharactersSubgroup]];;




