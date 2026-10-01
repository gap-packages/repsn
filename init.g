#################################################################
##          init.g file for the GAP 4 package - Repsn
##
##                      Vahid Dabbaghian
##

## the NC variants of the PreImages operations exist from GAP 4.17 on
if not IsBound( PreImagesNC ) then
    BindGlobal( "PreImagesNC", PreImages );
fi;
if not IsBound( PreImagesElmNC ) then
    BindGlobal( "PreImagesElmNC", PreImagesElm );
fi;
if not IsBound( PreImagesSetNC ) then
    BindGlobal( "PreImagesSetNC", PreImagesSet );
fi;
if not IsBound( PreImagesRepresentativeNC ) then
    BindGlobal( "PreImagesRepresentativeNC", PreImagesRepresentative );
fi;

## read the actual code.
ReadPackage( "repsn", "gap/func.g" );
ReadPackage( "repsn", "gap/data.g" );
ReadPackage( "repsn", "gap/repsn.g" );
