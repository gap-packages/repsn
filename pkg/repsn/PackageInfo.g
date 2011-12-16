#############################################################################
#############################################################################
##              PackageInfo.g for the GAP 4 package Repsn                  ##
##                                                                         ##
##                          Vahid Dabbaghian                               ##

SetPackageInfo( rec(
 PackageName    := "Repsn",
 MyVersion      := "3r0p2",
 MyWWWHome      := "http://www.sfu.ca/~vdabbagh",
 Subtitle       := "A GAP4 Package for constructing representations of finite groups",
 Version        := JoinStringsWithSeparator( SplitString( ~.MyVersion, "rp" ), "." ),
 Autoload       := false,
 Date           := "25/08/2011",
 PackageWWWHome := "http://www.sfu.ca/~vdabbagh/gap/repsn.html",
 ArchiveURL     := Concatenation( ~.MyWWWHome, "/gap/", LowercaseString( ~.PackageName ),~.MyVersion ),
 ArchiveFormats := "-win.zip",


Persons := [rec(
 LastName       := "Dabbaghian",
 FirstNames     := "Vahid",
 IsAuthor       := true,
 IsMaintainer   := true,
 Email          := "vdabbagh@sfu.ca",
 WWWHome        := ~.MyWWWHome,
 Place          := "Burnaby, Canada",
 Institution    := "Department of Mathematics, Simon Fraser University",
 PostalAddress  := Concatenation( [
       "Vahid Dabbaghian\n",
       "Department of Mathematics\n",
       "Simon Fraser University\n",
       "Burnaby, British Columbia\n",
       "V5A 1S6 Canada"] )  ), ],


 Status         := "accepted",
 CommunicatedBy := "Charles Wright (Eugene)",
 AcceptDate     := "05/2004",
 README_URL     := "http://www.sfu.ca/~vdabbagh/gap/README.repsn",
 PackageInfoURL := "http://www.sfu.ca/~vdabbagh/gap/PackageInfo.g",
 AbstractHTML   := Concatenation( [
                "The package provides <span class=\"pkgname\">GAP</span> functions ",
                "for computing characteristic zero matrix representations of finite groups."] ),

PackageDoc := rec(
 BookName         := "Repsn",
 ArchiveURLSubset := [ "doc" ],
 HTMLStart        := "doc/chap0.html",
 PDFFile          := "doc/manual.pdf",
 SixFile          := "doc/manual.six",
 LongTitle        := "Constructing matrix representations of finite groups",
 Autoload         := true ),

Dependencies := rec(
 GAP                    := ">= 4.3",
 NeededOtherPackages    := [ ],
 SuggestedOtherPackages := [ ],
 ExternalConditions     := [ ] ),

 AvailabilityTest := ReturnTrue,
# TestFile         := "", # "tst/testall.g",
 Keywords         := ["group representations", "matrix representations", "Dixon's method"],

BannerString := Concatenation(
   "-------------------------------------------------------\n",
   "Repsn for Constructing Representations of Finite Groups\n",
   "                  Version 3.0.2                        \n",
   "                                                       \n",
   "                    Written by                         \n",
   "                 Vahid Dabbaghian                      \n",
   "-------------------------------------------------------\n" ),
) );

