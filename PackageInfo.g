#############################################################################
#############################################################################
##              PackageInfo.g for the GAP 4 package Repsn                  ##
##                                                                         ##
##                      Vahid Dabbaghian-Abdoly                            ##

SetPackageInfo( rec(
 PackageName    := "Repsn",
 MyVersion      := "1r0p0",
 MyWWWHome      := "http://www.math.carleton.ca/~vdabbagh",
 Subtitle       := "A GAP4 Package for constructing representations of finite groups",
 Version        := JoinStringsWithSeparator( SplitString( ~.MyVersion, "rp" ), "." ),
 Autoload       := false,
 Date           := "25/04/2004",
 PackageWWWHome := "http://www.math.carleton.ca/~vdabbagh/gap/repsn.html",
 ArchiveURL     := Concatenation( ~.PackageWWWHome, "/", LowercaseString( ~.PackageName ),~.MyVersion ),
 ArchiveFormats := ".zip, .tar.gz, .tar.bz2, .zoo",


Persons := [rec(
 LastName       := "Dabbaghian-Abdoly",
 FirstNames     := "Vahid",
 IsAuthor       := true,
 IsMaintainer   := true,
 Email          := "vdabbagh@math.carleton.ca",
 WWWHome        := ~.MyWWWHome,
 Place          := "Ottawa, Canada",
 Institution    := "School of Mathematics and Statistics, Carleton University",
 PostalAddress  := Concatenation( [
       "Vahid Dabbaghian-Abdoly\n",
       "School of Mathematics and Statistics\n",
       "Carleton University\n",
       "Ottawa, Ontario\n",
       "K1S 5B6 Canada"] )  ), ],


 Status         := "accepted",
 CommunicatedBy := "Charles Wright",
 AcceptDate     := "05/2004",
 README_URL     := "http://www.math.carleton.ca/~vdabbagh/gap/README.repsn",
 PackageInfoURL := "http://www.math.carleton.ca/~vdabbagh/gap/PackageInfo.g",
 AbstractHTML   := Concatenation( [
                "The package provides <span class=\"pkgname\">GAP</span> functions ",
                "for computing characteristic zero matrix representations of finite groups."] ),

PackageDoc := rec(
 BookName         := "Repsn",
 ArchiveURLSubset := [ "doc", "html" ],
 HTMLStart        := "html/chap0.html",
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
 TestFile         := "", # "tst/testall.g",
 Keywords         := ["group representations", "matrix representations", "Dixon's method"],

BannerString := Concatenation(
   "---------------------------------------------------------------\n",
   "Loading repsn for Constructing Representations of Finite Groups\n",
   "                                                               \n",
   "            Written by Vahid Dabbaghian-Abdoly                 \n",
   "---------------------------------------------------------------\n" ),
) );

