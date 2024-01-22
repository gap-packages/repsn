#############################################################################
#############################################################################
##              PackageInfo.g for the GAP 4 package Repsn                  ##
##                                                                         ##
##                          Vahid Dabbaghian                               ##

SetPackageInfo( rec(

  PackageName    := "Repsn",
  Subtitle       := "Constructing representations of finite groups",
  Version        := "3.1.2",
  Date           := "22/01/2024", # dd/mm/yyyy format
  License        := "GPL-2.0-or-later",

 Persons := [
   rec(
     LastName       := "Dabbaghian",
     FirstNames     := "Vahid",
     IsAuthor       := true,
     IsMaintainer   := false,
     Email          := "vdabbagh@sfu.ca",
     WWWHome        := "https://www.sfu.ca/~vdabbagh",
     Place          := "Burnaby, Canada",
     Institution    := "Department of Mathematics, Simon Fraser University",
     PostalAddress  := Concatenation( [
           "Vahid Dabbaghian\n",
           "Department of Mathematics\n",
           "Simon Fraser University\n",
           "Burnaby, British Columbia\n",
           "V5A 1S6 Canada"] )
    ), 

    rec(
      LastName      := "GAP Team",
      FirstNames    := "The",
      IsAuthor      := false,
      IsMaintainer  := true,
      Email         := "support@gap-system.org",
    ),
  ],


  Status         := "accepted",
  CommunicatedBy := "Charles Wright (Eugene)",
  AcceptDate     := "05/2004",

  PackageWWWHome  := "https://gap-packages.github.io/repsn/",
  README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
  PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
  SourceRepository := rec(
      Type := "git",
      URL := "https://github.com/gap-packages/repsn",
  ),
  IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
  ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                   "/releases/download/v", ~.Version,
                                   "/repsn-", ~.Version ),
  ArchiveFormats := ".tar.gz",

  AbstractHTML   := Concatenation( [
                 "The package provides <span class=\"pkgname\">GAP</span> functions ",
                 "for computing characteristic zero matrix representations of finite groups."] ),

  PackageDoc := rec(
    BookName         := "Repsn",
    ArchiveURLSubset := [ "doc" ],
    HTMLStart        := "doc/chap0_mj.html",
    PDFFile          := "doc/manual.pdf",
    SixFile          := "doc/manual.six",
    LongTitle        := "Constructing representations of finite groups",
    Autoload         := true
  ),

  Dependencies := rec(
    GAP                    := ">= 4.8",
    NeededOtherPackages    := [ ],
    SuggestedOtherPackages := [ ],
    ExternalConditions     := [ ]
  ),

  AvailabilityTest := ReturnTrue,
  TestFile         := "tst/testall.g",
  Keywords         := ["group representations", "matrix representations", "Dixon's method"],

  AutoDoc := rec(
      TitlePage := rec(
          Acknowledgements := """
<P/>The first version of this package was obtained during my Ph.D.
studies at Carleton University. I would like to express deep gratitude
to my supervisor Professor John D. Dixon whose guidance and support were
crucial for the successful completion of this project. I also thank
Professor Charles Wright and referees for pointing out some important
comments to improve <Package>Repsn</Package>.

<P/>This documentation was prepared with the <Package>GAPDoc</Package>
package by Frank Lübeck and Max Neunhöffer.
""",
      ),
  ),

) );

