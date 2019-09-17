## Test environments
* Local: OS X install, R version 3.5.3
* AppVeyor: https://ci.appveyor.com/project/stineb/rpmodel
* using devtools::check_win_devel() returns error:
checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! Package inputenc Error: Unicode char − (U+2212)
(inputenc)                not set up for use with LaTeX.


## R CMD check results
### devtools::check() locally:
0 errors | 0 warnings | 0 notes

### devtools::check_rhub()
0 errors | 0 warnings | 1 note 

This is: 
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Benjamin Stocker <benjamin.stocker@gmail.com>'
New submission

Non-FOSS package license (file LICENSE)

Possibly mis-spelled words in DESCRIPTION:
  Stocker (7:114)
  al (7:84, 7:103, 7:125)
  et (7:81, 7:100, 7:122)

The Title field should be in title case. Current version is:
'P-model'
In title case that is:
'P-Model'

The Description field should not start with the package name,
  'This package' or similar.

Comment from B. Stocker: That's all ok as is.

## Downstream dependencies
* rlang, version 0.4.0 successfully installed.