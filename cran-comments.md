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
0 errors | 0 warnings | 0 notes

## Downstream dependencies
* rlang, version 0.4.0 successfully installed.